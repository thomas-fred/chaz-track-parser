from math import prod

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr

from chaz.models import ENV_PRESSURE, r_max_willoughby_2004, p_min_holland_1980


GENESIS_METHOD_CODE = {  # Short code for use in storm_id column
    "SD": "S",  # saturation deficit
    "CRH": "H",  # relative humidity
}
METERS_PER_SECOND_PER_KNOT = 0.51444
REFERENCE_DATE = "1950-01-01"  # netCDF 'time' variable is days since this date


def signed_longitude_to_strictly_positive(coords: tuple[float, float]) -> tuple[float, float]:
    return [(long + 360 if long < 0 else long, lat) for long, lat in coords]


def chaz_to_table(ds: xr.Dataset, genesis_method: str, sample_id: str) -> gpd.GeoDataFrame:
    """
    Create timestamps from reference and offsets
    Convert sparse datacube to dense table
    Create unique track_id
    Create geometry column
    """

    # Check ordering of dimensions
    for var in ("time", "longitude", "latitude"):
        assert ds[var].dims == ("lifelength", "stormID")
    assert ds["Mwspd"].dims == ("ensembleNum", "lifelength", "stormID")

    # Check we won't overflow our fixed-width strings in `storm_id`
    assert np.logical_and(0 <= ds.ensembleNum, ds.ensembleNum < 1E2).all()
    assert np.logical_and(0 <= ds.stormID, ds.stormID < 1E5).all()

    # Following from CHAZ_analysis.ipynb, create a timestamp
    ds['time_datetime'] = (
        ('lifelength', 'stormID'),
        pd.to_datetime(
            ds.time.values.ravel(), unit='D', origin=pd.Timestamp(REFERENCE_DATE)
        ).values.reshape(ds.time.values.shape)
    )

    data = []
    length = prod(ds.time_datetime.data.shape)
    storm = np.repeat(ds.stormID.data, ds.lifelength.shape)
    # N.B. We transpose the input arrays to have stormID as the first dim, then lifelength
    source_year = np.repeat(
        ds.time_datetime.data.T[:, 0].astype('datetime64[Y]').astype(int) + 1970,
        ds.lifelength.shape
    )
    timesteps = np.tile(range(len(ds.lifelength)), len(ds.stormID.data))
    timestamps = ds.time_datetime.data.T.reshape(length)
    longitude = ds.longitude.data.T.reshape(length)
    latitude = ds.latitude.data.T.reshape(length)

    # Extracting the data takes about 10s to generate 20M+ rows
    for i in ds.ensembleNum.data:
        data.append(
            pd.DataFrame(
                {
                    "time_utc": timestamps,
                    "source_year": source_year,
                    "storm": storm,
                    "sample": np.ones(length) * int(sample_id),
                    "ensemble": np.ones(length) * i,
                    "timestep": timesteps,
                    "longitude_deg": longitude,
                    "latitude_deg": latitude,
                    "max_wind_speed_ms": ds.Mwspd.data[i, :, :].T.reshape(length) * METERS_PER_SECOND_PER_KNOT,
                }
            )
        )
    df = pd.concat(data).dropna().set_index("time_utc", drop=True).astype({"sample": int, "ensemble": int})

    # Labeling with IDs takes about 35s for 17M rows
    df["track_id"] = \
        f"{GENESIS_METHOD_CODE[genesis_method]}" \
        + f"_{int(sample_id):03d}" \
        + df["source_year"].map(lambda x: f"_{x:04d}") \
        + df["storm"].map(lambda x: f"_{x:05d}") \
        + df["ensemble"].map(lambda x: f"_{x:02d}")

    return gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude_deg, df.latitude_deg), crs=4326)


def filter_by_year(df: pd.DataFrame, epoch: int, epoch_half_width_years: int) -> pd.DataFrame:
    return df[
        (epoch - epoch_half_width_years < df.source_year)
        & (df.source_year < epoch + epoch_half_width_years)
    ]


@np.vectorize
def saffir_simpson(wind_speed_ms: float):
    """
    Identify the Saffir-Simpson storm category given a wind speed in m/s.

    N.B. These classifications were developed for 1-minute sustained measurements.
    """

    if wind_speed_ms < 0:
        raise ValueError(f"{wind_speed_ms=} should be positive-valued")
    elif np.isnan(wind_speed_ms):
        return np.nan
    elif wind_speed_ms < 18:
        return -1
    elif wind_speed_ms < 33:
        return 0  # Tropical Storm
    elif wind_speed_ms < 43:
        return 1  # Category 1
    elif wind_speed_ms < 50:
        return 2  # Category 2
    elif wind_speed_ms < 58:
        return 3  # Category 3
    elif wind_speed_ms < 70:
        return 4  # Category 4
    elif wind_speed_ms >= 70:
        return 5  # Category 5


def tag_category(df: pd.DataFrame) -> pd.DataFrame:
    df["ss_category"] = saffir_simpson(df["max_wind_speed_ms"])
    return df


def tag_basin(df: gpd.GeoDataFrame, basins: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Tag each track point with the encompassing TC basin."""
    return df.to_crs(epsg=4326).sjoin(basins.to_crs(epsg=4326)).drop(columns="index_right")


def estimate_rmw(df: pd.DataFrame) -> pd.DataFrame:
    """Infer radius to maximum sustained winds with a model fit."""
    df["radius_to_max_winds_km"] = r_max_willoughby_2004(df.max_wind_speed_ms, df.latitude_deg)
    return df


def estimate_p_min(df: pd.DataFrame) -> pd.DataFrame:
    """Infer minimum eye pressure from a model fit."""
    df["min_pressure_hpa"] = p_min_holland_1980(
        df.basin_id.map(ENV_PRESSURE),
        df.max_wind_speed_ms,
        df.radius_to_max_winds_km * 1_000,
        df.geometry.y
    )
    # See Emanuel 1986 for a discussion of minimum possible pressures (circa 850 hPa)
    # https://doi.org/10.1175/1520-0469(1986)043%3C0585:AASITF%3E2.0.CO;2
    # Climate change may reduce this slighty, but very unlikely to pass 800 hPa
    df["min_pressure_hpa"] = np.where(
        (df["min_pressure_hpa"] < 800) | (df["min_pressure_hpa"] > 1050),
        np.nan,
        df["min_pressure_hpa"]
    )
    return df


def normalise_frequency(freq: pd.DataFrame, tracks: pd.DataFrame) -> pd.DataFrame:
    """
    Relabel years of tracks so that typical TC frequency (per-basin) is
    respected. Return a number of tracks and duration that is the largest
    possible while respecting the supplied long-term annual frequency
    constrant.

    Args:
        freq: Table of ´tc_per_year´, with ´basin_id´ index.
        tracks: Table of pre-processed track/eyeattributes, with DatetimeIndex.

    Returns:
        Mutated tracks table with normalised annual frequencies.
    """
    freq["tc_count"] = tracks.loc[:, ["track_id", "basin_id"]].groupby("basin_id").nunique()
    # Given the desired average TC frequency, how many years can we represent?
    freq["duration_year"] = np.round(freq["tc_count"] / freq["tc_per_year"], 0).astype(int)

    track_basin = tracks.loc[:, ["track_id", "basin_id"]].drop_duplicates().set_index(["track_id"], drop=True)
    track_year = []
    for basin_id in freq.index:
        basin = track_basin[track_basin.basin_id==basin_id].copy()
        basin["year"] = np.round(np.random.rand(len(basin)) * freq.loc[basin_id, "duration_year"], 0).astype(int)
        track_year.append(basin)
    track_year = pd.concat(track_year).sort_values("year")
    track_year = track_year[track_year.year < freq.duration_year.min()]
    tracks = tracks.join(track_year.year, on="track_id", how="inner")

    # Label with a tc_number (0 is first TC of the year)
    # Unique within a given year of SSP-GCM-genesis-method combination
    tracks["tc_number"] = -1
    for year in tracks.year.unique():
        mask = tracks.year == year
        tc_number, track_id = pd.factorize(tracks.loc[mask, "track_id"])
        tracks.loc[mask, "tc_number"] = tc_number
    assert -1 not in tracks["tc_number"].unique()

    tracks = tracks.sort_values(["year", "tc_number", "timestep"])
    return tracks.loc[
        :,
        [
            "year",
            "tc_number",
            "timestep",
            "track_id",
            "source_year",
            "sample",
            "ensemble",
            "basin_id",
            "ss_category",
            "max_wind_speed_ms",
            "radius_to_max_winds_km",
            "min_pressure_hpa",
            "geometry",
        ]
    ]

