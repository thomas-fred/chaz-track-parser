from math import prod

import numpy as np
import pandas as pd
import xarray as xr

from chaz.models import r_max


GENESIS_METHOD_CODE = {  # Short code for use in storm_id column
    "SD": "S",  # saturation deficit
    "CRH": "H",  # relative humidity
}
KNOTS_PER_METER_PER_SECOND = 0.051444
REFERENCE_DATE = "1950-01-01"  # netCDF 'time' variable is days since this date


def to_table(ds: xr.Dataset, genesis_method: str, sample_id: str) -> pd.DataFrame:
    """
    Create timestamps from reference and offsets
    Convert sparse datacube to dense table
    Infer radius to maximum winds with regression model
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

    # N.B. We transpose the input arrays to have stormID as the first dim, then lifelength
    data = []
    length = prod(ds.time_datetime.data.shape)
    storm_number = np.repeat(ds.stormID.data, ds.lifelength.shape)
    storm_start_year = np.repeat(
        ds.time_datetime.data.T[:, 0].astype('datetime64[Y]').astype(int) + 1970,
        ds.lifelength.shape
    )
    timestamps = ds.time_datetime.data.T.reshape(length)
    longitude = ds.longitude.data.T.reshape(length)
    latitude = ds.latitude.data.T.reshape(length)

    # Extracting the data takes about 10s to generate 20M+ rows
    for i in ds.ensembleNum.data:
        data.append(
            pd.DataFrame(
                {
                    "storm_number": storm_number,
                    "ensemble_number": np.ones(length) * i,
                    "storm_start_year": storm_start_year,
                    "datetime": timestamps,
                    "longitude_deg": longitude,
                    "latitude_deg": latitude,
                    "max_wind_speed_ms": ds.Mwspd.data[i, :, :].T.reshape(length) * KNOTS_PER_METER_PER_SECOND,
                }
            )
        )
    df = pd.concat(data)
    df["ensemble_number"] = df["ensemble_number"].astype(int)
    df = df.dropna().reset_index(drop=True)

    # Labeling with IDs takes about 35s for 17M rows
    df["storm_id"] = \
        f"G{GENESIS_METHOD_CODE[genesis_method]}" \
        + f"S{int(sample_id):03d}" \
        + df["storm_start_year"].map(lambda x: f"Y{x:04d}") \
        + df["storm_number"].map(lambda x: f"N{x:05d}") \
        + df["ensemble_number"].map(lambda x: f"E{x:02d}")

    df = df.loc[
        :,
        [
            "datetime",
            "storm_id",
            "storm_number",
            "ensemble_number",
            "longitude_deg",
            "latitude_deg",
            "max_wind_speed_ms",
        ]
    ]

    # Infer radius to maximum sustained winds with a model
    df["RMW_km"] = r_max(df.max_wind_speed_ms, df.latitude_deg)

    return df
