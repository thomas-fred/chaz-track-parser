import pandas as pd
from tqdm import tqdm
import xarray as xr

from chaz.constants import KNOTS_PER_METER_PER_SECOND
from chaz.models import r_max


def to_table(ds: xr.Dataset, genesis: str, sample: str) -> pd.DataFrame:
    """
    Create timestamps from reference and offsets
    Convert sparse datacube to dense table
    Infer radius to maximum winds with regression model
    """

    # Following from CHAZ_analysis.ipynb
    # Time variable seems to be days since January 1st 1950
    # Negative offsets -- prior to 1950 -- are for null values
    ds['time_datetime'] = (
        ('lifelength', 'stormID'),
        pd.to_datetime(
            ds.time.values.ravel(), unit='D', origin=pd.Timestamp('1950-01-01')
        ).values.reshape(ds.time.values.shape)
    )

    data = []
    assert len(ds.stormID) < 1E5  # we pad the storm ID to 5 digits
    for storm in tqdm(ds.stormID):
        df = pd.DataFrame(
            {
                "datetime": ds.time_datetime.sel(stormID=storm).to_pandas(),
                "longitude_deg": ds.longitude.sel(stormID=storm).to_pandas(),
                "latitude_deg": ds.latitude.sel(stormID=storm).to_pandas(),
                "max_wind_speed_ms": ds.Mwspd.sel(stormID=storm, ensembleNum=0).to_pandas() * KNOTS_PER_METER_PER_SECOND,
            }
        )
        df = df.dropna(subset="max_wind_speed_ms")
        if not df.empty:
            df["storm_id"] = f"S{genesis}{int(sample):03d}{df.datetime.iloc[0].year:d}{storm:05d}"
            data.append(df)

    data = pd.concat(data).reset_index()
    data = data.loc[:, ["storm_id", "datetime", "longitude_deg", "latitude_deg", "max_wind_speed_ms"]]

    data["RMW_km"] = r_max(data.max_wind_speed_ms, data.latitude_deg)

    return data
