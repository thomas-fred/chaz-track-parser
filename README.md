Parse CHAZ tropical cyclone tracks from netCDF to parquet format and infer some
additional variables.

## Process

- Read netCDF
- Create timestamps
- Convert to long format
- Estimate radius to maximum sustained winds (Willoughby 2004 log-linear model)
- Estimate minimum pressure (Holland 1980 profile with Vickery & Waldera 2008 fit for shape parameter)
- Save as parquet format
- Concatenate over samples
- Calibrate TC frequency to IBTrACS

## Installation

```shell
micromamba create -f environment.yaml -y
```

## Usage

To create tabular records for sample 0, all years of SSP 585 in GCM UKESM1-0-LL,
for the column relative humidity (CRH) genesis method:
```shell
snakemake -c1 data/out/genesis-CRH/SSP-585/GCM-UKESM1-0-LL/sample-000/tracks.gpq
```

And to produce tracks representative of the 2050 epoch, with calibrated frequencies:
```shell
snakemake -c1 data/out/genesis-CRH/SSP-585/GCM-UKESM1-0-LL/epoch-2050/tracks.gpq
```

N.B. For frequency calibration you will need a table of historic IBTrACS
observations at `data/in/IBTrACS.gpq`. The [open-gira](https://github.com/nismod/open-gira)
repository can produce one of these.

See `workflow/*.smk` for more examples.

## Output data schema

Each output file is for a given SSP, genesis method, GCM model and epoch combination.

Each file has the following fields:
```
Name                    Type      Description
year                    int64     (Resampled) year of TC. Use this for return period calculations.
tc_number               int64     Zero indexed counter of TCs within a year.
timestep                int64     Zero indexed counter of timesteps within a track.
track_id                str       Track identifier unique within a genesis/SSP/GCM combination (file).
source_year             int64     Year of GCM this TC was spawned from.
sample                  int64     Sample (0-9) in CHAZ output this track came from.
ensemble                int64     Ensemble member (0-39). Members share a path, but vary in intensity.
basin_id                str       Geographic basin a given track point is within.
ss_category             int64     Saffir-Simpson wind speed category: -1 unclassified, 0 Tropical Storm, 1-5 as usual.
max_wind_speed_ms       float64   Maximum rotational (no advection) wind speed estimated by CHAZ model [ms-1].
radius_to_max_winds_km  float64   Radius from eye at which the maximum is estimated to occur [km]. Post-processed.
min_pressure_hpa        float64   Minimum eye pressure [hPa]. Post-processed.
geometry                geometry  Location of TC eye centre, point. CRS is 4326.
```

Rows are synthetic tropical cyclone observations, each with the above fields.
They come on a datetime index for interpolation purposes. N.B. This index will
not tally with the `year`.
