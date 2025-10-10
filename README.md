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

To create tabular records for SSP 585 and GCM UKESM1-0-LL for the
saturation deficit (SD) genesis method:
```shell
snakemake -c1 data/out/ssp585/UKESM1-0-LL/SD/tracks.gpq
```

See `workflow/*.smk` for more examples.

## Output data schema

Each output file is for a given SSP, genesis method and GCM model combination.

Each file has the following fields:
```
year                         int64
tc_number                    int64
timestep                     int64
track_id                       str
source_year                  int64
sample                       int64
ensemble                     int64
basin_id                       str
max_wind_speed_ms          float64
radius_to_max_winds_km     float64
min_pressure_hpa           float64
geometry                  geometry
```
