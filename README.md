Parse CHAZ tropical cyclone tracks from netCDF to parquet format

## Process

- Read netCDF
- Create timestamps
- Convert to long format
- Estimate radius to maximum sustained winds (Willoughby 2004 log-linear model)
- Estimate minimum pressure (Holland 1980 profile with Vickery & Waldera 2008 fit for shape parameter)
- Save as parquet format
- Concatenate over samples

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
