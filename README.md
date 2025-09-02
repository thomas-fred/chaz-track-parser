Parse CHAZ tropical cyclone tracks from netCDF to parquet format

## Process

Read netCDF
Create timestamps
Convert to long format
Add radius to maximum sustained winds column with Willoughby 2004 log-linear model
Save as parquet format
Concatenate over samples

## Installation

```shell
micromamba create -f environment.yml -y
```

## Usage

To create tabular records for SSP 585 and GCM UKESM1-0-LL for the
saturation deficit (SD) genesis method:
```shell
snakemake -c1 data/out/ssp585/UKESM1-0-LL/SD/tracks.pq
```

See `workflow/parse.smk` for more examples.
