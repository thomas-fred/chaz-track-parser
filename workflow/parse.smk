rule netcdf_to_parquet:
    """
    Read netCDF
    Create timestamps from reference and offsets
    Convert sparse datacube to dense table
    Filter to epoch
    Infer radius to maximum winds with regression model
    Write as parquet

    Test with:
    snakemake -c1 data/out/ssp585/UKESM1-0-LL/2000/SD/000/tracks.pq
    """
    input:
        netcdf = f"{config['raw_data_dir']}/{{ssp}}/{{gcm}}_Global_2100_2ens{{sample}}_{{genesis}}_compressed.nc"
    output:
        parquet = "{data}/out/{ssp}/{gcm}/{epoch}/{genesis}/{sample}/tracks.pq"
    run:
        import xarray as xr

        import chaz.parse

        chaz.parse.to_table(
            xr.open_dataset(input.netcdf, engine='netcdf4'),
            wildcards.genesis,
            wildcards.sample,
            int(wildcards.epoch),
            int(config["epoch_half_width_years"]),
        ).to_parquet(output.parquet)


rule concat_samples:
    """
    Concatenate the samples for a given SSP, GCM and genesis method

    Test with:
    snakemake -c1 data/out/ssp585/UKESM1-0-LL/SD/tracks.pq
    """
    input:
        samples = expand("{{data}}/out/{{ssp}}/{{gcm}}/{{epoch}}/{{genesis}}/{sample}/tracks.pq", sample=SAMPLES)
    output:
        concat = "{data}/out/{ssp}/{gcm}/{epoch}/{genesis}/tracks.pq"
    run:
        import pandas as pd

        pd.concat([pd.read_parquet(path) for path in input.samples]).to_parquet(output.concat)


rule parse_ssp:
    input:
        expand(
            "data/out/{{ssp}}/{gcm}/{epoch}/{genesis}/tracks.pq",
            gcm=GCMS,
            epoch=EPOCHS,
            genesis=GENESIS_METHODS,
        )
    output:
        "data/out/{ssp}.flag"
    shell:
        "touch {output}"


rule parse_all:
    input:
        expand(
            "data/out/{ssp}.flag",
            ssp=SSPS,
        )
    output:
        "data/out.flag"
    shell:
        "touch {output}"