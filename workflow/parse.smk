rule netcdf_to_parquet:
    """
    Test with:
    snakemake -c1 data/output/ssp585/UKESM1-0-LL/000/SD/tracks.pq
    """
    input:
        netcdf = f"{config['raw_data_dir']}/{{SSP}}/{{GCM}}_Global_2100_2ens{{SAMPLE}}_{{GENESIS}}_compressed.nc"
    output:
        # N.B. Baseline and future tracks coexist in these files -- we should separate
        parquet = "{OUTPUT}/{SSP}/{GCM}/{SAMPLE}/{GENESIS}/tracks.pq"
    run:
        import chaz.parse

        chaz.parse.to_table(
            xr.open_dataset(input.netcdf, engine='netcdf4'),
            wildcards.GENESIS,
            wildcards.SAMPLE,
        ).to_parquet(output.parquet)
