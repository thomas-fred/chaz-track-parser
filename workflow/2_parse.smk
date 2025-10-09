rule netcdf_to_geoparquet:
    """
    Read netCDF
    Create timestamps from reference and offsets
    Convert sparse datacube to dense table
    Create unique track_id
    Create geometry column
    Infer radius to maximum winds with regression model
    Write as geoparquet

    ~12GB RAM per CPU

    Test with:
    snakemake -c1 data/out/genesis-SD/SSP-585/GCM-UKESM1-0-LL/sample-000/tracks.gpq
    """
    input:
        basins = rules.generate_basin_definition.output.basins,
        netcdf = f"{config['raw_data_dir']}/ssp{{ssp}}/{{gcm}}_Global_2100_2ens{{sample}}_{{genesis}}_compressed.nc"
    output:
        parquet = temp("{data}/out/genesis-{genesis}/SSP-{ssp}/GCM-{gcm}/sample-{sample}/tracks.gpq")
    run:
        import geopandas as gpd
        import xarray as xr

        from chaz.parse import chaz_to_table, tag_category, tag_basin, estimate_rmw, estimate_p_min


        estimate_p_min(
            estimate_rmw(
                tag_category(
                    tag_basin(
                        chaz_to_table(
                            xr.open_dataset(input.netcdf, engine="netcdf4"),
                            wildcards.genesis,
                            wildcards.sample
                        ),
                        gpd.read_file(input.basins),
                    )
                )
            )
        ).to_parquet(output.parquet)


rule concat_samples:
    """
    Concatenate the samples for a given SSP, GCM and genesis method

    Test with:
    snakemake -c1 data/out/genesis-SD/SSP-585/GCM-UKESM1-0-LL/tracks.gpq
    """
    input:
        samples = expand(
            "{{data}}/out/genesis-{{genesis}}/SSP-{{ssp}}/GCM-{{gcm}}/sample-{sample}/tracks.gpq",
            sample=SAMPLES
        )
    output:
        concat = "{data}/out/genesis-{genesis}/SSP-{ssp}/GCM-{gcm}/tracks.gpq"
    run:
        import geopandas as gpd
        import pandas as pd

        df = pd.concat([gpd.read_parquet(path) for path in input.samples])

        # Label with a tc_number (0 is first TC of the year)
        # Unique within a given year of SSP-GCM-genesis-method combination
        df["tc_number"] = -1
        for year in df.storm_start_year.unique():
            mask = df.storm_start_year == year
            tc_number, track_id = pd.factorize(df.loc[mask, "track_id"])
            df.loc[mask, "tc_number"] = tc_number
        assert -1 not in df["tc_number"].unique()

        df.loc[
            :,
            [
                "storm_start_year",
                "tc_number",
                "sample",
                "ensemble",
                "timestep",
                "basin_id",
                "track_id",
                "max_wind_speed_ms",
                "radius_to_max_winds_km",
                "geometry",
            ]
        ].to_parquet(output.concat)


rule filter_to_epoch:
    """
    Extract window of years that comprise an epoch

    Test with:
    snakemake -c1 data/out/genesis-SD/SSP-585/GCM-UKESM1-0-LL/epoch-2000/tracks.gpq
    """
    input:
        all_years = rules.concat_samples.output.concat
    output:
        epoch = "{data}/out/genesis-{genesis}/SSP-{ssp}/GCM-{gcm}/epoch-{epoch}/tracks.gpq"
    run:
        import geopandas as gpd

        from chaz.parse import filter_by_year

        filter_by_year(
            gpd.read_parquet(input.all_years),
            int(wildcards.epoch),
            int(config["epoch_half_width_years"])
        ).to_parquet(output.epoch)
