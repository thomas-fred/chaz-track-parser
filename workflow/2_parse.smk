rule netcdf_to_geoparquet:
    """
    Read netCDF
    Create timestamps from reference and offsets
    Convert sparse datacube to dense table
    Create unique track_id
    Create geometry column
    Infer radius to maximum winds with regression model
    Infer minimum pressure with a parametric model
    Write as geoparquet

    ~12GB RAM per CPU

    Test with:
    snakemake -c1 data/out/genesis-CRH/SSP-585/GCM-UKESM1-0-LL/sample-000/tracks.gpq
    """
    input:
        basins = rules.generate_basin_definition.output.basins,
        netcdf = f"{config['raw_data_dir']}/ssp{{ssp}}/{{gcm}}_Global_2100_2ens{{sample}}_{{genesis}}_compressed.nc"
    output:
        parquet = "{data}/out/genesis-{genesis}/SSP-{ssp}/GCM-{gcm}/sample-{sample}/tracks.gpq"
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


rule filter_to_epoch:
    """
    Extract window of years that comprise an epoch

    Test with:
    snakemake -c1 data/out/genesis-CRH/SSP-585/GCM-UKESM1-0-LL/epoch-2000/tracks-raw-freq.gpq
    """
    input:
        samples = expand(
            "{{data}}/out/genesis-{{genesis}}/SSP-{{ssp}}/GCM-{{gcm}}/sample-{sample}/tracks.gpq",
            sample=SAMPLES
        ),
    output:
        epoch = "{data}/out/genesis-{genesis}/SSP-{ssp}/GCM-{gcm}/epoch-{epoch}/tracks-raw-freq.gpq"
    run:
        import geopandas as gpd
        import pandas as pd

        from chaz.parse import filter_by_year

        pd.concat(
            [
                filter_by_year(
                    gpd.read_parquet(path),
                    int(wildcards.epoch),
                    int(config["epoch_half_width_years"])
                ) for path in input.samples
            ]
        ).to_parquet(output.epoch)


rule filter_to_epoch:
    """
    Extract window of years that comprise an epoch

    Test with:
    snakemake -c1 data/out/genesis-CRH/SSP-585/GCM-UKESM1-0-LL/epoch-2000/tracks-raw-freq.gpq
    """
    input:
        all_years = rules.concat_samples.output.concat
    output:
        epoch = "{data}/out/genesis-{genesis}/SSP-{ssp}/GCM-{gcm}/epoch-{epoch}/tracks-raw-freq.gpq"
    run:
        import geopandas as gpd

        from chaz.parse import filter_by_year

        filter_by_year(
            gpd.read_parquet(input.all_years),
            int(wildcards.epoch),
            int(config["epoch_half_width_years"])
        ).to_parquet(output.epoch)


rule normalise_frequency:
    """
    Relabel TC years per-basin so annual frequencies are plausible.

    Test with:
    snakemake -c1 data/out/genesis-CRH/SSP-585/GCM-UKESM1-0-LL/epoch-2000/tracks.gpq
    """
    input:
        epoch = rules.filter_to_epoch.output.epoch,
        historic_frequency = rules.historic_frequency.output.frequency,
    output:
        epoch = "{data}/out/genesis-{genesis}/SSP-{ssp}/GCM-{gcm}/epoch-{epoch}/tracks.gpq"
    run:
        import geopandas as gpd
        import pandas as pd

        from chaz.parse import normalise_frequency

        normalise_frequency(
            pd.read_csv(input.historic_frequency, na_filter=False).set_index("basin_id"),
            gpd.read_parquet(input.epoch)
        ).to_parquet(output.epoch)

