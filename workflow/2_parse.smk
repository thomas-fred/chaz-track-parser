rule netcdf_to_parquet:
    """
    Read netCDF
    Create timestamps from reference and offsets
    Convert sparse datacube to dense table
    Create unique track_id
    Infer radius to maximum winds with regression model
    Infer minimum pressure with a parametric model
    Write as parquet

    Test with:
    snakemake -c1 data/out/genesis-CRH/SSP-585/GCM-UKESM1-0-LL/sample-000/tracks.pq
    """
    input:
        basins = rules.generate_basin_definition.output.basins,
        netcdf = f"{config['raw_data_dir']}/ssp{{ssp}}/{{gcm}}_Global_2100_2ens{{sample}}_{{genesis}}_compressed.nc"
    output:
        parquet = "{data}/out/genesis-{genesis}/SSP-{ssp}/GCM-{gcm}/sample-{sample}/tracks.pq"
    run:
        import geopandas as gpd
        import xarray as xr

        from chaz.parse import chaz_to_table, tag_category, tag_basin, estimate_rmw, estimate_p_min

        ds = xr.open_dataset(input.netcdf)
        basins = gpd.read_file(input.basins)
        
        first_chunk = True
        for chunk in chaz_to_table(ds, wildcards.genesis, wildcards.sample):
            chunk = tag_category(chunk)
            chunk = tag_basin(chunk, basins)
            chunk = estimate_rmw(chunk)
            chunk = estimate_p_min(chunk)
            
            chunk.to_parquet(output.parquet, engine="fastparquet", append=not first_chunk)
            first_chunk = False


rule filter_to_epoch:
    """
    Extract window of years that comprise an epoch.

    Test with:
    snakemake -c1 data/out/genesis-CRH/SSP-585/GCM-UKESM1-0-LL/sample-000/epoch-2010.pq
    """
    input:
        sample="{data}/out/genesis-{genesis}/SSP-{ssp}/GCM-{gcm}/sample-{sample}/tracks.pq",
    output:
        epoch = "{data}/out/genesis-{genesis}/SSP-{ssp}/GCM-{gcm}/sample-{sample}/epoch-{epoch}.pq"
    run:
        import geopandas as gpd
        import pandas as pd

        from chaz.parse import filter_by_year

        filter_by_year(
            pd.read_parquet(input.sample),
            int(wildcards.epoch),
            int(config["epoch_half_width_years"])
        ).to_parquet(output.epoch)


rule concatenate_samples:
    """
    Concatenate all samples for a given epoch.

    Test with:
    snakemake -c1 data/out/genesis-CRH/SSP-585/GCM-UKESM1-0-LL/epoch-2010/tracks-raw-freq.pq
    """
    input:
        samples=expand(
            "{{data}}/out/genesis-{{genesis}}/SSP-{{ssp}}/GCM-{{gcm}}/sample-{sample}/epoch-{{epoch}}.pq",
            sample=SAMPLES,
        )
    output:
        epoch = "{data}/out/genesis-{genesis}/SSP-{ssp}/GCM-{gcm}/epoch-{epoch}/tracks-raw-freq.pq"
    run:
        import pandas as pd

        first_chunk = True
        for path in input.samples:
            chunk = pd.read_parquet(path)
            chunk.to_parquet(output.epoch, engine="fastparquet", append=not first_chunk)
            first_chunk = False


rule normalise_frequency:
    """
    Anchor per-basin annual TC frequency to observed rates, multiplied by some
    change between synthetic epochs. Relabel years of tracks so that expected TC
    frequencies are respected.

    Test with:
    snakemake -c1 data/out/genesis-CRH/SSP-585/GCM-UKESM1-0-LL/epoch-2010/tracks.pq
    """
    input:
        historic_frequency = rules.historic_frequency.output.frequency,
        baseline_epoch = "{data}/out/genesis-{genesis}/SSP-{ssp}/GCM-{gcm}/epoch-2010/tracks-raw-freq.pq",
        target_epoch = rules.concatenate_samples.output.epoch,
    output:
        epoch = "{data}/out/genesis-{genesis}/SSP-{ssp}/GCM-{gcm}/epoch-{epoch}/tracks.pq"
    run:
        import geopandas as gpd
        import pandas as pd

        from chaz.parse import normalise_frequency

        tracks = normalise_frequency(
            pd.read_csv(input.historic_frequency, na_filter=False).set_index("basin_id"),
            pd.read_parquet(input.baseline_epoch),
            pd.read_parquet(input.target_epoch),
        ).to_parquet(output.epoch)


rule parquet_to_geoparquet:
    """
    Convert a parquet file with `latitude_deg` and `longitude_deg` to a
    geoparquet file with Point geometries.

    Test with:
    snakemake -c1 data/out/genesis-CRH/SSP-585/GCM-UKESM1-0-LL/epoch-2010/tracks.gpq
    """
    input:
        parquet = "{path}.pq"
    output:
        geoparquet = "{path}.gpq"
    run:
        import geopandas as gpd
        import pandas as pd

        tracks = pd.read_parquet(input.parquet)
        gpd.GeoDataFrame(
            tracks,
            geometry=gpd.points_from_xy(tracks.longitude_deg, tracks.latitude_deg),
            crs=4326,
        ).drop(columns=["longitude_deg", "latitude_deg"]).to_parquet(output.geoparquet)


rule diagnostic_plot:
    """
    Draw diagnostic plots and maps of meteorological variables and TC frequency.

    Test with:
    snakemake -c1 data/out/genesis-CRH/SSP-585/GCM-UKESM1-0-LL/epoch-2010/plots
    """
    input:
        tracks = "{data}/out/genesis-{genesis}/SSP-{ssp}/GCM-{gcm}/epoch-{epoch}/tracks.gpq"
    output:
        plot_dir = directory("{data}/out/genesis-{genesis}/SSP-{ssp}/GCM-{gcm}/epoch-{epoch}/plots")
    run:
        from pathlib import Path

        import geopandas as gpd

        from chaz.plot import plot_joint_distributions, plot_tc_frequency, plot_tc_frequency_per_basin, plot_scatter_map

        plot_dir = Path(output.plot_dir)
        # While Snakemake will create required parent directories for output files
        # It will not create a directory if it is marked as an output itself
        plot_dir.mkdir()

        df = gpd.read_parquet(input.tracks)

        title = f"CHAZ: genesis-{wildcards.genesis}, SSP-{wildcards.ssp}, GCM-{wildcards.gcm}, epoch-{wildcards.epoch}"
        basin_ids = sorted(df.basin_id.unique())
        plot_joint_distributions(
            df,
            {
                "max_wind_speed_ms": ((9, 109), "linear"),
                "radius_to_max_winds_km": ((0, 124), "linear"),
                "min_pressure_hpa": ((859, 1049), "linear"),
            },
            title,
            hue_order=basin_ids
        ).savefig(plot_dir / "met_joint_dist.png")
        plot_tc_frequency(df, title).savefig(plot_dir / "tc_frequency.png")
        plot_tc_frequency_per_basin(df, title, hue_order=basin_ids).savefig(plot_dir / "tc_frequency_per_basin.png")
        plot_scatter_map(df, "max_wind_speed_ms", "Max wind speed [ms-1]", title).savefig(plot_dir / "wind_speed_map.png")

