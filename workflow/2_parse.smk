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
    snakemake -c1 data/out/genesis-CRH/SSP-585/GCM-UKESM1-0-LL/epoch-2010/norm-freq
    """
    input:
        historic_frequency = rules.historic_frequency.output.frequency,
        baseline_epoch = "{data}/out/genesis-{genesis}/SSP-{ssp}/GCM-{gcm}/epoch-2010/tracks-raw-freq.pq",
        target_epoch = rules.concatenate_samples.output.epoch,
    resources:
        mem_mb = 24000
    output:
        epoch = directory("{data}/out/genesis-{genesis}/SSP-{ssp}/GCM-{gcm}/epoch-{epoch}/norm-freq/")
    run:
        import logging
        from pathlib import Path

        import geopandas as gpd
        import pandas as pd
        from tqdm import tqdm

        from chaz.parse import normalise_frequency, relative_tc_frequency

        logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

        logging.info("Loading historic per-basin frequencies")
        target_freq = pd.read_csv(input.historic_frequency, na_filter=False).set_index("basin_id")
        synth_baseline_tracks = pd.read_parquet(
            input.baseline_epoch,
            columns=["basin_id", "track_id", "source_year"]
        )

        logging.info("Reading source tracks")
        synth_target_tracks = pd.read_parquet(input.target_epoch)
        target_freq["tc_per_year"] = target_freq["tc_per_year"] * relative_tc_frequency(
            synth_baseline_tracks,
            synth_target_tracks
        )

        logging.info("Normalising track frequencies")
        output_tracks: pd.DataFrame = normalise_frequency(target_freq, synth_target_tracks)

        logging.info("Writing out normalised frequency tracks")
        samples = sorted(output_tracks["sample"].unique())
        # Skip last sample, it will not be a complete millennia
        for sample in tqdm(samples[:-1], desc="Sample"):
            sample_dir = Path(output.epoch) / f"sample-{sample:03d}"
            sample_dir.mkdir(parents=True)
            output_tracks[output_tracks["sample"] == sample].to_parquet(sample_dir / "tracks.pq")


rule parquet_to_geoparquet:
    """
    Convert parquet with float longitude and latitude columns to geoparquet with
    Point geometries.

    Test with:
    snakemake -c1 data/out/genesis-CRH/SSP-585/GCM-UKESM1-0-LL/epoch-2010/norm-freq/sample-000/tracks.gpq
    """
    input:
        parquet = "{path}.pq"
    output:
        geoparquet = "{path}.gpq"
    run:
        import geopandas as gpd
        import pandas as pd

        df = pd.read_parquet(input.parquet)
        gdf = gpd.GeoDataFrame(
            df,
            geometry=gpd.points_from_xy(df.longitude_deg, df.latitude_deg),
            crs="EPSG:4326",
        )
        gdf.to_parquet(output.geoparquet)


rule diagnostic_plot:
    """
    Draw diagnostic plots and maps of meteorological variables and TC frequency.

    Test with:
    snakemake -c1 data/out/genesis-CRH/SSP-585/GCM-UKESM1-0-LL/epoch-2010/norm-freq/sample-000/plots
    """
    input:
        tracks = "{data}/out/genesis-{genesis}/SSP-{ssp}/GCM-{gcm}/epoch-{epoch}/norm-freq/sample-{sample}/tracks.gpq"
    output:
        plot_dir = directory("{data}/out/genesis-{genesis}/SSP-{ssp}/GCM-{gcm}/epoch-{epoch}/norm-freq/sample-{sample}/plots")
    run:
        from pathlib import Path

        import geopandas as gpd

        from chaz.plot import plot_joint_distributions, plot_tc_frequency, plot_tc_frequency_per_basin, plot_scatter_map

        plot_dir = Path(output.plot_dir)
        # While Snakemake will create required parent directories for output files
        # It will not create a directory if it is marked as an output itself
        plot_dir.mkdir(parents=True)

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
