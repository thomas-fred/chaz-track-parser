rule generate_basin_definition:
    """
    Create tropical cyclone basin polygons.

    Sources:
    - https://www.nature.com/articles/s41597-020-0381-2/tables/2
    - https://www.nature.com/articles/s41597-020-0381-2/figures/1

    Test with:
    snakemake -c1 data/basins.geojson
    """
    output:
        basins = "{data}/basins.geojson"
    run:
        import geopandas as gpd
        from shapely.geometry.polygon import Polygon

        def signed_longitude_to_strictly_positive(coords: tuple[float, float]) -> tuple[float, float]:
            return [(long + 360 if long < 0 else long, lat) for long, lat in coords]

        delta = 1E-3
        basins = (
            ("EP", ((-180, 60), (-180, 5), (-75, 5), (-75, 10), (-85, 10), (-85, 15), (-90, 15), (-90, 17.5), (-100, 17.5), (-100, 60))),
            ("NA", ((-100, 60), (-100, 17.5), (-90, 17.5), (-90, 15), (-85, 15), (-85, 10), (-75, 10), (-75, 5), (360 - delta, 5), (360 - delta, 60))),
            ("NI", ((30, 60), (30, 5), (100, 5), (100, 60))),
            ("SI", ((10, -5), (10, -60), (135, -60), (135, -5))),
            ("SP", ((135, -5), (135, -60), (-120, -60), (-120, -5))),
            ("WP", ((100, 60), (100, 5), (180, 5), (180, 60))),
        )
        basin_id, geometry = zip(*basins)
        gpd.GeoDataFrame(
            data={
                "basin_id": basin_id,
                "geometry": [Polygon(signed_longitude_to_strictly_positive(poly)) for poly in geometry]
            },
            crs=4326
        ).to_file(output.basins)


rule netcdf_to_geoparquet:
    """
    Read netCDF
    Create timestamps from reference and offsets
    Convert sparse datacube to dense table
    Create unique track_id
    Create geometry column
    Filter to epoch
    Infer radius to maximum winds with regression model
    Write as geoparquet

    ~12GB RAM per CPU

    Test with:
    snakemake -c1 data/out/ssp585/UKESM1-0-LL/2000/SD/000/tracks.gpq
    """
    input:
        basins = rules.generate_basin_definition.output.basins,
        netcdf = f"{config['raw_data_dir']}/{{ssp}}/{{gcm}}_Global_2100_2ens{{sample}}_{{genesis}}_compressed.nc"
    output:
        parquet = "{data}/out/{ssp}/{gcm}/{epoch}/{genesis}/{sample}/tracks.gpq"
    run:
        import geopandas as gpd
        import xarray as xr

        from chaz.parse import chaz_to_table, filter_by_year, tag_basin, append_missing_vars

        append_missing_vars(
            tag_basin(
                filter_by_year(
                    chaz_to_table(
                        xr.open_dataset(input.netcdf, engine='netcdf4'),
                        wildcards.genesis,
                        wildcards.sample
                    ),
                    int(wildcards.epoch),
                    int(config["epoch_half_width_years"]),
                ),
                gpd.read_file(input.basins),
            )
        ).to_parquet(output.parquet)


rule concat_samples:
    """
    Concatenate the samples for a given SSP, GCM and genesis method

    Test with:
    snakemake -c1 data/out/ssp585/UKESM1-0-LL/SD/tracks.gpq
    """
    input:
        samples = expand("{{data}}/out/{{ssp}}/{{gcm}}/{{epoch}}/{{genesis}}/{sample}/tracks.gpq", sample=SAMPLES)
    output:
        concat = "{data}/out/{ssp}/{gcm}/{epoch}/{genesis}/tracks.gpq"
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


rule parse_ssp:
    """
    Generate all tracks for a given SSP

    Test with:
    snakemake -c1 data/out/ssp585.flag
    """
    input:
        expand(
            "data/out/{{ssp}}/{gcm}/{epoch}/{genesis}/tracks.gpq",
            gcm=GCMS,
            epoch=EPOCHS,
            genesis=GENESIS_METHODS,
        )
    output:
        "data/out/{ssp}.flag"
    shell:
        "touch {output}"


rule parse_all:
    """
    Generate all tracks for all SSPs

    Test with:
    snakemake -c1 data/out.flag
    """
    input:
        expand(
            "data/out/{ssp}.flag",
            ssp=SSPS,
        )
    output:
        "data/out.flag"
    shell:
        "touch {output}"