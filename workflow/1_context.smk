rule generate_basin_definition:
    """
    Create tropical cyclone basin polygons.

    Sources:
    - https://www.nature.com/articles/s41597-020-0381-2/tables/2
    - https://www.nature.com/articles/s41597-020-0381-2/figures/1

    Test with:
    snakemake -c1 data/out/basins.geojson
    """
    output:
        basins = "{data}/out/basins.geojson"
    run:
        import geopandas as gpd
        from shapely.geometry.polygon import Polygon

        from chaz.parse import signed_longitude_to_strictly_positive

        delta = 1E-3
        basins = (
            ("EP", ((-180, 60), (-180, 0), (-75, 0), (-75, 10), (-85, 10), (-85, 15), (-90, 15), (-90, 17.5), (-100, 17.5), (-100, 60))),
            ("NA", ((-100, 60), (-100, 17.5), (-90, 17.5), (-90, 15), (-85, 15), (-85, 10), (-75, 10), (-75, 0), (360 - delta, 0), (360 - delta, 60))),
            ("NI", ((30, 60), (30, 0), (100, 0), (100, 60))),
            ("SI", ((10, 0), (10, -60), (135, -60), (135, 0))),
            ("SP", ((135, 0), (135, -60), (-120, -60), (-120, 0))),
            ("WP", ((100, 60), (100, 0), (180, 0), (180, 60))),
        )
        basin_id, geometry = zip(*basins)
        gpd.GeoDataFrame(
            data={
                "basin_id": basin_id,
                "geometry": [Polygon(signed_longitude_to_strictly_positive(poly)) for poly in geometry]
            },
            crs=4326
        ).to_file(output.basins)


rule historic_frequency:
    """
    Identify the frequency of TCs, by basin, from IBTrACS observational data

    IBTrACS input file created with https://github.com/nismod/open-gira/

    Test with -c1 data/out/historic_TC_freq.csv
    """
    input:
        basins = rules.generate_basin_definition.output.basins,
        ibtracs = "{data}/in/IBTrACS.gpq",
    output:
        frequency = "{data}/out/historic_TC_freq.csv"
    run:
        import geopandas as gpd
        import numpy as np

        from chaz.parse import signed_longitude_to_strictly_positive

        basins = gpd.read_file(input.basins)

        # Coerce longitudes to lie in 0 -> 360 range
        ibtracs = gpd.read_parquet(input.ibtracs)
        ibtracs.geometry = gpd.points_from_xy(
            np.where(
                ibtracs.geometry.x < 0,
                ibtracs.geometry.x + 360,
                ibtracs.geometry.x
            ),
            ibtracs.geometry.y
        )

        if "basin_id" in ibtracs.columns:
           ibtracs = ibtracs.drop(columns="basin_id")
        ibtracs = ibtracs.sjoin(basins).drop(columns="index_right")

        # Ignore data outside this window, it shows a markedly different frequency
        ibtracs_start = 2002
        ibtracs_end = 2023
        duration_years: int = ibtracs_end - ibtracs_start + 1

        tc_per_year = ibtracs.loc[
            (ibtracs.year >= ibtracs_start)
            & (ibtracs.year <= ibtracs_end),
            ["basin_id", "track_id"]
        ].groupby(["basin_id"]).nunique().rename(columns={"track_id": "tc_per_year"}) / duration_years

        tc_per_year.round(2).to_csv(output.frequency)

