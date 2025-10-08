rule parse_ssp:
    """
    Generate all tracks for a given SSP

    Test with:
    snakemake -c1 data/out/genesis-CRH/SSP-585.flag
    """
    input:
        expand(
            "data/out/genesis-{{genesis}}/SSP-{{ssp}}/GCM-{gcm}/epoch-{epoch}/tracks.gpq",
            gcm=GCMS,
            epoch=EPOCHS,
        )
    output:
        "data/out/genesis-{genesis}/SSP-{ssp}.flag"
    shell:
        "touch {output}"


rule parse_genesis_method:
    """
    Generate all tracks for a given genesis method

    Test with:
    snakemake -c1 data/out/genesis-CRH.flag
    """
    input:
        expand(
            "data/out/genesis-{{genesis}}/SSP-{ssp}/GCM-{gcm}/epoch-{epoch}/tracks.gpq",
            ssp=SSPS,
            gcm=GCMS,
            epoch=EPOCHS,
        )
    output:
        "data/out/genesis-{genesis}.flag"
    shell:
        "touch {output}"


rule parse_all:
    """
    Generate all tracks for all genesis methods

    Test with:
    snakemake -c1 data/out.flag
    """
    input:
        expand(
            "data/out/genesis-{genesis}.flag",
            genesis=GENESIS_METHODS,
        )
    output:
        "data/out.flag"
    shell:
        "touch {output}"
