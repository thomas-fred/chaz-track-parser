rule parse_ssp:
    """
    Generate all tracks for a given SSP

    Test with:
    snakemake -c1 data/out/genesis-CRH/SSP-585.flag
    """
    input:
        expand(
            "{{data}}/out/genesis-{{genesis}}/SSP-{{ssp}}/GCM-{gcm}/epoch-{epoch}/norm-freq",
            gcm=GCMS,
            epoch=EPOCHS,
        )
    output:
        "{data}/out/genesis-{genesis}/SSP-{ssp}.flag"
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
            "{{data}}/out/genesis-{{genesis}}/SSP-{ssp}/GCM-{gcm}/epoch-{epoch}/norm-freq",
            ssp=SSPS,
            gcm=GCMS,
            epoch=EPOCHS,
        )
    output:
        "{data}/out/genesis-{genesis}.flag"
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
            "{{data}}/out/genesis-{genesis}/SSP-{ssp}/GCM-{gcm}/epoch-{epoch}/norm-freq",
            genesis=GENESIS_METHODS,
            ssp=SSPS,
            gcm=GCMS,
            epoch=EPOCHS,
        )
    output:
        flag = "{data}/out.flag"
    shell:
        "touch {output}"


checkpoint write_manifest:
    """
    Trigger upstream work and write all output geoparquet file paths to a JSON
    file.

    N.B. We don't know ahead of time how many complete millennia of TC tracks
    each genesis-SSP-GCM-epoch will contain, (how many samples there will be
    for each), hence the checkpoint.
    """
    input:
        rules.parse_all.output.flag
    output:
        manifest = "{data}/out/manifest.json"
    run:
        import glob
        import json
        import os

        paths = sorted(
            glob.glob(
                f"{wildcards.data}/out/genesis-*/SSP-*/GCM-*/epoch-*/norm-freq/sample-*/tracks.pq",
                recursive=True
            )
        )
        paths = [path.replace(".pq", ".gpq") for path in paths]
        os.makedirs(os.path.dirname(output.manifest), exist_ok=True)
        with open(output.manifest, "w") as fp:
            json.dump(paths, fp, indent=2)


def read_manifest(wildcards):
    import json

    check = checkpoints.write_manifest.get(data=wildcards.data)
    with open(check.output.manifest) as fp:
        return json.load(fp)

rule archive:
    """
    Collate all of the normalised epochal results.

    Test with:
    snakemake -c1 data/out/CHAZ-normalised-freq/
    """
    input:
        manifest = read_manifest
    threads: 8
    output:
        archive_dir = directory("{data}/out/CHAZ-normalised-freq/")
    shell:
        """
        mkdir -p {output}

        for INPUT_FILEPATH in {input.manifest}; do

            OUTPUT_FILENAME=CHAZ_$(echo $INPUT_FILEPATH | grep -o 'genesis-.*' | sed 's#/#_#g' | sed 's/_norm-freq//g' | sed 's/_tracks//g')

            cp $INPUT_FILEPATH {output.archive_dir}/$OUTPUT_FILENAME &

            # Semaphore to run a maximum of `threads` concurrent copies
            while [ "$(jobs -rp | wc -l)" -ge "{threads}" ]; do
                sleep 0.1
            done

        done

        wait
        """

