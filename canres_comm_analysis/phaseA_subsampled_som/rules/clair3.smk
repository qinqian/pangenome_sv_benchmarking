rule clair3:
    threads: 24
    resources:
        mem_mb=96000
    input:
        cram = "output/align/{cell_line}_{platform}/{pair}/{assembly}.cram"
    output:
        outdir = directory("output/clair3/{cell_line}_{platform}/{pair}/{assembly}")
    conda: "clair3"
    shell:
        """
        if [[ {input.cram} =~ "hifi" ]]; then
            model_path=$(ls -d `pwd`/models/hifi_revio/)
            platform="hifi"
        else
            model_path=$(ls -d `pwd`/rerio/clair3_models/r1041_e82_400bps_sup_v430/)
            platform="ont"
        fi
        mkdir -p {output.outdir}
        run_clair3.sh \
        --bam_fn={input.cram} \
        --ref_fn={wildcards.assembly}.fa \
        --threads={threads} \
        --platform=${{platform}} \
        --model_path=${{model_path}} \
        --output={output.outdir} \
        --enable_phasing \
        --longphase_for_phasing
        """
