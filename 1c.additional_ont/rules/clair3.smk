rule clair3:
    threads: 24
    resources:
        mem_mb=96000
    input:
        cram = "{cell_line}_{assembly}.cram",
        crai = "{cell_line}_{assembly}.cram.crai"
    output:
        outdir = directory("output/clair3/{cell_line}_{assembly}")
    conda: "clair3"
    shell:
        """
        if [[ {input.cram} =~ "hifi" ]]; then
            model_path=$(ls -d ../1a.alignment_sv_tools/models/hifi_revio/)
            platform="hifi"
        else
            model_path=$(ls -d ../1a.alignment_sv_tools/rerio/clair3_models/r1041_e82_400bps_sup_v430/)
            platform="ont"
        fi
        mkdir -p {output.outdir}
        run_clair3.sh \
        --bam_fn={input.cram} \
        --ref_fn=../1a.alignment_sv_tools/{wildcards.assembly}.fa \
        --threads={threads} \
        --platform=${{platform}} \
        --model_path=${{model_path}} \
        --output={output.outdir} \
        --enable_phasing \
        --longphase_for_phasing
        """
