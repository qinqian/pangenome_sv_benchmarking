rule clair3:
    threads: 24
    input:
        cram = "output/align/{cell_line}/{assembly}.cram"
    output:
        outdir = directory("output/clair3/{cell_line}/{assembly}")
    conda: "clair3"
    shell:
        """
        model_path=$(ls -d `pwd`/models/hifi_revio/)
        platform="hifi"
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
