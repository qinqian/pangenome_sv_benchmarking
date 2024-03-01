version 1.0


workflow ClairWorkflow {
    input {
        File bam
        File bai
        String sample_id
        File assembly
        File fai
        String docker_image="hkubal/clair3:v1.0.5"
        Int preemptible=3
        Int boot_disk_size=20
        Int disk_space=300
        Int cpu = 10
        Int mem = 40
        String model_name = "hifi_revio"
        String platform = "hifi"
    }

    call clairTask {
        input:
            input_bam=bam,
            input_bai=bai,
            sample_id=sample_id,
            docker_image=docker_image,
            assembly=assembly,
            fai=fai,
            preemptible=preemptible,
            boot_disk_size=boot_disk_size,
            disk_space=disk_space,
            cpu=cpu,
            mem=mem,
            model_name=model_name,
            platform=platform
    }

    output {
        File phased_vcf = clairTask.output_vcf_gz
        File phased_vcf_tbi = clairTask.output_vcf_gz_tbi
    }
}


task clairTask {
    input {
        File input_bam
        File input_bai
        String sample_id
        File assembly
        File fai
        String docker_image
        Int preemptible=3
        Int boot_disk_size=10
        Int disk_space=20
        Int cpu = 10
        Int mem = 64
        String model_name = "hifi_revio"
        String platform = "hifi"
    }

    command <<<
	MODEL_NAME="~{model_name}"
        platform="~{platform}"

        if [[ $MODEL_NAME == "r1041_e82_400bps_sup_v430" ]]; then
            git clone https://github.com/nanoporetech/rerio
            python3 rerio/download_model.py --clair3
            #NOTE: https://github.com/nanoporetech/rerio
            #Our nanopore flowcell is R10.4.1 E8.2 (5kHz)
            #Dorado caller might be a different version
            #we used latest model provided by nanopore developer
            model_path=$(ls -d `pwd`/rerio/clair3_models/r1041_e82_400bps_sup_v430/)
        else
            model_path="/opt/models/${MODEL_NAME}"
        fi

        echo $model_path $platform $MODEL_NAME
        mkdir -p ~{sample_id}

        /opt/bin/run_clair3.sh \
        --bam_fn=~{input_bam} \
        --ref_fn=~{assembly} \
        --threads=~{cpu} \
        --platform=${platform} \
        --model_path=${model_path} \
        --output=~{sample_id} \
        --enable_phasing \
        --longphase_for_phasing
    >>>

    runtime {
        disks: "local-disk ~{disk_space} HDD"
        memory: "~{mem} GB"
        cpu: cpu
        preemptible: preemptible
        bootDiskSizeGb: boot_disk_size
        docker: docker_image
    }

    output {
        File output_vcf_gz = "~{sample_id}/phased_merge_output.vcf.gz"
        File output_vcf_gz_tbi = "~{sample_id}/phased_merge_output.vcf.gz.tbi"
    }
}
