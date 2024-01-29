version 1.0


workflow DeepVariantWorkflow {
    input {
        File bam
        String sample_id
        File assembly
        File par_region_bed
        String docker_image="google/deepvariant:1.6.0"
        Int preemptible=3
        Int boot_disk_size=20
        Int disk_space=300
        Int cpu = 48
        Int mem = 96
        String mode = "PACBIO"
    }

    call deepvariantTask {
        input:
            input_bam=bam,
            sample_id=sample_id,
            docker_image=docker_image,
            assembly=assembly,
            par_region_bed=par_region_bed,
            preemptible=preemptible,
            boot_disk_size=boot_disk_size,
            disk_space=disk_space,
            cpu=cpu,
            mem=mem,
            mode=mode
    }

    output {
        File vcf = deepvariantTask.vcf
    }
}


task deepvariantTask {
    input {
        File input_bam
        String mode
        String sample_id
        File assembly
        String docker_image
        Int preemptible=3
        Int boot_disk_size=10
        Int disk_space=20
        Int cpu = 10
        Int mem = 32
        File par_region_bed
    }

    command <<<
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=~{mode} \
        --ref=~{assembly} \
        --reads=~{input_bam} \
        --output_vcf=~{sample_id}.vcf \
        --num_shards=~{cpu} \
        --logging_dir=logs \
        --haploid_contigs="chrX,chrY" \
        --par_regions_bed=~{par_region_bed} \
        --dry_run=false
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
        File vcf = "~{sample_id}.vcf"
    }
}
