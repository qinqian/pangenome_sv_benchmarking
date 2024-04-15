version 1.0


workflow NanomonsvWorkflow {
    input {
        File tumor_bam_or_cram
        File tumor_bam_or_cram_index
        File fasta

        String tumor_sample_id

        File normal_bam_or_cram
        File normal_bam_or_cram_index
        String? normal_sample_id

        String docker_image="qianqin/nanomonsv"
        Int preemptible=3
        Int boot_disk_size=20
        Int disk_space=300
        Int cpu = 1
        Int mem = 96
        Boolean use_ssd = true
    }

    call nanomonsvTask as nanomonsvTumorNormal {
        input:
            tumor_cram=tumor_bam_or_cram,
            tumor_crai=tumor_bam_or_cram_index,
            normal_cram=normal_bam_or_cram,
            normal_crai=normal_bam_or_cram_index,
            fasta=fasta,

            sample_id=tumor_sample_id,
            docker_image=docker_image,
            preemptible=preemptible,
            boot_disk_size=boot_disk_size,
            disk_space=disk_space,
            cpu=cpu,
            mem=mem,
    }

    output {
        File vcf = nanomonsvTumorNormal.vcf
    }
}


task nanomonsvTask {
    input {
        File tumor_cram
        File tumor_crai
        File normal_cram
        File normal_crai
        File fasta

        String sample_id

        String docker_image
        Int preemptible=3
        Int boot_disk_size=10
        Int disk_space=20
        Int cpu = 10
        Int mem = 32
    }

    command <<<
        nanomonsv parse ~{tumor_cram} ~{sample_id}.tumor
        nanomonsv parse ~{normal_cram} ~{sample_id}.normal
        nanomonsv get ~{sample_id}.tumor ~{tumor_cram} ~{fasta} --control_prefix ~{sample_id}.normal --control_bam ~{normal_cram}
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
        File vcf = "~{sample_id}.tumor.nanomonsv.result.txt"
    }
}
