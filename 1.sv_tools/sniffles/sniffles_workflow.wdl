version 1.0


workflow SNFWorkflow {
    input {
        File cram
        File crai

        String sample_id
        String docker_image="qianqin/snf2.2"
        String? mosaic

        Int preemptible=3
        Int boot_disk_size=20
        Int disk_space=300
        Int cpu = 16
        Int mem = 96
    }

    call SNFTask {
        input:
            cram=cram,
            crai=crai,
            mosaic=mosaic,
            sample_id=sample_id,
            docker_image=docker_image,
            preemptible=preemptible,
            boot_disk_size=boot_disk_size,
            disk_space=disk_space,
            cpu=cpu,
            mem=mem
    }

    output {
        File vcf = SNFTask.vcf
    }
}


task SNFTask {
    input {
        File cram
        File crai
        String sample_id
        String? mosaic
        String docker_image
        Int preemptible=3
        Int boot_disk_size=10
        Int disk_space=20
        Int cpu = 10
        Int mem = 32
    }

    command <<<
        sniffles --threads ~{cpu} -i ~{cram} -v ~{sample_id}.vcf.gz --output-rnames --sample-id ~{sample_id} ~{" " + mosaic}
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
        File vcf = "~{sample_id}.vcf.gz"
        File tbi = "~{sample_id}.vcf.gz.tbi"
    }
}
