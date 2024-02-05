version 1.1

workflow SNFWorkflow {
    input {
        File tumor_bam_or_cram
        File tumor_bam_or_cram_index
        String tumor_sample_id

        File? normal_bam_or_cram
        File? normal_bam_or_cram_index
        String? normal_sample_id

        String? mosaic

        String docker_image="qianqin/snf2.2"

        Int preemptible=3
        Int boot_disk_size=20
        Int disk_space=300
        Int cpu = 24
        Int mem = 96
    }

    Array[String] sample_id = select_all([tumor_sample_id, normal_sample_id])
    Array[File] crams = select_all([tumor_bam_or_cram, normal_bam_or_cram])
    Array[File] cram_indexes = select_all([tumor_bam_or_cram_index, normal_bam_or_cram_index])

    if (length(crams) > 1) {
        scatter (idx in range(length(crams))) {
          call snfTask as snfTaskProc {
            input:
              cram = crams[idx],
              crai = cram_indexes[idx], 
              sample_id = sample_id[idx],
          }
        }
        call snifflesTask as snfTumorNormal {
            input:
                snfs = snfTaskProc.snf,

                sample_id=sample_id[0],
                docker_image=docker_image,
                preemptible=preemptible,
                boot_disk_size=boot_disk_size,
                disk_space=disk_space,
                cpu=cpu,
                mem=mem,
        }
    }


    if (length(crams) == 1) {
        call SNFTask as snfTumor {
            input:
                cram=crams[0],
                crai=cram_indexes[0],
                mosaic=mosaic,

                sample_id=sample_id[0],
                docker_image=docker_image,
                preemptible=preemptible,
                boot_disk_size=boot_disk_size,
                disk_space=disk_space,
                cpu=cpu,
                mem=mem,
        }
    }

    output {
        File snf_vcf = select_first([snfTumor.vcf, snfTumorNormal.vcf])
        File snf_vcftbi = select_first([snfTumor.tbi, snfTumorNormal.tbi])
    }
}


task snfTask {
    input {
        File cram
        File crai

        String sample_id
        String docker_image="qianqin/snf2.2"
        Int preemptible=3
        Int boot_disk_size=10
        Int disk_space=20
        Int cpu = 10
        Int mem = 32
    }

    command <<<
	sniffles --input ~{cram} --snf ~{sample_id}.snf --threads ~{cpu}
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
        File snf = "~{sample_id}.snf"
    }
}


task snifflesTask {
    input {
        Array[File] snfs

        String sample_id

        String docker_image
        Int preemptible=3
        Int boot_disk_size=10
        Int disk_space=20
        Int cpu = 10
        Int mem = 32
    }

    command <<<
        sniffles --input ~{sep(' ', snfs)} --output-rnames --vcf ~{sample_id}-multisample.vcf.gz
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
        File vcf = "~{sample_id}-multisample.vcf.gz"
        File tbi = "~{sample_id}-multisample.vcf.gz.tbi"
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
