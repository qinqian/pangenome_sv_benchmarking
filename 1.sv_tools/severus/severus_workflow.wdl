version 1.0


workflow SeverusWorkflow {
    input {
        File tumor_bam_or_cram
        File tumor_bam_or_cram_index
        File? normal_bam_or_cram
        File? normal_bam_or_cram_index
        String sample_id

        File assembly
        String docker_image="qianqin/severus"

        Int preemptible=3
        Int boot_disk_size=20
        Int disk_space=300
        Int cpu = 24
        Int mem = 96
    }

    if (defined(normal_bam_or_cram) && defined(normal_bam_or_cram_index)) {
       Array[File] crams = [tumor_bam_or_cram, normal_bam_or_cram]
       Array[File] cram_indexes = [tumor_bam_or_cram_index, normal_bam_or_cram_index]
    } else {
       Array[File] crams = [tumor_bam_or_cram]
       Array[File] cram_indexes = [tumor_bam_or_cram_index]
    }

    scatter (sample in crams) {
      call whatshapTask {
        input:
          cram = cram
      }
    }

    call severusTask {
        input:
            input_bam=whatshapTask.bam,
            input_bai=whatshapTask.bai,
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


task whatshapTask {
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
        File vntr
        File phased_vcf
    }

    command <<<
	whatshap haplotag --reference ~{assembly} ~{phased_vcf} ~{input_bam} -o ~{sample_id}.haplotagged.bam --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads=~{cpu}
	samtools index ~{sample_id}.haplotagged.bam
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
        File haplobam = "~{sample_id}.haplotagged.bam"
        File haplobam = "~{sample_id}.haplotagged.bam.bai"
    }
}


task severusTask {
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
        File vntr
        File phased_vcf
    }

    command <<<
	whatshap haplotag --reference ~{assembly} ~{phased_vcf} ~{input_bam} -o ~{sample_id}.haplotagged.bam --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads=~{cpu}

	whatshap haplotag --reference ~{assembly} ~{phased_vcf} ~{input_bam} -o ~{sample_id}.haplotagged.bam --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads=~{cpu}

	samtools index ~{sample_id}.haplotagged.bam

        severus --target-bam ~{sample_id}.haplotagged.bam --out-dir ~{sample_id}_severus_out -t ~{cpu} --phasing-vcf ~{phased_vcf} \
        	    --vntr-bed ~{vntr}
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
        File haplobam = "~{sample_id}.haplotagged.bam"
    }
}
