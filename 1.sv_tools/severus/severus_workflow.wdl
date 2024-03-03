version 1.0


workflow SeverusWorkflow {
    input {
        File tumor_bam_or_cram
        File tumor_bam_or_cram_index
        String tumor_sample_id

        File? normal_bam_or_cram
        File? normal_bam_or_cram_index
        String? normal_sample_id

        # single sample mode
        File phased_tumor_vcfgz
        File phased_tumor_vcf_index

        # tumor-normal pair mode
        File? phased_normal_vcfgz
        File? phased_normal_vcf_index

        File vntr

        File assembly
        File assembly_index
        String docker_image="qianqin/severus"

        Int preemptible=3
        Int boot_disk_size=20
        Int disk_space=300
        Int cpu = 24
        Int mem = 96
        Boolen use_ssd = True
    }

    Array[String] sample_id = select_all([tumor_sample_id, normal_sample_id])
    Array[File] crams = select_all([tumor_bam_or_cram, normal_bam_or_cram])
    Array[File] cram_indexes = select_all([tumor_bam_or_cram_index, normal_bam_or_cram_index])
    File phased_vcfgz = select_first([phased_normal_vcfgz, phased_tumor_vcfgz])
    File phased_vcfgz_index = select_first([phased_normal_vcf_index, phased_tumor_vcf_index])

    scatter (idx in range(length(crams))) {
      call whatshapTask {
        input:
          cram = crams[idx],
          crai = cram_indexes[idx], 
          phased_vcfgz = phased_vcfgz,
          phased_vcfgz_index = phased_vcfgz_index,
          sample_id = sample_id[idx],
          assembly = assembly,
          assembly_index = assembly_index,
          cpu = cpu,
          mem = mem,
          use_ssd = use_ssd,
          preemptible=preemptible,
          boot_disk_size=boot_disk_size,
          disk_space=disk_space
      }
    }

    Array[File] haplocrams = select_all(whatshapTask.haplocram)
    Array[File] haplocrais = select_all(whatshapTask.haplocrai)

    if (length(haplocrams) > 1) {
        call severusTask as severusTumorNormal {
            input:
                tumor_cram=haplocrams[0],
                tumor_crai=haplocrais[0],
                normal_cram=haplocrams[1],
                normal_crai=haplocrais[1],
                phased_vcf=phased_vcfgz,
                phased_vcf_index=phased_vcfgz_index,
                vntr=vntr,

                sample_id=sample_id[0],
                docker_image=docker_image,
                preemptible=preemptible,
                boot_disk_size=boot_disk_size,
                disk_space=disk_space,
                cpu=cpu,
                mem=mem,
        }
    } 

    if (length(haplocrams) == 1) {
        call severusTask as severusTumor {
            input:
                tumor_cram=haplocrams[0],
                tumor_crai=haplocrais[0],
                phased_vcf=phased_vcfgz,
                phased_vcf_index=phased_vcfgz_index,
                vntr=vntr,

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
        #File vcf_all = select_first([severusTumor.vcf_all, severusTumorNormal.vcf_all])
        #File vcf_som = select_first([severusTumor.vcf_som, severusTumorNormal.vcf_som])
        Array[File] vcf_all = select_first([severusTumor.vcf, severusTumorNormal.vcf])
    }
}


task whatshapTask {
    input {
        File cram
        File crai
        File phased_vcfgz
        File phased_vcfgz_index

        String sample_id
        File assembly
        File assembly_index
        String docker_image="qianqin/severus"
        Int preemptible=3
        Int boot_disk_size=10
        Int disk_space=20
        Int cpu = 10
        Int mem = 64
        Boolen use_ssd = False
    }

    command <<<
	whatshap haplotag --reference ~{assembly} ~{phased_vcfgz} ~{cram} -o ~{sample_id}.haplotagged.cram --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads=~{cpu}
	samtools index ~{sample_id}.haplotagged.cram
    >>>

    runtime {
        disks: "local-disk " + ceil(size(cram, "GB)+size(crai, "GB)+size(phased_vcfgz, "GB")+size(phased_vcfgz_index, "GB")+size(assembly, "GB")+disk_space) + " " + (if use_ssd then "SSD" else "HDD")
        memory: "~{mem} GB"
        cpu: cpu
        preemptible: preemptible
        bootDiskSizeGb: boot_disk_size
        docker: docker_image
        maxRetries: 3
    }

    output {
        File haplocram = "~{sample_id}.haplotagged.cram"
        File haplocrai = "~{sample_id}.haplotagged.cram.crai"
    }
}


task severusTask {
    input {
        File tumor_cram
        File tumor_crai

        File? normal_cram
        File? normal_crai

        File phased_vcf
        File phased_vcf_index

        String sample_id
        File vntr

        String docker_image
        Int preemptible=3
        Int boot_disk_size=10
        Int disk_space=20
        Int cpu = 10
        Int mem = 32
    }

    command <<<
        severus --target-bam ~{tumor_cram} ~{"--control-bam " + normal_cram}  --out-dir ~{sample_id}_severus_out -t ~{cpu} --phasing-vcf ~{phased_vcf} \
                ~{"--vntr-bed " + vntr}
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
        #File vcf_all = glob("~{sample_id}_severus_out/all_SVs/*vcf")[0]
        #File vcf_som = glob("~{sample_id}_severus_out/somatic_SVs/*vcf")[0]
        Array[File] vcf = flatten([glob("~{sample_id}_severus_out/all_SVs/*vcf"), glob("~{sample_id}_severus_out/somatic_SVs/*vcf")])
    }
}
