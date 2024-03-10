version 1.0


workflow LineargenomeAlignment {
    input {
        File cram
        File cram_reference
        String sample_id
        File assembly
        String docker_image="docker.io/qianqin/minimap2_ds"
        Int preemptible=3
        Int boot_disk_size=20
        Int disk_space=300
        Int cpu = 24
        Int mem = 64
        String mode = "map-hifi"
        Boolean use_ssd = true
    }

    call minimapTask {
        input:
            input_cram=cram,
            input_cram_reference=cram_reference,
            sample_id=sample_id,
            docker_image=docker_image,
            assembly=assembly,
            preemptible=preemptible,
            boot_disk_size=boot_disk_size,
            disk_space=disk_space,
            cpu=cpu,
            mem=mem,
            mode=mode,
					  use_ssd=use_ssd
    }

    output {
        File paf = minimapTask.paf
    }
}


task minimapTask {
    input {
        File input_cram
        File input_cram_reference
        String mode
        String sample_id
        File assembly
        String docker_image
        Int preemptible=3
        Int boot_disk_size=10
        Int disk_space=20
        Int cpu = 10
        Int mem = 32
        Boolean use_ssd
    }

    command <<<
        samtools collate -@ ~{cpu} -Oun128 --reference ~{input_cram_reference} ~{input_cram} | samtools fastq -@ ~{cpu} - | minimap2 --ds -c -x ~{mode} -t ~{cpu} ~{assembly} - > ~{sample_id}.paf 
    >>>

    runtime {
        disks: "local-disk ~{disk_space}" + " " + (if use_ssd then "SSD" else "HDD")
        memory: "~{mem} GB"
        cpu: cpu
        preemptible: preemptible
        bootDiskSizeGb: boot_disk_size
        docker: docker_image
    }

    output {
        File paf = "~{sample_id}.paf"
    }
}
