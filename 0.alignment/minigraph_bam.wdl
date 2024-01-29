version 1.0


workflow PangenomeAlignment {
    input {
        File bam
        String sample_id
        File assembly
        String docker_image="qianqin/minigraph:v0.20"
        Int preemptible=3
        Int boot_disk_size=20
        Int disk_space=300
        Int cpu = 24
        Int mem = 64
    }

    call minigraphTask {
        input:
            input_bam=bam,
            sample_id=sample_id,
            docker_image=docker_image,
            assembly=assembly,
            preemptible=preemptible,
            boot_disk_size=boot_disk_size,
            disk_space=disk_space,
            cpu=cpu,
            mem=mem,
    }

    output {
        File gaf = minigraphTask.gaf
    }
}


task minigraphTask {
    input {
        File input_bam
        String sample_id
        File assembly
        String docker_image
        Int preemptible=3
        Int boot_disk_size=10
        Int disk_space=20
        Int cpu = 10
        Int mem = 32
    }

    command <<<
        samtools fasta -@ ~{cpu} ~{input_bam} | minigraph -c -t ~{cpu} ~{assembly} - > ~{sample_id}.gaf
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
        File gaf = "~{sample_id}.gaf"
    }
}
