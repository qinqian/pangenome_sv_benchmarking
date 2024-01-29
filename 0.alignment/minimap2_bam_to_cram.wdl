version 1.0


workflow LineargenomeAlignment {
    input {
        File bam
        String sample_id
        File assembly
        String docker_image="trinityctat/pbfusion:v0.3.1"
        Int preemptible=3
        Int boot_disk_size=20
        Int disk_space=300
        Int cpu = 2
        Int mem = 16
        String mode = "map-hifi"
    }

    call minimapTask {
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
            mode=mode
    }

    output {
        File cram = minimapTask.cram
        File crai = minimapTask.crai
    }
}


task minimapTask {
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
    }

    command <<<
        samtools fasta -@ ~{cpu} ~{input_bam} | minimap2 -ax ~{mode} -t ~{cpu} ~{assembly} - | samtools sort -l 1 -@ ~{cpu} --output-fmt CRAM --reference ~{assembly} > ~{sample_id}.cram
        samtools index ~{sample_id}.cram
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
        File cram = "~{sample_id}.cram"
        File crai = "~{sample_id}.cram.crai"
    }
}
