from mgenepy import terra
import dalmatian as dm
from depmap_omics_upload.mgenepy import terra as terra_cleanup
from mgenepy.utils import helper as h 
import subprocess

pangenome_sv_workspace = "broad-firecloud-ccle/pangenome"
wm = dm.WorkspaceManager(pangenome_sv_workspace)

# Workspace column name to corresponding assembly
# This minigraph_cram will output pafs
minigraph_cram_output_to_assemblies = {"minigraph_chm13v2_gaf": "\"gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/chm13v2.0.fa\"", 
                                       "minigraph_chm13graph_gaf": "\"gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/chm13-90c.r518.gfa\"",
                                       "minigraph_grch38graph_gaf": "\"gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/GRCh38-90c.r518.gfa\"",
                                       "minigraph_grch38_noalt_gaf": "\"gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna\"", 
                                       "minigraph_grch37graph_gaf": "\"gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/GRCh37-91c.r559.gfa.gz\"",
                                       "minigraph_grch37linear_gaf": "\"gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/hs37d5.fa\""}
#TODO: automatically generate this in the future
remaining_task = ["minigraph_chm13graph_gaf", "minigraph_grch38graph_gaf", "minigraph_grch37linear_gaf"]


# Workspace column name to corresponding assembly
# This minigraph_cram_to_cram will output pafs
minimap2_cram_output_to_assemblies = {("minimap2_chm13v2_cram", "minimap2_chm13v2_crai"): "\"gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/chm13v2.0.fa\"", 
                                      ("minimap2_grch37_noalt_cram", "minimap2_grch37_noalt_crai"): "\"gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/hs37d5.fa\""}


# Clair3 for severus
clair3_phasing = {("minimap2_chm13v2_cram", "minimap2_chm13v2_crai", "\"gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/chm13v2.0.fa\"", "\"gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/chm13v2.0.fa.fai\""): ("clair_chm13v2_phased_vcf", "clair_chm13v2_phased_vcf_tbi"),
                  ("minimap2_grch37_noalt_cram", "minimap2_grch37_noalt_crai", "\"gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/hs37d5.fa\"", "\"gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/hs37d5.fa.fai\""): ("clair_grch37_phased_vcf", "clair_grch37_phased_vcf_tbi")}

# NOTE: GRCh38 distributed in different columns
clair3_phasing_grch38 = {("cram", "crai", "nanopore_merged_techreps", "\"gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna\"", "\"gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai\""): ("clair_phased_vcf", "clair_phased_vcf_tbi", '\"ont\"', '\"r1041_e82_400bps_sup_v430\"'),
                         ("bam", "bai", "all_pacbio_hg002_4cancerpairedcells", "\"gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna\"", "\"gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai\""): ("clair_phased_vcf", "clair_phased_vcf_tbi", '\"hifi\"', '\"hifi_revio\"')}


# Sniffles2 single sample columns
# for GRCh37 and Chm13v2 of all samples
sniffles2_single_mode_t2t_37 = {("minimap2_grch37_noalt_cram", "minimap2_grch37_noalt_crai"): ("sniffles_grch37_vcf", "sniffles_grch37_tbi"),
                                ("minimap2_chm13v2_cram", "minimap2_chm13v2_crai"): ("sniffles_chm13v2_vcf", "sniffles_chm13v2_tbi")}

# NOTE: GRCh38 distributed in different columns
sniffles2_single_mode_hg38 = {("cram", "crai", "nanopore_merged_techreps"): ("sniffles_grch38_noalt_vcf", "sniffles_grch38_noalt_tbi"),
                              ("bam", "bai", "all_pacbio_hg002_4cancerpairedcells"): ("sniffles_grch38_noalt_vcf", "sniffles_grch38_noalt_tbi")}

# severus for GRCh37 and Chm13v2
# for both tumor-only and tumor-normal pair mode
severus_tumor_only_or_pair_grch37 = {("minimap2_grch37_noalt_cram", "minimap2_grch37_noalt_crai", "clair_grch37_phased_vcf", "clair_grch37_phased_vcf_tbi", "\"gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/human_hs37d5.trf.bed\"", "\"gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/hs37d5.fa\"", "\"gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/hs37d5.fa.fai\""): ("severus_grch37_singlelsample_vcf_all", "all_samples"),  # tumor only use phased tumor vcf
                                     ("minimap2_grch37_noalt_cram", "minimap2_grch37_noalt_crai", "paired_normal_grch37_cram", "paired_normal_grch37_crai", "clair_phased_grch37_normal_vcf", "clair_phased_grch37_normal_vcf_index", "\"gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/human_hs37d5.trf.bed\"", "\"gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/hs37d5.fa\"", "\"gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/hs37d5.fa.fai\""): ("severus_grch37_tumor_normal_pair_vcf_all", "all_tumor_normal_pairs")}  # tumor-normal pair only use phased normal vcf
                            # ("minimap2_chm13v2_cram", "minimap2_chm13v2_crai", "clair_chm13v2_phased_vcf", "clair_chm13v2_phased_vcf_tbi", "\"gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/chm13v2.0.fa\"", "\"gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/chm13v2.0.fa.fai\""): "severus_chm13v2_singlelsample_vcf_all"}


def submit_minigraph_gaf_jobs(sample_set_id, minigraph_wdl="minigraph_cram", use_callcache=False):
    submission_id = wm.create_submission(minigraph_wdl, sample_set_id, 'sample_set', expression='this.samples', use_callcache=use_callcache)

def submit_minimap2_cram_jobs(sample_set_id, minimap2_wdl="minimap2_cram_to_cram", use_callcache=False):
    submission_id = wm.create_submission(minimap2_wdl, sample_set_id, 'sample_set', expression='this.samples', use_callcache=use_callcache)

def submit_sniffles2_cram_jobs(sample_set_id, minimap2_wdl="sniffles_workflow", use_callcache=True):
    submission_id = wm.create_submission(minimap2_wdl, sample_set_id, 'sample_set', expression='this.samples', use_callcache=use_callcache)

def submit_clair3_cram_jobs(sample_set_id, minimap2_wdl="clair3_bam", use_callcache=False):
    submission_id = wm.create_submission(minimap2_wdl, sample_set_id, 'sample_set', expression='this.samples', use_callcache=use_callcache)

def submit_severus_cram_jobs(sample_set_id, minimap2_wdl="severus_workflow", use_callcache=True):
    submission_id = wm.create_submission(minimap2_wdl, sample_set_id, 'sample_set', expression='this.samples', use_callcache=use_callcache)

def clean_up(to_be_clean_submission_ids, dry_run=True):
    print("cleaning workspaces")
    subprocess.call("rm -f clean.logs", shell=True)
    for submission_id in to_be_clean_submission_ids:
        if dry_run:
            subprocess.call(f"gsutil du -sh gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/{submission_id} >> clean.logs", shell=True)
        else:
            subprocess.call(f"gsutil rm -r gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/{submission_id}", shell=True)

def main():
    samples_df = wm.get_samples()
    #wm.update_sample_set('all_samples', samples_df.index)
    #wm.update_sample_set('all_tumor_normal_pairs', ["HCC1395", "COLO829", "COLO829_ONT"])
    print(wm.get_sample_sets())

    #for output_col in remaining_task:
    #    assembly = minigraph_cram_output_to_assemblies[output_col]
    #    print(output_col, assembly)
    #    old_config = wm.get_config("minigraph_cram")
    #    print(old_config)
    #    old_config["inputs"]['PangenomeAlignment.assembly'] = assembly
    #    print(old_config)
    #    old_config['outputs']['PangenomeAlignment.gaf'] = f'this.{output_col}'
    #    new_config = wm.update_config(old_config)
    #    submit_minigraph_gaf_jobs("nanopore_merged_techreps")
    #    #break

    #for output_col, assembly in minimap2_cram_output_to_assemblies.items():
    #    old_config = wm.get_config("minimap2_cram_to_cram")
    #    print('-------')
    #    print(old_config)
    #    old_config["inputs"]['LineargenomeAlignment.assembly'] = assembly
    #    old_config['outputs']['LineargenomeAlignment.cram'] = f'this.{output_col[0]}'
    #    old_config['outputs']['LineargenomeAlignment.crai'] = f'this.{output_col[1]}'
    #    print(old_config)
    #    new_config = wm.update_config(old_config)
    #    submit_minimap2_cram_jobs("nanopore_merged_techreps")

    #for ((cram, crai), (vcf, tbi)) in sniffles2_single_mode_t2t_37.items():
    #    print(cram, crai)
    #    old_config = wm.get_config("sniffles_workflow")
    #    old_config["inputs"]["SNFWorkflow.boot_disk_size"] = '100'
    #    old_config["inputs"]["SNFWorkflow.disk_space"] = '200'
    #    old_config["inputs"]["SNFWorkflow.mosaic"] = ''
    #    print('-------')
    #    print(old_config)
    #    old_config["inputs"]['SNFWorkflow.cram'] = f'this.{cram}'
    #    old_config["inputs"]['SNFWorkflow.crai'] = f'this.{crai}'
    #    old_config['outputs']['SNFWorkflow.vcf'] = f'this.{vcf}'
    #    old_config['outputs']['SNFWorkflow.tbi'] = f'this.{tbi}'
    #    print(old_config)
    #    new_config = wm.update_config(old_config)
    #    submit_sniffles2_cram_jobs("all_samples")

    #for ((cram, crai, sample_set), (vcf, tbi)) in sniffles2_single_mode_hg38.items():
    #    print(cram, crai)
    #    old_config = wm.get_config("sniffles_workflow")
    #    old_config["inputs"]["SNFWorkflow.boot_disk_size"] = '100'
    #    old_config["inputs"]["SNFWorkflow.disk_space"] = '200'
    #    old_config["inputs"]["SNFWorkflow.mosaic"] = ''
    #    print('-------')
    #    print(old_config)
    #    old_config["inputs"]['SNFWorkflow.cram'] = f'this.{cram}'
    #    old_config["inputs"]['SNFWorkflow.crai'] = f'this.{crai}'
    #    old_config['outputs']['SNFWorkflow.vcf'] = f'this.{vcf}'
    #    old_config['outputs']['SNFWorkflow.tbi'] = f'this.{tbi}'
    #    print(old_config)
    #    new_config = wm.update_config(old_config)
    #    submit_sniffles2_cram_jobs(sample_set)

    #for ((cram, crai, fa, fai), (vcf, tbi)) in clair3_phasing.items():
    #    print(cram, crai)
    #    old_config = wm.get_config("clair3_bam")
    #    print(old_config)
    #    old_config["inputs"]["ClairWorkflow.assembly"] = fa
    #    old_config["inputs"]["ClairWorkflow.fai"] = fai
    #    old_config["inputs"]['ClairWorkflow.bam'] = f'this.{cram}'
    #    old_config["inputs"]['ClairWorkflow.bai'] = f'this.{crai}'

    #    old_config["inputs"]['ClairWorkflow.model_name'] = 'hifi_revio'
    #    old_config["inputs"]['ClairWorkflow.platform'] = 'hifi'

    #    old_config['outputs']['ClairWorkflow.phased_vcf'] = f'this.{vcf}'
    #    old_config['outputs']['ClairWorkflow.phased_vcf_tbi'] = f'this.{tbi}'
    #    print(old_config)
    #    new_config = wm.update_config(old_config)
    #    submit_clair3_cram_jobs("all_pacbio_hg002_4cancerpairedcells")

    for ((cram, crai, sample_set, fa, fai), (vcf, tbi, platform, model_name)) in clair3_phasing_grch38.items():
        print(cram, crai)
        old_config = wm.get_config("clair3_bam")
        print(old_config)
        old_config["inputs"]["ClairWorkflow.assembly"] = fa
        old_config["inputs"]["ClairWorkflow.fai"] = fai
        old_config["inputs"]['ClairWorkflow.bam'] = f'this.{cram}'
        old_config["inputs"]['ClairWorkflow.bai'] = f'this.{crai}'
        old_config['outputs']['ClairWorkflow.phased_vcf'] = f'this.{vcf}'
        old_config['outputs']['ClairWorkflow.phased_vcf_tbi'] = f'this.{tbi}'
        if "hifi" in platform : # skip pacbio which we have right model
            continue
        old_config["inputs"]['ClairWorkflow.model_name'] = model_name
        old_config["inputs"]['ClairWorkflow.platform'] = platform
        print(old_config)
        new_config = wm.update_config(old_config)
        submit_clair3_cram_jobs(sample_set)

    for ((cram, crai, fa, fai), (vcf, tbi)) in clair3_phasing.items():
        print(cram, crai)
        old_config = wm.get_config("clair3_bam")
        print(old_config)
        old_config["inputs"]["ClairWorkflow.assembly"] = fa
        old_config["inputs"]["ClairWorkflow.fai"] = fai
        old_config["inputs"]['ClairWorkflow.bam'] = f'this.{cram}'
        old_config["inputs"]['ClairWorkflow.bai'] = f'this.{crai}'

        old_config["inputs"]['ClairWorkflow.model_name'] = '\"r1041_e82_400bps_sup_v430\"'
        old_config["inputs"]['ClairWorkflow.platform'] = '\"ont\"'

        old_config['outputs']['ClairWorkflow.phased_vcf'] = f'this.{vcf}'
        old_config['outputs']['ClairWorkflow.phased_vcf_tbi'] = f'this.{tbi}'
        print(old_config)
        new_config = wm.update_config(old_config)
        submit_clair3_cram_jobs("nanopore_merged_techreps")

    #for ((cram, crai, sample_set, fa, fai), (vcf, tbi)) in clair3_phasing_grch38.items():
    #    print(cram, crai)
    #    old_config = wm.get_config("clair3_bam")
    #    print(old_config)
    #    old_config["inputs"]["ClairWorkflow.assembly"] = fa
    #    old_config["inputs"]["ClairWorkflow.fai"] = fai
    #    old_config["inputs"]['ClairWorkflow.bam'] = f'this.{cram}'
    #    old_config["inputs"]['ClairWorkflow.bai'] = f'this.{crai}'
    #    old_config['outputs']['ClairWorkflow.phased_vcf'] = f'this.{vcf}'
    #    old_config['outputs']['ClairWorkflow.phased_vcf_tbi'] = f'this.{tbi}'
    #    print(old_config)
    #    new_config = wm.update_config(old_config)
    #    submit_clair3_cram_jobs(sample_set)

    #for (inputs_severus, (vcf, sample_set)) in severus_tumor_only_or_pair_grch37.items():
    #    old_config = wm.get_config("severus_workflow")
    #    print(inputs_severus)
    #    if len(inputs_severus) == 7: # tumor-only mode
    #        old_config["inputs"]["SeverusWorkflow.tumor_bam_or_cram"] = f'this.{inputs_severus[0]}'
    #        old_config["inputs"]["SeverusWorkflow.tumor_bam_or_cram_index"] = f'this.{inputs_severus[1]}'

    #        old_config["inputs"]["SeverusWorkflow.phased_tumor_vcfgz"] = f'this.{inputs_severus[2]}'
    #        old_config["inputs"]["SeverusWorkflow.phased_tumor_vcf_index"] = f'this.{inputs_severus[3]}'
    #        old_config["inputs"]["SeverusWorkflow.phased_normal_vcfgz"] = f''
    #        old_config["inputs"]["SeverusWorkflow.phased_normal_vcf_index"] = f''

    #        old_config["inputs"]["SeverusWorkflow.vntr"] = inputs_severus[4]
    #        old_config["inputs"]["SeverusWorkflow.assembly"] = inputs_severus[5]
    #        old_config["inputs"]["SeverusWorkflow.assembly_index"] = inputs_severus[6]
    #        old_config["inputs"]["SeverusWorkflow.boot_disk_size"] = '350'
    #        old_config["inputs"]["SeverusWorkflow.disk_space"] = '500'
    #        old_config['outputs']['SeverusWorkflow.vcf_all'] = f'this.{vcf}'
    #        old_config["inputs"]["SeverusWorkflow.normal_bam_or_cram"] = ''
    #        old_config["inputs"]["SeverusWorkflow.normal_bam_or_cram_index"] = ''
    #        old_config["inputs"]["SeverusWorkflow.normal_sample_id"] = ''
    #        print(old_config)
    #    else: # tumor-normal pair
    #        continue
    #        #old_config["inputs"]["SeverusWorkflow.tumor_bam_or_cram"] = f'this.{inputs_severus[0]}'
    #        #old_config["inputs"]["SeverusWorkflow.tumor_bam_or_cram_index"] = f'this.{inputs_severus[1]}'
    #        #old_config["inputs"]["SeverusWorkflow.normal_bam_or_cram"] = f'this.{inputs_severus[2]}'
    #        #old_config["inputs"]["SeverusWorkflow.normal_bam_or_cram_index"] = f'this.{inputs_severus[3]}'
    #        #old_config["inputs"]["SeverusWorkflow.normal_sample_id"] = f'this.normal_sample_id'
    #        #old_config["inputs"]["SeverusWorkflow.boot_disk_size"] = '200'
    #        #old_config["inputs"]["SeverusWorkflow.disk_space"] = '200'

    #        #old_config["inputs"]["SeverusWorkflow.phased_normal_vcfgz"] = f'this.{inputs_severus[4]}'
    #        #old_config["inputs"]["SeverusWorkflow.phased_normal_vcf_index"] = f'this.{inputs_severus[5]}'

    #        #old_config["inputs"]["SeverusWorkflow.vntr"] = inputs_severus[6]
    #        #old_config["inputs"]["SeverusWorkflow.assembly"] = inputs_severus[7]
    #        #old_config["inputs"]["SeverusWorkflow.assembly_index"] = inputs_severus[8]
    #        #old_config['outputs']['SeverusWorkflow.vcf_all'] = f'this.{vcf}'
    #    new_config = wm.update_config(old_config)
    #    submit_severus_cram_jobs(sample_set)

    status = wm.get_submission_status(filter_active=False)
    status.loc[:, 'clean_up'] = False
    # remove all jobs that has not supported ds:Z tag correctly
    status.iloc[-31:, -1] = True
    status.to_csv("job_status.csv")

    #clean_up(status.submission_id[status.clean_up])
    #clean_up(status.submission_id[status.clean_up], dry_run=False)


if __name__ == "__main__":
    main()
