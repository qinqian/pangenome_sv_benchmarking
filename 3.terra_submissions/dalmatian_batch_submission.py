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


# Sniffles2 single sample columns
# for GRCh37 and Chm13v2 of all samples
sniffles2_single_mode_t2t_37 = {("minimap2_grch37_noalt_cram", "minimap2_grch37_noalt_crai"): ("sniffles_grch37_vcf", "sniffles_grch37_tbi"),
                                ("minimap2_chm13v2_cram", "minimap2_chm13v2_crai"): ("sniffles_chm13v2_vcf", "sniffles_chm13v2_tbi")}

# NOTE: GRCh38 distributed in different columns
# sniffles2_single_mode_hg38 = {("cram", "crai"): "sample_set"}


def submit_minigraph_gaf_jobs(sample_set_id, minigraph_wdl="minigraph_cram", use_callcache=False):
    submission_id = wm.create_submission(minigraph_wdl, sample_set_id, 'sample_set', expression='this.samples', use_callcache=use_callcache)

def submit_minimap2_cram_jobs(sample_set_id, minimap2_wdl="minimap2_cram_to_cram", use_callcache=False):
    submission_id = wm.create_submission(minimap2_wdl, sample_set_id, 'sample_set', expression='this.samples', use_callcache=use_callcache)


def submit_sniffles2_cram_jobs(sample_set_id, minimap2_wdl="sniffles_workflow", use_callcache=True):
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

    for ((cram, crai), (vcf, tbi)) in sniffles2_single_mode_t2t_37.items():
        print(cram, crai)
        old_config = wm.get_config("sniffles_workflow")
        old_config["inputs"]["SNFWorkflow.boot_disk_size"] = '100'
        old_config["inputs"]["SNFWorkflow.disk_space"] = '200'
        old_config["inputs"]["SNFWorkflow.mosaic"] = ''
        print('-------')
        print(old_config)
        old_config["inputs"]['SNFWorkflow.cram'] = f'this.{cram}'
        old_config["inputs"]['SNFWorkflow.crai'] = f'this.{crai}'
        old_config['outputs']['SNFWorkflow.vcf'] = f'this.{vcf}'
        old_config['outputs']['SNFWorkflow.tbi'] = f'this.{tbi}'
        print(old_config)
        new_config = wm.update_config(old_config)
        submit_sniffles2_cram_jobs("all_samples")

    status = wm.get_submission_status(filter_active=False)
    status.loc[:, 'clean_up'] = False
    # remove all jobs that has not supported ds:Z tag correctly
    status.iloc[-31:, -1] = True
    status.to_csv("job_status.csv")

    #clean_up(status.submission_id[status.clean_up])
    #clean_up(status.submission_id[status.clean_up], dry_run=False)


if __name__ == "__main__":
    main()
