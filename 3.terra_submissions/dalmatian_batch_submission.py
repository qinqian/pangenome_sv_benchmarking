from mgenepy import terra
import dalmatian as dm
from depmap_omics_upload.mgenepy import terra as terra_cleanup
from mgenepy.utils import helper as h 

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


def submit_minigraph_gaf_jobs(sample_set_id, minigraph_wdl="minigraph_cram", use_callcache=False):
    submission_id = wm.create_submission(minigraph_wdl, sample_set_id, 'sample_set', expression='this.samples', use_callcache=use_callcache)
    #await terra.waitForSubmission(pangenome_sv_workspace, submission_id)


#def clean_up():
#    if doCleanup:
#        print("cleaning workspaces")
#        torm = await terra_cleanup.deleteHeavyFiles("broad-firecloud-ccle/DEV_DepMap_WGS_CN")
#        h.parrun(['gsutil rm '+i for i in torm], cores=8)
#        terra_cleanup.removeFromFailedWorkflows("broad-firecloud-ccle/DEV_DepMap_WGS_CN", dryrun=False)


def main():
    print(wm.get_samples())
    print(wm.get_sample_sets())

    #for output_col, assembly in minigraph_cram_output_to_assemblies.items():
    for output_col in remaining_task:
        assembly = minigraph_cram_output_to_assemblies[output_col]
        print(output_col, assembly)
        old_config = wm.get_config("minigraph_cram")
        print(old_config)
        old_config["inputs"]['PangenomeAlignment.assembly'] = assembly
        print(old_config)
        old_config['outputs']['PangenomeAlignment.gaf'] = f'this.{output_col}'
        new_config = wm.update_config(old_config)
        submit_minigraph_gaf_jobs("nanopore_merged_techreps")
        #break

    status = wm.get_submission_status()
    status.to_csv("job_status.csv")


if __name__ == "__main__":
    main()

