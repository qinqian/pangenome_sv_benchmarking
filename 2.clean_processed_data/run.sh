#!/bin/bash

clean_alignment_data() {
    mkdir -p data/minigraph/chm13graph
    mkdir -p data/minigraph/chm13linear
    mkdir -p data/minimap2/chm13
    # add tumor chm13 minimap2 here
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/1c46d4f2-b4a3-4c01-8d27-0381d4ddb008/LineargenomeAlignment/e0866200-76e7-45a0-896c-71b9e28411df/call-minimapTask/COLO829.paf data/minimap2/chm13
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/99435f4f-ae67-4823-a734-f9babc979669/LineargenomeAlignment/39f23cb3-e819-407c-9927-1ede6eee9b12/call-minimapTask/COLO829_ONT.paf data/minimap2/chm13
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/1c46d4f2-b4a3-4c01-8d27-0381d4ddb008/LineargenomeAlignment/09cad466-4206-4113-94d4-116aec698cc8/call-minimapTask/HCC1395.paf data/minimap2/chm13

    #mv ../../phaseA3_L1_TSD_polyA_SVA/minigraph_chm13_graph/*gaf data/minigraph/chm13graph
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/d1aff087-1873-4e94-b148-77945403f2a2/PangenomeAlignment/2ebdab8b-d82e-43be-8507-24dc7926c847/call-minigraphTask/HG002_PACBIO_REVIO.gaf .

    #mv ../../phaseA3_L1_TSD_polyA_SVA/minigraph_chm13v2_linear/*gaf data/minigraph/chm13linear
    #gsutil cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/78e13b29-e087-4fc7-8dcf-8e71d2b71935/PangenomeAlignment/\*/call-minigraphTask/\*.gaf data/minigraph/chm13linear
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/39030d53-2a6a-4a47-914b-733e66074512/PangenomeAlignment/299a0846-3c64-4fdf-b784-f366eb21026d/call-minigraphTask/HG002_PACBIO_REVIO.gaf data/minigraph/chm13linear

    #gsutil cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/5494a02f-ebb1-4522-aba6-6cae7f4938fe/PangenomeAlignment/ab73d380-e61c-42a5-97af-f855a44edae7/call-minigraphTask/COLO829_ONT.gaf data/minigraph/chm13linear
    gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/f66ecd5a-a9f0-4a18-9260-fb1f8a1b8bfe/PangenomeAlignment/f20e0964-5400-481f-b97f-4c7010a8f706/call-minigraphTask/attempt-4/COLO829_ONT.gaf data/minigraph/chm13graph

    #mkdir -p data/minigraph/grch38graph
    #mkdir -p data/minigraph/grch38linear
    #mkdir -p data/minimap2/grch38

    #for gaf in gs://broad_pangenome_sv_alvin/minigraph_GRCh38_graph_add_corrected_ds_Z/COLO829.gaf gs://broad_pangenome_sv_alvin/minigraph_GRCh38_graph_add_corrected_ds_Z/HCC1395.gaf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/3a622caf-59a8-4fc6-8ee8-901f42038b86/PangenomeAlignment/6f07348b-bced-48f5-9424-eaff0e11fd88/call-minigraphTask/attempt-4/COLO829_ONT.gaf; do
    #	gsutil -m cp $gaf data/minigraph/grch38graph
    #done

    #for gaf in gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/d03abbfa-a842-49cb-93bf-d738cdb5ea4d/PangenomeAlignment/28145827-64a0-480a-a7e2-54b9843c23c1/call-minigraphTask/attempt-4/COLO829_ONT.gaf; do
    #    gsutil -m cp $gaf data/minigraph/grch38linear
    #done
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/2c5d0d65-d575-4491-a311-2cf9a1fd695f/PangenomeAlignment/\*/call-minigraphTask/\*.gaf data/minigraph/grch38linear
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/f9c8deb0-66f4-4745-ae68-ca1f1fbd431c/PangenomeAlignment/5bd6b167-e107-4395-9235-fbf447bac80f/call-minigraphTask/attempt-4/HG002_PACBIO_REVIO.gaf data/minigraph/grch38linear

    #for gaf in gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/a788a6e1-7696-4b0e-bcb7-ecadb7e57741/LineargenomeAlignment/a04216ad-e4fe-4beb-b81f-ab74f4512812/call-minimapTask/attempt-3/COLO829.paf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/a788a6e1-7696-4b0e-bcb7-ecadb7e57741/LineargenomeAlignment/62378f57-1603-48b1-98bd-5a028e4b8c23/call-minimapTask/HCC1395.paf; do
    #    gsutil -m cp $gaf data/minimap2/grch38
    #done

    mkdir -p data/minigraph/grch37graph
    mkdir -p data/minigraph/grch37linear
    mkdir -p data/minimap2/grch37

    #for gaf in gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/be96f542-c25d-4155-9b44-905e018ee705/PangenomeAlignment/8dfd6372-088c-4c00-b557-6562e41771bc/call-minigraphTask/COLO829.gaf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/4669590c-9fb2-43c2-b1fe-01bfe25156f0/PangenomeAlignment/07ae938e-beb8-4109-a6ab-27c086a55d9f/call-minigraphTask/attempt-2/COLO829_ONT.gaf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/be96f542-c25d-4155-9b44-905e018ee705/PangenomeAlignment/56199a96-db98-4eaa-a048-d9ead6d3dac0/call-minigraphTask/attempt-3/HCC1395.gaf; do
    #    gsutil -m cp $gaf data/minigraph/grch37graph
    #done

    #for gaf in gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/e5a3c44b-65cd-4f87-9b54-8a9a2053b35f/PangenomeAlignment/ae7a4b21-0ca6-4e75-9e19-31cab31aafc3/call-minigraphTask/COLO829.gaf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/0fc46dfe-068c-485e-85bb-ff9de9f5de6c/PangenomeAlignment/87993fc5-eb48-41de-8d58-ae86e773dd9e/call-minigraphTask/attempt-2/COLO829_ONT.gaf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/e5a3c44b-65cd-4f87-9b54-8a9a2053b35f/PangenomeAlignment/b5727ffa-d367-49ac-8e2e-9f3387d1ac0b/call-minigraphTask/HCC1395.gaf; do
    #    gsutil -m cp $gaf data/minigraph/grch37linear
    #done

    #for paf in gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/5c90b169-2fcf-49c6-99b7-256584a425c4/LineargenomeAlignment/9a3e3256-758f-4eab-ac95-6da910c04be8/call-minimapTask/COLO829.paf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/5c90b169-2fcf-49c6-99b7-256584a425c4/LineargenomeAlignment/a3ed5ed3-8a4c-427f-b34c-41ab12267a18/call-minimapTask/HCC1395.paf; do
    #    gsutil -m cp $paf data/minimap2/grch37
    #done
}

clean_sniffles() {
    mkdir -p data/sniffles2/grch37
    mkdir -p data/sniffles2/grch38

    #for vcf in gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/3a88f5a3-929f-439b-a5eb-2e43e3135e0d/SNFWorkflow/70ab025d-86cb-45b6-82d5-249a753f9136/call-SNFTask/COLO829.vcf.gz\* gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/3a88f5a3-929f-439b-a5eb-2e43e3135e0d/SNFWorkflow/8b548aea-be9b-40ea-a18d-eed14d415c1e/call-SNFTask/COLO829_ONT.vcf.gz\* gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/3a88f5a3-929f-439b-a5eb-2e43e3135e0d/SNFWorkflow/6ae133a0-5f03-46a4-b21f-432f70ee0eb7/call-SNFTask/HCC1395.vcf.gz\*; do
    #    gsutil -m cp $vcf data/sniffles2/grch37
    #done

    #for vcf in gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/27981bb8-033a-456f-96c0-133c5a32f5dc/SNFWorkflow/871f2809-41d5-448a-b56d-5f5eb9a6a008/call-SNFTask/COLO829.vcf.gz\* gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/4ad27225-2824-47b5-a2a6-899185b3c298/SNFWorkflow/4a032a73-c785-4259-be3c-f6f8eaaff3ce/call-SNFTask/attempt-2/COLO829_ONT.vcf.gz\* gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/27981bb8-033a-456f-96c0-133c5a32f5dc/SNFWorkflow/e1843bed-2711-45c0-9931-daef90aa1d4d/call-SNFTask/HCC1395.vcf.gz\*; do
    #    gsutil -m cp $vcf data/sniffles2/grch38
    #done
}

clean_severus() {
    mkdir -p data/severus/grch37
    mkdir -p data/severus/grch38

    gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/728dae11-5ac3-49c8-820a-79f40e502be9/SeverusWorkflow/254b4ab3-03d6-45f2-8b1f-94ba34fedfd9/call-severusTumor/attempt-2/glob-443a79c6c426960de8526eea118fdfe5/severus_all.vcf data/severus/grch37/COLO829.vcf
    gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/728dae11-5ac3-49c8-820a-79f40e502be9/SeverusWorkflow/d6aefc7a-0ef2-43d1-ba22-730fb4df486a/call-severusTumor/attempt-2/glob-966d0305651d70bfbb44ea0eb6b76fc4/severus_all.vcf data/severus/grch37/COLO829_ONT.vcf
    gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/728dae11-5ac3-49c8-820a-79f40e502be9/SeverusWorkflow/2b0beb1e-9482-41ac-b3aa-a977a61cc60b/call-severusTumor/glob-c54210335a4577a1db5ddcf420283162/severus_all.vcf data/severus/grch37/HCC1395.vcf

    gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/feb01427-357a-4f8b-9ab2-47789d75b054/SeverusWorkflow/4d5dda7f-2d54-47e8-aee8-f6f96fd802f2/call-severusTumor/attempt-2/glob-443a79c6c426960de8526eea118fdfe5/severus_all.vcf data/severus/grch38/COLO829.vcf
    gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/e9d77840-2603-40fd-94a5-d01575046496/SeverusWorkflow/884e8dc5-fdcf-403b-81a7-6e9d898a2101/call-severusTumor/attempt-4/glob-966d0305651d70bfbb44ea0eb6b76fc4/severus_all.vcf data/severus/grch38/COLO829_ONT.vcf
    gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/feb01427-357a-4f8b-9ab2-47789d75b054/SeverusWorkflow/7b3d881e-ad6a-4284-a6bc-f342ff8c4daf/call-severusTumor/attempt-4/glob-c54210335a4577a1db5ddcf420283162/severus_all.vcf data/severus/grch38/HCC1395.vcf

}

add_reference() {
    mkdir -p reference
    gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/human_GRCh38_no_alt_analysis_set.trf.bed .
    gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/human_hs37d5.trf.bed .
    gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/chm13v2.0.trf.bed .
}

main() {
    clean_alignment_data
    #clean_sniffles
    #clean_severus
    #add_reference
}

main

