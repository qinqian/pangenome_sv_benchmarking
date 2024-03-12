#!/bin/bash

clean_alignment_data() {
    mkdir -p data/minigraph/chm13graph
    mkdir -p data/minigraph/chm13linear
    mkdir -p data/minimap2/chm13

    mkdir -p data/minigraph/grch38graph
    mkdir -p data/minigraph/grch38linear
    mkdir -p data/minimap2/grch38

    #for gaf in gs://broad_pangenome_sv_alvin/minigraph_GRCh38_graph_add_corrected_ds_Z/COLO829.gaf gs://broad_pangenome_sv_alvin/minigraph_GRCh38_graph_add_corrected_ds_Z/HCC1395.gaf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/3a622caf-59a8-4fc6-8ee8-901f42038b86/PangenomeAlignment/6f07348b-bced-48f5-9424-eaff0e11fd88/call-minigraphTask/attempt-4/COLO829_ONT.gaf; do
    #	gsutil -m cp $gaf data/minigraph/grch38graph
    #done

    #for gaf in gs://broad_pangenome_sv_alvin/minigraph_GCA_000001405.15_GRCh38_no_alt_linear/COLO829.gaf gs://broad_pangenome_sv_alvin/minigraph_GCA_000001405.15_GRCh38_no_alt_linear/HCC1395.gaf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/d03abbfa-a842-49cb-93bf-d738cdb5ea4d/PangenomeAlignment/28145827-64a0-480a-a7e2-54b9843c23c1/call-minigraphTask/attempt-4/COLO829_ONT.gaf; do
    #    gsutil -m cp $gaf data/minigraph/grch38linear
    #done

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

main() {
    #clean_alignment_data
    #clean_sniffles
    #clean_severus
}

main

