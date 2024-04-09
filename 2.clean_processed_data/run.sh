#!/bin/bash

clean_alignment_data() {
    mkdir -p data/minigraph/chm13graph
    mkdir -p data/minigraph/chm13linear
    mkdir -p data/minimap2/chm13
    mkdir -p data/minimap2/hg002

    # HG002 aligned to hg002
    # gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/71d5381d-1de9-4887-8cb5-003f7998990a/LineargenomeAlignment/3e89517a-e6fe-48ab-a540-662bc6ae0510/call-minimapTask/HG002_PACBIO_REVIO.paf data/minimap2/hg002
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/b103af91-ee30-4ded-9cd3-b428e8185313/LineargenomeAlignment/0fd735eb-ca04-458f-a48d-12745fd0042a/call-minimapTask/attempt-2/COLO829.paf data/minimap2/hg002
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/b103af91-ee30-4ded-9cd3-b428e8185313/LineargenomeAlignment/3f229caa-2ef1-4e74-980a-ff8b684b959f/call-minimapTask/COLO829-BL.paf data/minimap2/hg002

    #for paf in gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/b103af91-ee30-4ded-9cd3-b428e8185313/LineargenomeAlignment/0fd735eb-ca04-458f-a48d-12745fd0042a/call-minimapTask/attempt-2/COLO829.paf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/b103af91-ee30-4ded-9cd3-b428e8185313/LineargenomeAlignment/3f229caa-2ef1-4e74-980a-ff8b684b959f/call-minimapTask/COLO829-BL.paf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/b103af91-ee30-4ded-9cd3-b428e8185313/LineargenomeAlignment/85f933f7-cd23-432e-8c61-ed9f070283a3/call-minimapTask/HCC1395.paf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/b103af91-ee30-4ded-9cd3-b428e8185313/LineargenomeAlignment/3c02e5bc-bd4f-43f7-9ebc-d10d947c01f3/call-minimapTask/HCC1395-BL.paf; do
    #    gsutil -m cp $paf data/minimap2/hg002
    #done

    # chm13 linear minimap2
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/1c46d4f2-b4a3-4c01-8d27-0381d4ddb008/LineargenomeAlignment/e0866200-76e7-45a0-896c-71b9e28411df/call-minimapTask/COLO829.paf data/minimap2/chm13
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/99435f4f-ae67-4823-a734-f9babc979669/LineargenomeAlignment/39f23cb3-e819-407c-9927-1ede6eee9b12/call-minimapTask/COLO829_ONT.paf data/minimap2/chm13
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/1c46d4f2-b4a3-4c01-8d27-0381d4ddb008/LineargenomeAlignment/09cad466-4206-4113-94d4-116aec698cc8/call-minimapTask/HCC1395.paf data/minimap2/chm13
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/1c46d4f2-b4a3-4c01-8d27-0381d4ddb008/LineargenomeAlignment/19cb4750-a9f4-4425-9f16-8751e7f724fa/call-minimapTask/COLO829-BL.paf data/minimap2/chm13
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/99435f4f-ae67-4823-a734-f9babc979669/LineargenomeAlignment/489422f9-d954-4329-b4d1-f3dc4ff38316/call-minimapTask/attempt-4/COLO829BL_ONT.paf data/minimap2/chm13
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/1c46d4f2-b4a3-4c01-8d27-0381d4ddb008/LineargenomeAlignment/5542935a-1761-4e02-afdd-7c26bf7a3ca2/call-minimapTask/HCC1395-BL.paf data/minimap2/chm13
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/99435f4f-ae67-4823-a734-f9babc979669/LineargenomeAlignment/a84285e0-ffc8-46cb-9374-7cd59200a9ac/call-minimapTask/attempt-4/HG002_ONT_sup.paf data/minimap2/chm13
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/b5b1b507-867a-437b-b320-e6083f034725/LineargenomeAlignment/fc95b387-5a60-4ef4-8517-46b46f20aec7/call-minimapTask/attempt-2/HG002_PACBIO_REVIO.paf data/minimap2/chm13

    # chm13 graph minigraph: most important alignment triple checks
    #mv ../../phaseA3_L1_TSD_polyA_SVA/minigraph_chm13_graph/*gaf data/minigraph/chm13graph
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/d1aff087-1873-4e94-b148-77945403f2a2/PangenomeAlignment/2ebdab8b-d82e-43be-8507-24dc7926c847/call-minigraphTask/HG002_PACBIO_REVIO.gaf data/minigraph/chm13graph
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/f66ecd5a-a9f0-4a18-9260-fb1f8a1b8bfe/PangenomeAlignment/75171014-f801-4207-8b3a-90ee10b5c0bc/call-minigraphTask/attempt-4/HG002_ONT_sup.gaf data/minigraph/chm13graph
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/f66ecd5a-a9f0-4a18-9260-fb1f8a1b8bfe/PangenomeAlignment/f20e0964-5400-481f-b97f-4c7010a8f706/call-minigraphTask/attempt-4/COLO829_ONT.gaf data/minigraph/chm13graph
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/f66ecd5a-a9f0-4a18-9260-fb1f8a1b8bfe/PangenomeAlignment/d34ff229-54de-402b-8d7d-1fe00ce49beb/call-minigraphTask/COLO829BL_ONT.gaf data/minigraph/chm13graph

    # chm13 linear minigraph
    #mv ../../phaseA3_L1_TSD_polyA_SVA/minigraph_chm13v2_linear/*gaf data/minigraph/chm13linear
    #gsutil cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/78e13b29-e087-4fc7-8dcf-8e71d2b71935/PangenomeAlignment/\*/call-minigraphTask/\*.gaf data/minigraph/chm13linear
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/39030d53-2a6a-4a47-914b-733e66074512/PangenomeAlignment/299a0846-3c64-4fdf-b784-f366eb21026d/call-minigraphTask/HG002_PACBIO_REVIO.gaf data/minigraph/chm13linear
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/5494a02f-ebb1-4522-aba6-6cae7f4938fe/PangenomeAlignment/ce66b6f2-2001-4161-989d-489a0fe11b64/call-minigraphTask/attempt-2/HG002_ONT_sup.gaf data/minigraph/chm13linear
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/5494a02f-ebb1-4522-aba6-6cae7f4938fe/PangenomeAlignment/5e8990ae-e989-413b-88fd-047b7c765cfb/call-minigraphTask/attempt-2/COLO829BL_ONT.gaf data/minigraph/chm13linear
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/5494a02f-ebb1-4522-aba6-6cae7f4938fe/PangenomeAlignment/ab73d380-e61c-42a5-97af-f855a44edae7/call-minigraphTask/COLO829_ONT.gaf data/minigraph/chm13linear

    mkdir -p data/minigraph/grch38graph
    mkdir -p data/minigraph/grch38linear
    mkdir -p data/minimap2/grch38

    gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/d6f44772-38d8-4a28-95b4-851472c1662f/LineargenomeAlignment/6fcd2814-a615-42df-b465-83300b6a1c71/call-minimapTask/attempt-2/HG002_PACBIO_REVIO.paf data/minimap2/grch38

    #mv ../../phaseA3_L1_TSD_polyA_SVA/minigraph_grch38_graph/*gaf data/minigraph/grch38graph
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/0af12b88-d354-4448-81aa-a49078569d59/PangenomeAlignment/\*/call-minigraphTask/\*.gaf data/minigraph/grch38graph
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/0af12b88-d354-4448-81aa-a49078569d59/PangenomeAlignment/\*/call-minigraphTask/\*/\*.gaf data/minigraph/grch38graph
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/48bbec13-bced-4532-a0ba-371d5fd46c52/PangenomeAlignment/e7880bfe-9db4-44ec-9718-a09b4bd3b1e8/call-minigraphTask/attempt-3/HG002_PACBIO_REVIO.gaf data/minigraph/grch38graph
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/3a622caf-59a8-4fc6-8ee8-901f42038b86/PangenomeAlignment/339d35bb-de52-40eb-b366-04369bc38487/call-minigraphTask/attempt-4/HG002_ONT_sup.gaf data/minigraph/grch38graph

    #for gaf in gs://broad_pangenome_sv_alvin/minigraph_GRCh38_graph_add_corrected_ds_Z/COLO829.gaf gs://broad_pangenome_sv_alvin/minigraph_GRCh38_graph_add_corrected_ds_Z/HCC1395.gaf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/3a622caf-59a8-4fc6-8ee8-901f42038b86/PangenomeAlignment/6f07348b-bced-48f5-9424-eaff0e11fd88/call-minigraphTask/attempt-4/COLO829_ONT.gaf; do
    #for gaf in gs://broad_pangenome_sv_alvin/minigraph_GRCh38_graph_add_corrected_ds_Z/COLO829-BL.gaf gs://broad_pangenome_sv_alvin/minigraph_GRCh38_graph_add_corrected_ds_Z/HCC1395-BL.gaf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/3a622caf-59a8-4fc6-8ee8-901f42038b86/PangenomeAlignment/30382589-e524-4fd9-afe1-c3173add872b/call-minigraphTask/COLO829BL_ONT.gaf; do
    #	gsutil -m cp $gaf data/minigraph/grch38graph
    #done

    ##for gaf in gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/d03abbfa-a842-49cb-93bf-d738cdb5ea4d/PangenomeAlignment/28145827-64a0-480a-a7e2-54b9843c23c1/call-minigraphTask/attempt-4/COLO829_ONT.gaf; do
    #for gaf in gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/d03abbfa-a842-49cb-93bf-d738cdb5ea4d/PangenomeAlignment/f4850519-b9d0-477c-9f53-12698212de69/call-minigraphTask/attempt-3/COLO829BL_ONT.gaf; do
    #    gsutil -m cp $gaf data/minigraph/grch38linear
    #done

    ##gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/2c5d0d65-d575-4491-a311-2cf9a1fd695f/PangenomeAlignment/\*/call-minigraphTask/\*.gaf data/minigraph/grch38linear
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/f9c8deb0-66f4-4745-ae68-ca1f1fbd431c/PangenomeAlignment/5bd6b167-e107-4395-9235-fbf447bac80f/call-minigraphTask/attempt-4/HG002_PACBIO_REVIO.gaf data/minigraph/grch38linear
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/d03abbfa-a842-49cb-93bf-d738cdb5ea4d/PangenomeAlignment/336e2b35-c784-42c5-ae1c-9e432d40e4ce/call-minigraphTask/attempt-3/HG002_ONT_sup.gaf data/minigraph/grch38linear

    #for gaf in gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/a788a6e1-7696-4b0e-bcb7-ecadb7e57741/LineargenomeAlignment/a04216ad-e4fe-4beb-b81f-ab74f4512812/call-minimapTask/attempt-3/COLO829.paf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/a788a6e1-7696-4b0e-bcb7-ecadb7e57741/LineargenomeAlignment/62378f57-1603-48b1-98bd-5a028e4b8c23/call-minimapTask/HCC1395.paf; do
    #for gaf in gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/a788a6e1-7696-4b0e-bcb7-ecadb7e57741/LineargenomeAlignment/87c9bf55-4e75-4784-be31-f8704565b5cb/call-minimapTask/COLO829-BL.paf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/a788a6e1-7696-4b0e-bcb7-ecadb7e57741/LineargenomeAlignment/e73a79c2-3ee8-4abf-b4ab-dd431e287927/call-minimapTask/attempt-4/HCC1395-BL.paf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/5309edd4-995d-4675-9237-75e45da93548/LineargenomeAlignment/45901e47-b619-4637-80e5-3beb88e91fb0/call-minimapTask/attempt-4/COLO829_ONT.paf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/5309edd4-995d-4675-9237-75e45da93548/LineargenomeAlignment/6a24537b-c793-4ebe-b252-e15c07ae901a/call-minimapTask/attempt-4/COLO829BL_ONT.paf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/5309edd4-995d-4675-9237-75e45da93548/LineargenomeAlignment/e047b882-f43f-43df-bd13-123677afff6a/call-minimapTask/attempt-4/HG002_ONT_sup.paf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/d6f44772-38d8-4a28-95b4-851472c1662f/LineargenomeAlignment/6fcd2814-a615-42df-b465-83300b6a1c71/call-minimapTask/attempt-2/HG002_PACBIO_REVIO.paf;do
    #    gsutil -m cp $gaf data/minimap2/grch38
    #done

    mkdir -p data/minigraph/grch37graph
    mkdir -p data/minigraph/grch37linear
    mkdir -p data/minimap2/grch37

    #gsutil cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/5408b7e3-a5b9-4ece-9551-1cfe66ae0da6/LineargenomeAlignment/02a3ce45-4d99-495c-9a04-920241909358/call-minimapTask/attempt-4/HG002_PACBIO_REVIO.paf data/minimap2/grch37

    #mv ../../phaseA3_L1_TSD_polyA_SVA/minigraph_grch37_graph/*gaf data/minigraph/grch37graph
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/be96f542-c25d-4155-9b44-905e018ee705/PangenomeAlignment/8dfd6372-088c-4c00-b557-6562e41771bc/call-minigraphTask/COLO829.gaf .
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/be96f542-c25d-4155-9b44-905e018ee705/PangenomeAlignment/6a337792-fa50-4522-83de-011632a6175d/call-minigraphTask/COLO829-BL.gaf .
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/4669590c-9fb2-43c2-b1fe-01bfe25156f0/PangenomeAlignment/7c3d6447-4eaf-41e7-9413-91961e666ddf/call-minigraphTask/COLO829BL_ONT.gaf .
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/4669590c-9fb2-43c2-b1fe-01bfe25156f0/PangenomeAlignment/07ae938e-beb8-4109-a6ab-27c086a55d9f/call-minigraphTask/attempt-2/COLO829_ONT.gaf .
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/be96f542-c25d-4155-9b44-905e018ee705/PangenomeAlignment/56199a96-db98-4eaa-a048-d9ead6d3dac0/call-minigraphTask/attempt-3/HCC1395.gaf .
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/be96f542-c25d-4155-9b44-905e018ee705/PangenomeAlignment/f83e59bc-afab-4477-b4a1-be165653d497/call-minigraphTask/attempt-2/HCC1395-BL.gaf .
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/0ea4741c-eb60-46be-911c-8f3813336647/PangenomeAlignment/e88b381a-abec-4daf-9a65-dbb94e460bf9/call-minigraphTask/attempt-3/HG002_PACBIO_REVIO.gaf .
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/4669590c-9fb2-43c2-b1fe-01bfe25156f0/PangenomeAlignment/a96f9741-849d-427b-8270-b986a0cc1462/call-minigraphTask/attempt-2/HG002_ONT_sup.gaf .
    #for gaf in gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/be96f542-c25d-4155-9b44-905e018ee705/PangenomeAlignment/8dfd6372-088c-4c00-b557-6562e41771bc/call-minigraphTask/COLO829.gaf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/4669590c-9fb2-43c2-b1fe-01bfe25156f0/PangenomeAlignment/07ae938e-beb8-4109-a6ab-27c086a55d9f/call-minigraphTask/attempt-2/COLO829_ONT.gaf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/be96f542-c25d-4155-9b44-905e018ee705/PangenomeAlignment/56199a96-db98-4eaa-a048-d9ead6d3dac0/call-minigraphTask/attempt-3/HCC1395.gaf; do
    #    gsutil -m cp $gaf data/minigraph/grch37graph
    #done

    #mv ../../phaseA3_L1_TSD_polyA_SVA/minigraph_grch37_linear/*gaf data/minigraph/grch37linear
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/0fc46dfe-068c-485e-85bb-ff9de9f5de6c/PangenomeAlignment/8d69ed77-988c-483d-bb6a-eed4ab733f20/call-minigraphTask/HG002_ONT_sup.gaf .
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/e5a3c44b-65cd-4f87-9b54-8a9a2053b35f/PangenomeAlignment/78aa1bee-204d-476b-a03e-574c0880dc03/call-minigraphTask/HG002_PACBIO_REVIO.gaf .
    #for gaf in gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/e5a3c44b-65cd-4f87-9b54-8a9a2053b35f/PangenomeAlignment/ae7a4b21-0ca6-4e75-9e19-31cab31aafc3/call-minigraphTask/COLO829.gaf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/0fc46dfe-068c-485e-85bb-ff9de9f5de6c/PangenomeAlignment/87993fc5-eb48-41de-8d58-ae86e773dd9e/call-minigraphTask/attempt-2/COLO829_ONT.gaf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/e5a3c44b-65cd-4f87-9b54-8a9a2053b35f/PangenomeAlignment/b5727ffa-d367-49ac-8e2e-9f3387d1ac0b/call-minigraphTask/HCC1395.gaf; do
    #for gaf in gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/e5a3c44b-65cd-4f87-9b54-8a9a2053b35f/PangenomeAlignment/26a916b6-4359-43ed-939d-4733ebea2497/call-minigraphTask/COLO829-BL.gaf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/0fc46dfe-068c-485e-85bb-ff9de9f5de6c/PangenomeAlignment/b1941c6f-ef12-4c7d-8af7-72fe558a4ccf/call-minigraphTask/COLO829BL_ONT.gaf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/e5a3c44b-65cd-4f87-9b54-8a9a2053b35f/PangenomeAlignment/10000eec-5bf7-48a2-b8fe-464ee47439c4/call-minigraphTask/HCC1395-BL.gaf;do
    #    gsutil -m cp $gaf data/minigraph/grch37linear &
    #done

    #for paf in gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/5c90b169-2fcf-49c6-99b7-256584a425c4/LineargenomeAlignment/9a3e3256-758f-4eab-ac95-6da910c04be8/call-minimapTask/COLO829.paf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/5c90b169-2fcf-49c6-99b7-256584a425c4/LineargenomeAlignment/a3ed5ed3-8a4c-427f-b34c-41ab12267a18/call-minimapTask/HCC1395.paf; do
    #for paf in gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/5c90b169-2fcf-49c6-99b7-256584a425c4/LineargenomeAlignment/054be214-9e35-4bb5-9af0-488476feaa3b/call-minimapTask/COLO829-BL.paf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/5c90b169-2fcf-49c6-99b7-256584a425c4/LineargenomeAlignment/f7202af5-3ccc-4e30-bd75-a300d75976a3/call-minimapTask/HCC1395-BL.paf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/9a15f177-5ac0-4f94-8a7e-1430ef77a554/LineargenomeAlignment/288e8e07-211c-4cad-a64f-46d3e38e506d/call-minimapTask/attempt-4/COLO829BL_ONT.paf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/9a15f177-5ac0-4f94-8a7e-1430ef77a554/LineargenomeAlignment/71fe178c-5465-4272-ac7a-fe10bd3d7806/call-minimapTask/attempt-4/COLO829_ONT.paf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/9a15f177-5ac0-4f94-8a7e-1430ef77a554/LineargenomeAlignment/759b3694-f4bf-455e-8cbb-731a135f7916/call-minimapTask/HG002_ONT_sup.paf gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/b5dab3de-e2de-4324-beb2-b87945cd9bc4/LineargenomeAlignment/a268ccae-8596-4939-988f-e1602bcf5db4/call-minimapTask/HG002_PACBIO_REVIO.paf; do
    #    gsutil -m cp $paf data/minimap2/grch37
    #done
}

clean_sniffles() {
    mkdir -p data/sniffles2/grch37
    mkdir -p data/sniffles2/grch38
    mkdir -p data/sniffles2/chm13

    #tumor only mode
    #for vcf in gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/3a88f5a3-929f-439b-a5eb-2e43e3135e0d/SNFWorkflow/70ab025d-86cb-45b6-82d5-249a753f9136/call-SNFTask/COLO829.vcf.gz\* gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/3a88f5a3-929f-439b-a5eb-2e43e3135e0d/SNFWorkflow/8b548aea-be9b-40ea-a18d-eed14d415c1e/call-SNFTask/COLO829_ONT.vcf.gz\* gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/3a88f5a3-929f-439b-a5eb-2e43e3135e0d/SNFWorkflow/6ae133a0-5f03-46a4-b21f-432f70ee0eb7/call-SNFTask/HCC1395.vcf.gz\*; do
    #    gsutil -m cp $vcf data/sniffles2/grch37
    #done

    ##hg002
    #for vcf in gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/3a88f5a3-929f-439b-a5eb-2e43e3135e0d/SNFWorkflow/8039e150-df3b-4391-9267-d778800a150d/call-SNFTask/attempt-2/HG002_ONT_sup.vcf.gz gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/3a88f5a3-929f-439b-a5eb-2e43e3135e0d/SNFWorkflow/8039e150-df3b-4391-9267-d778800a150d/call-SNFTask/attempt-2/HG002_ONT_sup.vcf.gz.tbi gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/3a88f5a3-929f-439b-a5eb-2e43e3135e0d/SNFWorkflow/e1a2ad0d-7e13-4f4e-b297-fb9adebe61f6/call-SNFTask/HG002_PACBIO_REVIO.vcf.gz.tbi gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/3a88f5a3-929f-439b-a5eb-2e43e3135e0d/SNFWorkflow/e1a2ad0d-7e13-4f4e-b297-fb9adebe61f6/call-SNFTask/HG002_PACBIO_REVIO.vcf.gz; do
    #    gsutil -m cp $vcf data/sniffles2/grch37
    #done

    #for vcf in gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/27981bb8-033a-456f-96c0-133c5a32f5dc/SNFWorkflow/871f2809-41d5-448a-b56d-5f5eb9a6a008/call-SNFTask/COLO829.vcf.gz\* gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/4ad27225-2824-47b5-a2a6-899185b3c298/SNFWorkflow/4a032a73-c785-4259-be3c-f6f8eaaff3ce/call-SNFTask/attempt-2/COLO829_ONT.vcf.gz\* gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/27981bb8-033a-456f-96c0-133c5a32f5dc/SNFWorkflow/e1843bed-2711-45c0-9931-daef90aa1d4d/call-SNFTask/HCC1395.vcf.gz\*; do
    #    gsutil -m cp $vcf data/sniffles2/grch38
    #done

    #hg002
    #for vcf in gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/4ad27225-2824-47b5-a2a6-899185b3c298/SNFWorkflow/2347cd57-0ddc-4b93-859f-87f045d505f6/call-SNFTask/HG002_ONT_sup.vcf.gz.tbi gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/4ad27225-2824-47b5-a2a6-899185b3c298/SNFWorkflow/2347cd57-0ddc-4b93-859f-87f045d505f6/call-SNFTask/HG002_ONT_sup.vcf.gz gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/27981bb8-033a-456f-96c0-133c5a32f5dc/SNFWorkflow/9e6bac59-88f0-4ca4-a5ce-793603ede424/call-SNFTask/HG002_PACBIO_REVIO.vcf.gz.tbi gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/27981bb8-033a-456f-96c0-133c5a32f5dc/SNFWorkflow/9e6bac59-88f0-4ca4-a5ce-793603ede424/call-SNFTask/HG002_PACBIO_REVIO.vcf.gz; do
    #    gsutil -m cp $vcf data/sniffles2/grch38
    #done

    #hg002
    for vcf in gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/e64ded27-805b-443a-84db-1541c99c1dcf/SNFWorkflow/36eaf157-ab9c-4c4d-a977-cc0690297188/call-SNFTask/COLO829.vcf.gz gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/e64ded27-805b-443a-84db-1541c99c1dcf/SNFWorkflow/36eaf157-ab9c-4c4d-a977-cc0690297188/call-SNFTask/COLO829.vcf.gz.tbi gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/e64ded27-805b-443a-84db-1541c99c1dcf/SNFWorkflow/6c876fc0-6971-4347-b6ca-874685f64f9f/call-SNFTask/COLO829_ONT.vcf.gz gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/e64ded27-805b-443a-84db-1541c99c1dcf/SNFWorkflow/6c876fc0-6971-4347-b6ca-874685f64f9f/call-SNFTask/COLO829_ONT.vcf.gz.tbi gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/e64ded27-805b-443a-84db-1541c99c1dcf/SNFWorkflow/c6513cda-9366-4f0b-95e5-cc65af4655cb/call-SNFTask/attempt-2/HCC1395.vcf.gz gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/e64ded27-805b-443a-84db-1541c99c1dcf/SNFWorkflow/c6513cda-9366-4f0b-95e5-cc65af4655cb/call-SNFTask/attempt-2/HCC1395.vcf.gz.tbi; do
	   gsutil -m cp $vcf data/sniffles2/chm13
    done

    #for vcf in gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/e64ded27-805b-443a-84db-1541c99c1dcf/SNFWorkflow/5432c093-83d9-443b-b91d-ddfebf76fdd7/call-SNFTask/HG002_ONT_sup.vcf.gz gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/e64ded27-805b-443a-84db-1541c99c1dcf/SNFWorkflow/5432c093-83d9-443b-b91d-ddfebf76fdd7/call-SNFTask/HG002_ONT_sup.vcf.gz.tbi gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/e64ded27-805b-443a-84db-1541c99c1dcf/SNFWorkflow/e62925ea-3b69-499f-8401-ee3c16e5f569/call-SNFTask/HG002_PACBIO_REVIO.vcf.gz gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/e64ded27-805b-443a-84db-1541c99c1dcf/SNFWorkflow/e62925ea-3b69-499f-8401-ee3c16e5f569/call-SNFTask/HG002_PACBIO_REVIO.vcf.gz.tbi; do
    #    gsutil -m cp $vcf data/sniffles2/chm13
    #done
}

clean_severus() {
    mkdir -p data/severus/grch37
    mkdir -p data/severus/grch38
    mkdir -p data/severus/chm13

    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/728dae11-5ac3-49c8-820a-79f40e502be9/SeverusWorkflow/254b4ab3-03d6-45f2-8b1f-94ba34fedfd9/call-severusTumor/attempt-2/glob-443a79c6c426960de8526eea118fdfe5/severus_all.vcf data/severus/grch37/COLO829.vcf
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/728dae11-5ac3-49c8-820a-79f40e502be9/SeverusWorkflow/d6aefc7a-0ef2-43d1-ba22-730fb4df486a/call-severusTumor/attempt-2/glob-966d0305651d70bfbb44ea0eb6b76fc4/severus_all.vcf data/severus/grch37/COLO829_ONT.vcf
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/728dae11-5ac3-49c8-820a-79f40e502be9/SeverusWorkflow/2b0beb1e-9482-41ac-b3aa-a977a61cc60b/call-severusTumor/glob-c54210335a4577a1db5ddcf420283162/severus_all.vcf data/severus/grch37/HCC1395.vcf
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/728dae11-5ac3-49c8-820a-79f40e502be9/SeverusWorkflow/43adc864-544a-42f9-9aff-68a7c1ed2bf5/call-severusTumor/glob-44bab286e11b18dfd5a120503d966eba/severus_all.vcf data/severus/grch37/HG002_ONT.vcf
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/728dae11-5ac3-49c8-820a-79f40e502be9/SeverusWorkflow/852f85c8-8788-4ef9-8819-68841bc8ebb8/call-severusTumor/attempt-4/glob-0aacef18f32e6d689d769b045b20774a/severus_all.vcf data/severus/grch37/HG002_PACBIO_REVIO.vcf

    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/feb01427-357a-4f8b-9ab2-47789d75b054/SeverusWorkflow/4d5dda7f-2d54-47e8-aee8-f6f96fd802f2/call-severusTumor/attempt-2/glob-443a79c6c426960de8526eea118fdfe5/severus_all.vcf data/severus/grch38/COLO829.vcf
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/e9d77840-2603-40fd-94a5-d01575046496/SeverusWorkflow/884e8dc5-fdcf-403b-81a7-6e9d898a2101/call-severusTumor/attempt-4/glob-966d0305651d70bfbb44ea0eb6b76fc4/severus_all.vcf data/severus/grch38/COLO829_ONT.vcf
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/feb01427-357a-4f8b-9ab2-47789d75b054/SeverusWorkflow/7b3d881e-ad6a-4284-a6bc-f342ff8c4daf/call-severusTumor/attempt-4/glob-c54210335a4577a1db5ddcf420283162/severus_all.vcf data/severus/grch38/HCC1395.vcf
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/e9d77840-2603-40fd-94a5-d01575046496/SeverusWorkflow/65049b5e-f2f3-47a9-abc9-c0e16be58759/call-severusTumor/glob-44bab286e11b18dfd5a120503d966eba/severus_all.vcf data/severus/grch38/HG002_ONT.vcf
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/feb01427-357a-4f8b-9ab2-47789d75b054/SeverusWorkflow/4cbf4628-49ed-416f-a808-61f57d15ef10/call-severusTumor/attempt-4/glob-0aacef18f32e6d689d769b045b20774a/severus_all.vcf data/severus/grch38/HG002_PACBIO_REVIO.vcf

    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/222d23a0-8298-48b7-8481-9d7260350604/SeverusWorkflow/c427af73-9350-480d-96e3-f9303ba8474d/call-severusTumor/glob-443a79c6c426960de8526eea118fdfe5/severus_all.vcf data/severus/chm13/COLO829.vcf
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/222d23a0-8298-48b7-8481-9d7260350604/SeverusWorkflow/00b0043c-c9c6-4ba6-b027-97289be041e5/call-severusTumor/glob-966d0305651d70bfbb44ea0eb6b76fc4/severus_all.vcf data/severus/chm13/COLO829_ONT.vcf
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/222d23a0-8298-48b7-8481-9d7260350604/SeverusWorkflow/e9d6144e-7a11-4fd3-b20d-793e620063a2/call-severusTumor/glob-c54210335a4577a1db5ddcf420283162/severus_all.vcf data/severus/chm13/HCC1395.vcf
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/222d23a0-8298-48b7-8481-9d7260350604/SeverusWorkflow/77e31fbf-cd42-45cd-a980-bd898a561db8/call-severusTumor/attempt-3/glob-44bab286e11b18dfd5a120503d966eba/severus_all.vcf data/severus/chm13/HG002_ONT.vcf
    #gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/222d23a0-8298-48b7-8481-9d7260350604/SeverusWorkflow/d1ab1b81-bd30-4867-9a91-8d80423f7651/call-severusTumor/glob-0aacef18f32e6d689d769b045b20774a/severus_all.vcf data/severus/chm13/HG002_PACBIO_REVIO.vcf

    #tumor normal pair
    mkdir -p data/severus/grch37_pair
    mkdir -p data/severus/grch38_pair
    mkdir -p data/severus/chm13_pair

    gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/7db51db5-5aeb-4acb-bc6f-9034f2a2b97b/SeverusWorkflow/f4dc0e0c-3856-4bd9-8a10-c1bcdc1afe0b/call-severusTumorNormal/attempt-2/glob-b9e1d563654cd7134d51f5a06670f7e5/severus_somatic.vcf data/severus/chm13_pair/COLO829_pair.vcf
    gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/7db51db5-5aeb-4acb-bc6f-9034f2a2b97b/SeverusWorkflow/ebf8db9e-2f4b-4dff-813c-0015e903782a/call-severusTumorNormal/attempt-3/glob-6e0e890dbf8220d68888f7613d6f747e/severus_somatic.vcf data/severus/chm13_pair/COLO829_ONT_pair.vcf
    gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/7db51db5-5aeb-4acb-bc6f-9034f2a2b97b/SeverusWorkflow/4a6b65be-1ae1-44de-b12e-9ec5a147c954/call-severusTumorNormal/glob-78df5ddc04f04cc6bc7990afe3e0604e/severus_somatic.vcf data/severus/chm13_pair/HCC1395_pair.vcf

    gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/595cdb4d-bfcc-475e-838f-4c98d00957d3/SeverusWorkflow/e6897378-3792-4960-bd4a-5b9d4d30dcc8/call-severusTumorNormal/attempt-2/glob-b9e1d563654cd7134d51f5a06670f7e5/severus_somatic.vcf data/severus/grch37_pair/COLO829_pair.vcf
    gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/595cdb4d-bfcc-475e-838f-4c98d00957d3/SeverusWorkflow/0627a1ea-0764-4379-861f-f7a2056570e6/call-severusTumorNormal/attempt-2/glob-6e0e890dbf8220d68888f7613d6f747e/severus_somatic.vcf data/severus/grch37_pair/COLO829_ONT_pair.vcf
    gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/595cdb4d-bfcc-475e-838f-4c98d00957d3/SeverusWorkflow/2937053c-71ae-4483-bbb6-263ef8b5f833/call-severusTumorNormal/attempt-2/glob-78df5ddc04f04cc6bc7990afe3e0604e/severus_somatic.vcf data/severus/grch37_pair/HCC1395_pair.vcf

    gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/f11d6d7f-9bae-47fc-8dfd-9231fb2ff86d/SeverusWorkflow/6abd87a5-d269-4df0-a15d-7e87ec830b57/call-severusTumorNormal/attempt-2/glob-b9e1d563654cd7134d51f5a06670f7e5/severus_somatic.vcf data/severus/grch38_pair/COLO829_pair.vcf
    gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/77a1635f-17f0-4541-8d08-2f86e7c4ed25/SeverusWorkflow/c14b4415-9cab-4785-8433-602816588fc3/call-severusTumorNormal/attempt-4/glob-6e0e890dbf8220d68888f7613d6f747e/severus_somatic.vcf data/severus/grch38_pair/COLO829_ONT_pair.vcf
    gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/f11d6d7f-9bae-47fc-8dfd-9231fb2ff86d/SeverusWorkflow/b5ff4d09-7924-417e-8fc9-46cd562f1c45/call-severusTumorNormal/glob-78df5ddc04f04cc6bc7990afe3e0604e/severus_somatic.vcf data/severus/grch38_pair/HCC1395_pair.vcf
}

add_reference() {
    mkdir -p reference
    gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/human_GRCh38_no_alt_analysis_set.trf.bed .
    gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/human_hs37d5.trf.bed .
    gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/chm13v2.0.trf.bed .
}

main() {
    clean_alignment_data

    #clean_severus
    #clean_sniffles
    #add_reference
}

main

