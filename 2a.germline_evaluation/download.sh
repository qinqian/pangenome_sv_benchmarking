#!/bin/bash -ex

gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/9a5cada4-8841-47a8-bdec-d0c9497e81a8/SeverusWorkflow/7af7bbde-5022-4f65-b794-b49833b8a258/call-severusTumor/glob-0aacef18f32e6d689d769b045b20774a/severus_all.vcf severus_HG002_PACBIO_REVIO.vcf

bcftools sort severus_HG002_PACBIO_REVIO.vcf > severus_HG002_PACBIO_REVIO_sort.vcf
bgzip -f severus_HG002_PACBIO_REVIO_sort.vcf
tabix severus_HG002_PACBIO_REVIO_sort.vcf.gz

#gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/3a88f5a3-929f-439b-a5eb-2e43e3135e0d/SNFWorkflow/e1a2ad0d-7e13-4f4e-b297-fb9adebe61f6/call-SNFTask/HG002_PACBIO_REVIO.vcf.gz sniffles2_HG002_PACBIO_REVIO.vcf.gz
#gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/3a88f5a3-929f-439b-a5eb-2e43e3135e0d/SNFWorkflow/e1a2ad0d-7e13-4f4e-b297-fb9adebe61f6/call-SNFTask/HG002_PACBIO_REVIO.vcf.gz.tbi sniffles2_HG002_PACBIO_REVIO.vcf.gz.tbi

#gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/3a88f5a3-929f-439b-a5eb-2e43e3135e0d/SNFWorkflow/8039e150-df3b-4391-9267-d778800a150d/call-SNFTask/attempt-2/HG002_ONT_sup.vcf.gz sniffles2_HG002_ONT_sup.vcf.gz
#gsutil -m cp gs://fc-secure-062e6633-7a72-4623-9394-491e0ce6c324/submissions/3a88f5a3-929f-439b-a5eb-2e43e3135e0d/SNFWorkflow/8039e150-df3b-4391-9267-d778800a150d/call-SNFTask/attempt-2/HG002_ONT_sup.vcf.gz.tbi sniffles2_HG002_ONT_sup.vcf.gz.tbi
