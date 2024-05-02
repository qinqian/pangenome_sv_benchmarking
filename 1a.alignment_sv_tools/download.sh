#wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
#wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz -O grch38.fa.gz # NOTE: need unzip and bgzip -c 
#wget -b -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.0.1.fasta.gz -O hg002v1.0.1.fa.gz

wget -b -c https://zenodo.org/records/6983934/files/chm13-90c.r518.gfa.gz?download=1 -O chm13.gfa.gz
wget -b -c https://zenodo.org/records/6983934/files/GRCh38-90c.r518.gfa.gz?download=1 -O grch38.gfa.gz

