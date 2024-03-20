Install snakemake
=================

mamba create -n snakemake python=3.11 snakemake 

Install truvari
=================

conda activate snakemake
python3 -m pip install Truvari 


liftover bed file
====================
pip3 install CrossMap 
pip install liftover
# picard will miss some variants when liftover
mamba install  bioconda::picard

