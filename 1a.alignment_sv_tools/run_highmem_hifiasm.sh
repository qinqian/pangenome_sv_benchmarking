#!/bin/bash -ex

cd output/align/NCI2009_hifi1/

#output/align/NCI2009_hifi1/chm13_mixdown.fastq.gz

#~/software/hifiasm/hifiasm -t32 -o chm13_mixdown.asm chm13_mixdown.fastq.gz > NCI2009_hifi1.logs 2>&1 
~/software/hifiasm/hifiasm -t32 -o chm13_mixdown.asm chm13_mixdown.fastq.gz > NCI2009_hifi1.stdout 2> NCI2009_hifi1.err
