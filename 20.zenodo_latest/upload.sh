#!/bin/bash -ex

export ZENODO_TOKEN=wy6tdhNY6i1nblSpu68yY1jeLuhKe8oCeGjiA7dNCSTJNx8lfMhcX4qRBYuv


main() {
####for folder in  COLO829_hifi1        COLO829_ont1        COLO829_truth_grch38  HCC1395_hifi1_mixed  HCC1937_hifi1_mixed  HCC1954_hifi1_mixed  NCI1437_hifi1_mixed  NCI2009_hifi1_mixed COLO829_hifi1_mixed  COLO829_ont1_mixed  ensembl               HCC1395_hifi1     HCC1937_hifi1        HCC1954_hifi1        NCI1437_hifi1        NCI2009_hifi1        normal_cell; do

for folder in  COLO829_hifi1_mixed COLO829_hifi1        COLO829_ont1        COLO829_truth_grch38  ensembl               HCC1395_hifi1     HCC1937_hifi1        HCC1954_hifi1        NCI1437_hifi1        NCI2009_hifi1        normal_cell; do
    echo $folder
    tar -cvzf ${folder}.tar.gz $folder
    zenodo-upload/zenodo_upload.sh 17317401  ${folder}.tar.gz
done
}

main
