#!/bin/bash -ex

export ZENODO_TOKEN=wy6tdhNY6i1nblSpu68yY1jeLuhKe8oCeGjiA7dNCSTJNx8lfMhcX4qRBYuv


main() {
####for folder in  COLO829_hifi1        COLO829_ont1        COLO829_truth_grch38  HCC1395_hifi1_mixed  HCC1937_hifi1_mixed  HCC1954_hifi1_mixed  NCI1437_hifi1_mixed  NCI2009_hifi1_mixed COLO829_hifi1_mixed  COLO829_ont1_mixed  ensembl               HCC1395_hifi1     HCC1937_hifi1        HCC1954_hifi1        NCI1437_hifi1        NCI2009_hifi1        normal_cell; do
#for folder in asm; do   
#for folder in asm_downsample  COLO829_hifi1  COLO829_hifi1_mixed  COLO829_ont1  COLO829_truth_grch38  ensembl  ensembl_downsample   HCC1395_hifi1  HCC1937_hifi1  HCC1954_hifi1  NCI1437_hifi1  NCI2009_hifi1; do
#    echo $folder
#    tar -cvzf ${folder}.tar.gz $folder
#    ../20.zenodo_latest/zenodo-upload/zenodo_upload.sh 18989691  ${folder}.tar.gz
#done

# ========== 配置 ==========
SRC_DIR="."                          # 源目录（包含所有样本文件夹）
FLAT_DIR="zenodo_flat_upload"        # 扁平化输出目录
# ==========================

#mkdir -p "$FLAT_DIR"
#
#folders="asm asm_downsample COLO829_hifi1 COLO829_hifi1_mixed COLO829_ont1 COLO829_truth_grch38 ensembl ensembl_downsample HCC1395_hifi1 HCC1937_hifi1 HCC1954_hifi1 NCI1437_hifi1 NCI2009_hifi1"
#
#for top in $folders; do
#    [ -d "$SRC_DIR/$top" ] || { echo "Skip missing: $top"; continue; }
#    echo "=== Processing: $top ==="
#    
#    find "$SRC_DIR/$top" -type f | while read -r filepath; do
#        # 去掉顶层目录前缀，得到相对路径
#        relpath="${filepath#$SRC_DIR/$top/}"
#        
#        # 把路径中的 / 全部替换成 __
#        # 例如: nanomonsv/result.vcf  ->  nanomonsv__result.vcf
#        flatname=$(echo "$relpath" | sed 's|/|__|g')
#        
#        # 最终文件名: 顶层__相对路径
#        newname="${top}__${flatname}"
#        dest="$FLAT_DIR/$newname"
#        
#        # 冲突检测：如果重名，自动加 _dup1, _dup2...
#        if [ -e "$dest" ]; then
#            base="${newname%.*}"
#            ext="${newname##*.}"
#            counter=1
#            if [ "$base" = "$newname" ]; then
#                # 无扩展名的情况
#                while [ -e "$dest" ]; do
#                    dest="$FLAT_DIR/${newname}_dup${counter}"
#                    counter=$((counter + 1))
#                done
#            else
#                while [ -e "$dest" ]; do
#                    dest="$FLAT_DIR/${base}_dup${counter}.${ext}"
#                    counter=$((counter + 1))
#                done
#            fi
#            echo "WARNING: Name conflict, renamed to $(basename "$dest")"
#        fi
#        
#        cp -v "$filepath" "$dest"
#    done
#done
#
#echo ""
#echo "✅ Done. Flattened files are in: ./$FLAT_DIR"
#echo "Total files: $(find "$FLAT_DIR" -type f | wc -l)"
#echo ""
#echo "Preview (first 20):"
#ls -1 "$FLAT_DIR" | head -20

#for i in ${FLAT_DIR}/*; do
#    ../20.zenodo_latest/zenodo-upload/zenodo_upload.sh 21852833 ${i}
#done

#for i in zenodo_flat_upload/*l+s*3_filtered.vcf; do
#    if [[ "$filename" == *downsample* ]]; then
#        echo "Skipping downsample file: $i"
#        continue
#    fi
#    ../20.zenodo_latest/zenodo-upload/zenodo_upload.sh 21852833 ${i}
#done

#for i in zenodo_flat_upload/*l+g_3_filtered.vcf; do
#    if [[ "$filename" == *downsample* ]]; then
#        echo "Skipping downsample file: $i"
#        continue
#    fi
#    ../20.zenodo_latest/zenodo-upload/zenodo_upload.sh 21852833 ${i}
#done

echo "Hello world"
../20.zenodo_latest/zenodo-upload/zenodo_upload.sh 21852833 COLO829_truth_grch38/COLO829.truth.hs38.vcf

}

main
