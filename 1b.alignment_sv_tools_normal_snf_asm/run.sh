
main() {
 minisv.js view -c 2 -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed -l 100 ../1b.alignment_sv_tools_normal/output/sniffles_mosaic/HG002/grch38.vcf.gz   | wc -l

#minisv filterasm [OPTIONS] READIDTSV MSVASM OUTSTAT VCFFILE
#conda activate msvpy
 minisv filterasm -a -c 0 -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed -l 100 ../1b.alignment_sv_tools_normal/output/sniffles_mosaic/HG002/grch38.vcf.gz /hlilab/hli/gafcall/normal_v2/HG002.self.Q0.gsv.gz out.stat ../1b.alignment_sv_tools_normal/output/sniffles_mosaic/HG002/grch38.vcf.gz | grep -v "#" 
}
main
