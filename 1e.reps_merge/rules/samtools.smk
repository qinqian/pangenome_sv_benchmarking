import os


rule samtools_merge:
     input:
         cram = expand("../1a.alignment_sv_tools/output/align/{{cell_line}}_{{platform}}{rep}/{{pair}}/{{assembly}}.cram", rep=[1,2]),
         crai = expand("../1a.alignment_sv_tools/output/align/{{cell_line}}_{{platform}}{rep}/{{pair}}/{{assembly}}.cram.crai", rep=[1,2])
     output:
         cram = os.path.abspath("../1a.alignment_sv_tools/output/align/{cell_line}_{platform}/{pair}/{assembly}.cram"),
         crai = os.path.abspath("../1a.alignment_sv_tools/output/align/{cell_line}_{platform}/{pair}/{assembly}.cram.crai")
     threads: 4
     shell:
         """
../1a.alignment_sv_tools/samtools/samtools merge --threads {threads} -l 1 --write-index --output-fmt CRAM --reference ../1a.alignment_sv_tools/{wildcards.assembly}.fa -o {output.cram} {input.cram} 
        """
