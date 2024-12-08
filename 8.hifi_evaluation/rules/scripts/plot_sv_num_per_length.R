library(ggplot2)
library(tidyverse)

do_bar_chart <- function(data_path, out_path, threads, myparam) {
    # R code
    print(data_path)
    res = read_tsv(data_path)
    print(res)
    # remove inversions and all SVs
    res = res[,c(-5,-6)]
    res[,4] = res[,4] - res[,3]
    res[,3] = res[,3] - res[,2]
    res[,2] = res[,2] - res[,1] 
    # remove >20k
    res = res[,-4]
    #res = res %>% pivot_longer(cols=c(`translocation`, `>1M`, `>100k`, `>20k`))
    res = res %>% pivot_longer(cols=c(`translocation`, `>1M`, `>100k`))
    #res$file = gsub("_mergeflt_c5s0.gsv.gz", "", basename(res$file))
    res$tool = str_extract(res$file, "(GRCh38|CHM13v2|gafcall|minisv_pair|minisv_mosaic|nanomonsv|savana|severus|sniffles_mosaic|sniffles|svision|cutesv)", group=1)
    res$cell_line = str_extract(res$file, "(HG002rep1|HG002rep2)", group=1)
    res$cell_line[is.na(res$cell_line)] = 'truth'
    res$tool[res$tool == 'gafcall'] = 'minisv'
    res$params = str_extract(res$file, "(l_g|l_s|l_ts|l_tg|l_tgs|l)_(join_merge|merge)", group=1)
    res$params[is.na(res$params)] = ''
    res = res %>% filter(params!='l')
    #res = as_tibble(res) %>% filter(tool == 'gafcall')
    res$genome = ifelse(grepl("chm13", res$file, ignore.case=T), "chm13", "hg38")
    res$file = gsub("chm13|hg38", "", basename(res$file))
    res$tool = gsub("GRCh38|CHM13v2", "truthset", res$tool)
    res$comb = ifelse(res$params=="", res$cell_line, paste0(res$cell_line, ":", res$params))

    ggplot(res, aes(x=reorder(comb, value), y=value, fill=name)) + geom_bar(position='stack', stat='identity') + facet_wrap(~genome+tool, ncol=6, scales='free_x')+theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + ylab(expression("#SV ">="5 Supporting reads")) + xlab("") #+ ylim(0, 50)
    ggsave(out_path[['pdf']], width=9.5, height=6)
    write_tsv(res, out_path[['tsv']])
}

do_bar_chart(snakemake@input[[1]], snakemake@output, snakemake@threads, snakemake@config)

