library(ggplot2)
library(tidyverse)

do_bar_chart <- function(data_path, out_path, threads, myparam) {
    # R code
    print(data_path)
    res = read_tsv(data_path)
    print(res)
    res =res[,c(-5,-6)]
    res[,4] = res[,4] - res[,3]
    res[,3] = res[,3] - res[,2]
    res[,2] = res[,2] - res[,1] 
    res = res %>% pivot_longer(cols=c(`translocation`, `>1M`, `>100k`, `>20k`))
    res$file = gsub("_mergeflt_c5s0.gsv.gz", "", basename(res$file))
    res$genome = ifelse(grepl("chm13", res$file), "chm13", "hg38")
    res$file = gsub("chm13|hg38", "", basename(res$file))
    ggplot(res, aes(x=reorder(file, value), y=value, fill=name)) + geom_bar(position='stack', stat='identity') + facet_wrap(~genome, ncol=2, scales='free_x')+theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + ylab(expression("#SV ">="5 Supporting reads")) + xlab("")
    ggsave(out_path, width=7, height=6)
}

do_bar_chart(snakemake@input[[1]], snakemake@output[[1]], snakemake@threads, snakemake@config)
