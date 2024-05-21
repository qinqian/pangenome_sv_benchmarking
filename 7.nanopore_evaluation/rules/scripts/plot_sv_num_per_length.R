library(ggplot2)
library(tidyverse)

do_bar_chart <- function(data_path, out_path, threads, myparam) {
    # R code
    print(data_path)
    res = read_tsv(data_path)
    res =res[,-5]
    res = res %>% pivot_longer(cols=c(`>100`, `>1M`, `>100k`, `>20k`,  inv))
    res$file = gsub("_mergeflt_c5s0.gsv.gz", "", basename(res$file))
    res$genome = ifelse(grepl("chm13", res$file), "chm13", "hg38")
    ggplot(res, aes(x=file, y=value, fill=name)) + geom_bar(position='stack', stat='identity') + facet_wrap(~genome, ncol=2, scales='free')+theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + ylab(expression("#SV ">="5 Supporting reads")) + xlab("")
    ggsave(out_path)
}

do_bar_chart(snakemake@input[[1]], snakemake@output[[1]], snakemake@threads, snakemake@config)
