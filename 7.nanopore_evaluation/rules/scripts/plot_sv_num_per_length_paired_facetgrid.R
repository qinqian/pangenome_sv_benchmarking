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
    res = res %>% filter(grepl(myparam[['select']], file))
    # remove >20k
    res = res[,-4]

    res = res %>% pivot_longer(cols=c(`translocation`, `>1M`, `>100k`))
    res$file = gsub("_merge_c5s0.msv.gz", "", basename(res$file))

    res$cell_line = ifelse(grepl("ont", res$file), 
                           apply(matrix(res$file), 1, function(x) paste(unlist(strsplit(x, '_'))[1:3], collapse="")),
                           apply(matrix(res$file), 1, function(x) paste(unlist(strsplit(x, '_'))[1:2], collapse="")))
    print(res$cell_line)

    res$genome = ifelse(grepl("chm13", res$file), "chm13", "hg38")
    res$file = gsub("chm13|grch38", "", basename(res$file))
    #res$file = gsub("ont\\d", "", basename(res$file))

    ggplot(res, aes(x=reorder(file, value), y=value, fill=name)) + geom_bar(position='stack', stat='identity') + facet_grid(genome~cell_line, scales='free')+theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + ylab(expression("#SV ">="5 Supporting reads")) + xlab("")
    ggsave(out_path, width=11.5, height=5)
}

do_bar_chart(snakemake@input[[1]], snakemake@output[[1]], snakemake@threads, snakemake@params)
