library(ggplot2)
library(tidyverse)
library(ggthemes)
library(patchwork)


do_bar_chart <- function(data_path, out_path, threads, myparam) {
    # R code
    print(myparam)
    print(data_path)
    res = read_tsv(data_path)
    res =res[,c(-5,-6)]
    # subtract translocation
    res[,4] = res[,4] - res[,3]
    res[,3] = res[,3] - res[,2]
    res[,2] = res[,2] - res[,1] 
    # remove >20k
    res = res[,-4]
    res = res %>% pivot_longer(cols=c(`translocation`, `>1M`, `>100k`), names_to='Size')
    res$file = gsub(".c3s0.msv", "", basename(res$file))
    res$genome = ifelse(grepl("chm13", res$file), "chm13", "hg38")
    res$file = ifelse(grepl("chm13", res$file) & (!grepl("l\\+x", res$file)), gsub("l\\+", "l\\+t", res$file), res$file)

    res.normal = res %>% filter(!grepl(myparam[['select']], file))
    res = res %>% filter(grepl(myparam[['select']], file))

    res$cell_line = apply(matrix(res$file), 1, function(x) unlist(strsplit(x, '\\.'))[1])
    res.normal$cell_line = apply(matrix(res.normal$file), 1, function(x) unlist(strsplit(x, '\\.'))[1])
    print(table(res$cell_line))
    res$file = gsub("hg38|chm13", "", res$file)
    res.normal$file = gsub("hg38|chm13", "", res.normal$file)

    ggplot(res, aes(x=reorder(file, value), y=value, fill=Size)) + geom_bar(position='stack', stat='identity') + facet_grid(genome~cell_line, scales='free')+theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + ylab(expression("#SV ">="5 Supporting reads")) + xlab("")
    ggsave(out_path[['paired']], width=9.5, height=5)

    ggplot(res.normal, aes(x=reorder(file, value), y=value, fill=Size)) + geom_bar(position='stack', stat='identity') + facet_grid(genome~cell_line, scales='free')+theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + ylab(expression("#SV ">="5 Supporting reads")) + xlab("")
    ggsave(out_path[['normal']], width=9.5, height=5)

    defined_theme = theme_clean() + theme(legend.title=element_text(size=7), strip.text=element_text(size=7), legend.text=element_text(size=7), axis.title.x=element_text(size=7), axis.title.y=element_text(size=7), axis.text.y=element_text(size=7), axis.text.x=element_text(size=7, angle=90, hjust=1, vjust=0.5)) 

    res.hg38 = res %>% filter(genome == 'hg38')
    res.chm13 = res %>% filter(genome == 'chm13')
    ggplot(res.hg38, aes(x=reorder(file, value), y=value, fill=Size)) + geom_bar(position='stack', stat='identity') + facet_grid(genome~cell_line, scales='free')+ ylab(expression("#SV ">="5 Supporting reads")) + xlab("") + defined_theme
    ggsave(out_path[['hg38paired']], width=6.5, height=3.5)

    ggplot(res.chm13, aes(x=reorder(file, value), y=value, fill=Size)) + geom_bar(position='stack', stat='identity') + facet_grid(genome~cell_line, scales='free')+defined_theme+ ylab(expression("#SV ">="5 Supporting reads")) + xlab("")
    ggsave(out_path[['chm13paired']], width=6.5, height=3.5)

    res.normal.hg38 = res.normal %>% filter(genome == 'hg38')
    res.normal.chm13 = res.normal %>% filter(genome == 'chm13')
    ggplot(res.normal.hg38, aes(x=reorder(file, value), y=value, fill=Size)) + geom_bar(position='stack', stat='identity') + facet_grid(genome~cell_line, scales='free')+defined_theme + ylab(expression("#SV ">="5 Supporting reads")) + xlab("")
    ggsave(out_path[['hg38normal']], width=6.5, height=3.5)

    ggplot(res.normal.chm13, aes(x=reorder(file, value), y=value, fill=Size)) + geom_bar(position='stack', stat='identity') + facet_grid(genome~cell_line, scales='free') + defined_theme + ylab(expression("#SV ">="5 Supporting reads")) + xlab("")
    ggsave(out_path[['chm13normal']], width=6.5, height=3.5)
}


do_bar_chart(snakemake@input[[1]], snakemake@output, snakemake@threads, snakemake@params)
