library(ggplot2)
library(ggbreak)
library(tidyverse)
library(ggthemes)
library(patchwork)

custom_colors <- c(
  ">100k" = rgb(245, 106, 100, maxColorValue = 255),
  ">1M" = rgb(24, 176, 61, maxColorValue = 255),
  "translocation" = rgb(85, 148, 249, maxColorValue = 255),
  "all_except_inv" = rgb(21, 183, 187, maxColorValue = 255),
  ">20k" = rgb(221, 180, 104, maxColorValue = 255)
)


defined_theme = theme_clean() + theme(legend.title=element_text(size=5), strip.text=element_text(size=5), legend.text=element_text(size=5), axis.title.x=element_text(size=5), axis.title.y=element_text(size=5), axis.text.y=element_text(size=5), axis.text.x=element_text(size=5, angle=90, hjust=1, vjust=0.5), panel.background = element_blank()) 


preprocess_data <- function(x, count_thresh=2) {
    res = read_tsv(x)
    res = res %>% filter(count == count_thresh)
    # columns: translocation   >1M     >100k   >20k    all_except_inv  inv     file
    res.small =res[, c(4,5,7)]
    # 4: >20k, 5, all except inv
    # 5 -> 100bp - 20kb
    res.small[, 2] = res.small[, 2] - res[, 1]
    # 4 -> 20k - 100k
    res.small[, 1] = res.small[, 1] - res[, 3]

    res = res[, c(-5,-6)]
    # subtract translocation
    res[,4] = res[,4] - res[,3]
    res[,3] = res[,3] - res[,2]
    res[,2] = res[,2] - res[,1] 
    # remove >20k
    res = res[,-4]
    res = res %>% pivot_longer(cols=c(`translocation`, `>1M`, `>100k`), names_to='Size')
    res.small = res.small %>% pivot_longer(cols=c('all_except_inv', '>20k'), names_to='Size')
    return(list(res, res.small))
}


clean_meta <- function(x, param) {
    # identify the facet genome information
    x$genome = ifelse(grepl("chm13", x$file), "chm13", "hg38")

    # get the file name
    x$file = ifelse(grepl("severus|sniffle", x$file), dirname(x$file), gsub(".c2s0.msv", "", basename(x$file)))

    # trim file name as label to align chm13 and hg38 minisv names
    x$file = ifelse(x$genome == 'chm13' & (!grepl("l\\+x", x$file)), gsub("l\\+", "l\\+t", x$file), x$file)

    x.normal = x %>% filter(!grepl(param[['select']], file))
    x = x %>% filter(grepl(param[['select']], file))

    # match cell line in the severus
    x$cell_line = ifelse(grepl("severus", x$file),
                         gsub("hifi1", "BL", gsub("_", "", str_extract(x$file, '\\w+_BL'))),
                         ifelse(grepl("sniffle", x$file), gsub("hifi1", "BL", gsub("_", "", str_extract(x$file, '\\w+_hifi1'))), 
                                apply(matrix(x$file), 1, function(x) unlist(strsplit(x, '\\.'))[1])))

    x.normal$cell_line = ifelse(grepl("severus|sniffle", x.normal$file), 
                                str_extract(x.normal$file, "HG002|HG00099|HG01192|NA18983|HG03225"), 
                                apply(matrix(x.normal$file), 1, function(x) unlist(strsplit(x, '\\.'))[1]))

    x$file = gsub("hg38|chm13", "", x$file)
    x.normal$file = gsub("hg38|chm13", "", x.normal$file)

    x$file[grepl('severus', x$file)] = 'Severus'
    x$file[grepl('sniffles\\/', x$file)] = 'snf'
    x$file[grepl('sniffles_mosaic\\/', x$file)] = 'snf_mosaic'

    x.normal$file[grepl('severus', x.normal$file)] = 'Severus'
    x.normal$file[grepl('sniffles\\/', x.normal$file)] = 'snf'
    x.normal$file[grepl('sniffles_mosaic\\/', x.normal$file)] = 'snf_mosaic'

    write_tsv(x, file="paired_normal_hifi_sv_length.tsv")
    write_tsv(x.normal, file="hprc_normal_hifi_sv_length.tsv")
    list(x, x.normal)
}


plot_bar <- function(x, is_germline=F, add_break=F) {
    if (is_germline) {
        p=ggplot(x, aes(x=reorder(file, value), y=value, fill=Size)) + geom_bar(position='stack', stat='identity') + facet_grid(genome~cell_line, scales='free')+ylab(expression("#germline SV ">="5 Supporting reads")) + xlab("") + defined_theme
    } else {
        p=ggplot(x, aes(x=reorder(file, value), y=value, fill=Size)) + geom_bar(position='stack', stat='identity') + facet_grid(genome~cell_line, scales='free')+ylab(expression("#mosaic SV ">="2 Supporting reads")) + xlab("") + defined_theme
    }
    #if (add_break) {
    #    p = p + scale_y_break(c(100, 1000), scales="free")
    #}
    p
}


do_bar_chart <- function(data_path, out_path, threads, myparam) {
    # R code
    print(myparam)
    print(data_path)
    parsed_res = preprocess_data(data_path, count_thresh=2)
    res = parsed_res[[1]]
    res.small = parsed_res[[2]]

    parsed_res_germline = preprocess_data(data_path, count_thresh=5)
    res.germline = parsed_res_germline[[1]]
    res.small.germline = parsed_res_germline[[2]]

    cleaned_res = clean_meta(res, myparam)
    # paired normal
    res = cleaned_res[[1]]
    # normal hprc
    res.normal = cleaned_res[[2]]

    parsed_res_small = clean_meta(res.small, myparam)
    res.small = parsed_res_small[[1]]
    res.small.normal = parsed_res_small[[2]]

    clean_res_germline = clean_meta(res.germline, myparam)
    # germline SV parser
    res.germline = clean_res_germline[[1]]
    res.normal.germline = clean_res_germline[[2]]

    clean_res_small_germline = clean_meta(res.small.germline, myparam)
    res.small.germline = clean_res_small_germline[[1]]
    res.small.germline.normal = clean_res_small_germline[[2]]

    print(table(res$cell_line))
    print(table(res.normal$cell_line))

    res.small = res.small %>% filter(!grepl('l\\+x$|l\\+t$', res.small$file))

    res.small = res.small %>% filter(!grepl('Severus', res.small$file))
    res.small = res.small %>% filter(!grepl('^snf$', res.small$file))

    # remove caller with too high SV numbers
    res.small.normal = res.small.normal %>% filter(!grepl('l\\+x$|l\\+t$', res.small.normal$file))
    res.small.germline.normal = res.small.germline.normal %>% filter(!grepl('^snf_mosaic$', res.small.germline.normal$file))
    res.normal.germline = res.normal.germline %>% filter(!grepl('^snf_mosaic$', res.normal.germline$file))
    res.small.normal = res.small.normal %>% filter(!grepl('Severus', res.small.normal$file))
    res.small.normal = res.small.normal %>% filter(!grepl('^snf$', res.small.normal$file))

    pdf(out_path[['paired']], width=8.6, height=3)
    p_big = plot_bar(res)
    p_small = plot_bar(res.small)
    print(p_big + p_small+plot_layout(guides = "collect") + scale_fill_manual(values = custom_colors))
    dev.off()

    pdf(out_path[['normal']], width=8.6, height=3)
    p_big = plot_bar(res.normal)
    p_small = plot_bar(res.small.normal)
    print(p_big + p_small+plot_layout(guides = "collect") + scale_fill_manual(values = custom_colors))
    dev.off()

    pdf(out_path[['normal_germline']], width=8.6, height=3)
    p_big = plot_bar(res.normal.germline, is_germline=T)
    p_small = plot_bar(res.small.germline.normal, is_germline=T)
    print(p_big + p_small + plot_layout(guides = "collect") + scale_fill_manual(values = custom_colors))
    dev.off()

    res.hg38 = res %>% filter(genome == 'hg38')
    res.chm13 = res %>% filter(genome == 'chm13')
    res.small.hg38 = res.small %>% filter(genome == 'hg38')
    res.small.chm13 = res.small %>% filter(genome == 'chm13')

    pdf(out_path[['hg38paired']], width=8.6, height=3)
    p_big = plot_bar(res.hg38)
    p_small = plot_bar(res.small.hg38)
    print(p_big + p_small+plot_layout(guides = "collect") + scale_fill_manual(values = custom_colors))
    dev.off()

    pdf(out_path[['chm13paired']], width=8.6, height=3)
    p_big = plot_bar(res.chm13)
    p_small = plot_bar(res.small.chm13)
    print(p_big + p_small+plot_layout(guides = "collect") + scale_fill_manual(values = custom_colors))
    dev.off()

    res.normal.hg38 = res.normal %>% filter(genome == 'hg38')
    res.normal.chm13 = res.normal %>% filter(genome == 'chm13')
    res.small.normal.hg38 = res.small.normal %>% filter(genome == 'hg38')
    res.small.normal.chm13 = res.small.normal %>% filter(genome == 'chm13')

    res.normal.germline.hg38 = res.normal.germline %>% filter(genome == 'hg38')
    res.normal.germline.chm13 = res.normal.germline %>% filter(genome == 'chm13')
    res.small.normal.germline.hg38 = res.small.germline.normal %>% filter(genome == 'hg38')
    res.small.normal.germline.chm13 = res.small.germline.normal %>% filter(genome == 'chm13')

    pdf(out_path[['hg38normal']], width=8.6, height=3.2)
    p_big = plot_bar(res.normal.hg38)
    p_small = plot_bar(res.small.normal.hg38)
    print(p_big + p_small+plot_layout(guides = "collect") + scale_fill_manual(values = custom_colors))
    dev.off()

    pdf(out_path[['hg38normal_germline']], width=8.6, height=3.2)
    p_big   = plot_bar(res.normal.germline.hg38, is_germline=T)
    p_small = plot_bar(res.small.normal.germline.hg38, is_germline=T)
    print(p_big + p_small+plot_layout(guides = "collect") + scale_fill_manual(values = custom_colors))
    dev.off()

    pdf(out_path[['chm13normal']], width=8.6, height=3.2)
    p_big = plot_bar(res.normal.chm13)
    p_small = plot_bar(res.small.normal.chm13)
    print(p_big + p_small+plot_layout(guides = "collect") + scale_fill_manual(values = custom_colors))
    dev.off()
}

do_bar_chart(snakemake@input[[1]], snakemake@output, snakemake@threads, snakemake@params)
