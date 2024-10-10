library(ggplot2)
library(ggthemes)
library(stringr)
library(tidyverse)
library(patchwork)

custom_colors <- c(
  "msv:tg" = rgb(249, 134, 130, maxColorValue = 255),
  "minisvl+tg" = rgb(249, 134, 130, maxColorValue = 255),
  "minisvl+tgs" = rgb(187, 142, 33, maxColorValue = 255),
  "msv:tgsnonormal" = rgb(22, 183, 139, maxColorValue = 255),
  "msv:tgs" = rgb(187, 142, 33, maxColorValue = 255),
  "minisvl+gs" = rgb(187, 142, 33, maxColorValue = 255),
  "minisvl+g" = rgb(249, 134, 130, maxColorValue = 255),
  "nanomonsv" = rgb(147, 203, 118, maxColorValue = 255),
  "savana" = rgb(22, 183, 139, maxColorValue = 255),
  "severus_mosaic" = rgb(16, 174, 228, maxColorValue = 255),
  "severus_somatic" = "orange",
  "severus_lowaf" = rgb(16, 174, 228, maxColorValue = 255),
  "sniffles" = rgb(154, 130, 251, maxColorValue = 255),
  "sniffles_somatic" = rgb(248, 89, 206, maxColorValue = 255),
  "snf_mosaic" = rgb(154, 130, 251, maxColorValue = 255),
  "svision" = rgb(248, 89, 206, maxColorValue = 255),
  "severus_mosaic_filterasm" = "black",
  "severus_mosaic_rawasm" = "gray"
)

get_theme <- function(size=12, angle=0) {
    defined_theme = theme_bw(base_size=size) + theme(legend.title=element_text(size=size), strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05)) #, legend.position="bottom", legend.box = "horizontal") 
    defined_theme
}

load_metrics <- function(data_path) {
    df_list = list()
}

do_bar_chart <- function(input, out_path, threads, myparam) {
    data_path = unlist(input[['stat']])
    data_list = list()
    for (file in data_path) {
        data_list[[file]] = read_tsv(file, col_names=T)
        data_list[[file]] = data_list[[file]] %>% mutate(cell_line = gsub("_hifi1_somatic_generation2_filterasm.stat", "", basename(file)))
    }

    metrics = bind_rows(data_list)
    print(head(metrics))
    write_tsv(metrics, out_path[['stat']])

    metrics$asm_support[metrics$asm_support >= quantile(metrics$asm_support, 0.99)] = quantile(metrics$asm_support, 0.99)
    metrics$read_name_number[metrics$read_name_number >= quantile(metrics$read_name_number, 0.99)] = quantile(metrics$read_name_number, 0.99)

    pdf(out_path[['pdf']], width=9.5, height=6)
    p = ggplot() + geom_point(data=metrics, aes(x=read_name_number, y=asm_support, colour=svlen, shape=svlen), alpha=0.6) + xlab("Before filtering") + ylab("After filtering") + 
        facet_wrap(~cell_line, ncol=3, nrow=2, scales = "free") + get_theme()
    print(p)
    dev.off()
}

do_bar_chart(snakemake@input, snakemake@output, snakemake@threads, snakemake@params)


