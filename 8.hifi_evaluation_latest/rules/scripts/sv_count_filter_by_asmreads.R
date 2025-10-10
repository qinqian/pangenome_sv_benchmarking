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

do_bar_chart <- function(input, out_path, threads, myparam) {
    data_path1 = unlist(input[['severus_stat']])
    data_path2 = unlist(input[['msv_stat']])
    data_path3 = unlist(input[['snf_stat']])
    data_path4 = unlist(input[['nano_stat']])

    data_list = list()
    for (file in c(data_path1, data_path2, data_path3, data_path4)) {
        data_list[[file]] = read_tsv(file, col_names=T)
        data_list[[file]] = data_list[[file]] %>% mutate(cell_line = gsub("_hifi1_somatic_generation2_filterasm.stat", "", basename(file)))
    }

    metrics = bind_rows(data_list)
    metrics = metrics %>% mutate(tool = ifelse(
      grepl("severus", cell_line),
      "severus", 
      ifelse(grepl("snf", cell_line),
             "snf",
	     ifelse(grepl("nano", cell_line),
		    "nanomonsv",
             "msv")
      )
    ))

    write_tsv(metrics, out_path[['stat']])

    metrics = metrics %>% mutate(cell_line = gsub("msv_|snf_|severus_|nanomonsv_", "", cell_line))
    print(table(metrics$tool))

    metrics_summarize.before = list()
    for (cutoff in seq(2, 10)) {
        metrics_sub = metrics %>% filter(read_name_number>=cutoff) %>% group_by(tool, cell_line) %>% summarize(total=n()) %>% ungroup()
        metrics_sub = metrics_sub %>% mutate(cutoff=cutoff)
        metrics_summarize.before[[cutoff]] = metrics_sub
    }
    metrics_summarize.before = bind_rows(metrics_summarize.before)

    metrics_summarize = list()
    for (cutoff in seq(2, 10)) {
        metrics_sub = metrics %>% filter(asm_support>=cutoff) %>% group_by(tool, cell_line) %>% summarize(total=n()) %>% ungroup()
        metrics_sub = metrics_sub %>% mutate(cutoff=cutoff)
        metrics_summarize[[cutoff]] = metrics_sub
    }

    metrics_summarize = bind_rows(metrics_summarize)
 
    metrics_summarize.before = metrics_summarize.before %>% mutate(step='beforeasm')
    metrics_summarize = metrics_summarize %>% mutate(step='afterasm')
    metrics_summarize = bind_rows(metrics_summarize.before, metrics_summarize)
    write_tsv(metrics_summarize, out_path[['stat_sv_num']])

    pdf(out_path[['curve_pdf']], width=9.5, height=6)
    p = ggplot(data=metrics_summarize, aes(x=cutoff, y=total, colour=factor(tool), shape=factor(tool))) + geom_line(aes(linetype=factor(step))) + geom_point(alpha=0.6) + xlab("Self-assembly overlapped read cutoff") + ylab("SV numbers") + facet_wrap(~cell_line, ncol=3, nrow=2, scales = "free") + get_theme() 
    print(p)
    dev.off()

}

do_bar_chart(snakemake@input, snakemake@output, snakemake@threads, snakemake@params)

