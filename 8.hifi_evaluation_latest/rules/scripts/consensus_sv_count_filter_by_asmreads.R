library(ggplot2)
library(ggthemes)
library(ggpattern)
library(stringr)
library(tidyverse)
library(patchwork)

custom_colors <- c(
  "RP" = rgb(46, 108, 128, maxColorValue = 255),
  "RN" = rgb(233, 127, 90, maxColorValue = 255)
)


get_theme <- function(size=12, angle=0) {
    defined_theme = theme_clean(base_size=size) + theme(legend.title=element_text(size=size), strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05)) #, legend.position="bottom", legend.box = "horizontal") 
    defined_theme
}


do_bar_chart <- function(input, out_path, threads, myparam) {
    data_path1 = unlist(input[['stat']])
    data_list = list()
    for (file in c(data_path1)) {
        res = read_tsv(file, col_names=F)
        count = str_extract(file, "count(\\d+)", group=1)
        res = res %>% mutate(cell_line = gsub("_hifi1_grch38_count\\d_eval.tsv", "", basename(file)))
        res = res %>% mutate(cell_line = gsub("_hifi1_grch38_count\\d_eval_100k.tsv", "", cell_line))
        res = res %>% mutate(count = count)
        res = res %>% mutate(eval_file = file)
        data_list[[file]] = res
    }

    metrics = do.call(rbind, data_list)
    metrics$count = as.numeric(metrics$count)
    metrics = metrics %>% filter(count >= 2)
    print(head(metrics))

    metrics = bind_rows(data_list)
    print(metrics)
    metrics = metrics %>% mutate(tool = ifelse(
      grepl("severus_lowaf25", X5),
      "severus", 
      ifelse(grepl("sniffles|snf_", X5),
             "snf",
	     ifelse(grepl("l\\+t_mosaic", X5),
             "msv_lt",
	     ifelse(grepl("l\\+t\\+g\\+s_mosaic", X5),
             "msv_ltgs",
	     ifelse(grepl("l\\+t\\+g_mosaic", X5),
             "msv_ltg",
	     ifelse(grepl("l\\+t\\+s_mosaic", X5),
             "msv_lts",
             ifelse(
             grepl("_msv_", X5),
             "msv_ltgs_asmonlyname",
	     NA
             )
	     ))
	     )
	     )
      )
    ))

    metrics = metrics %>% mutate(group=paste0(ifelse(grepl("asm\\.vcf", X5), "caller+asm", "caller")))
    metrics = metrics %>% mutate(group2=ifelse(grepl("readname_asm", X5), "_onlyname", ""))
    metrics = metrics %>% mutate(group=paste0(group, group2))

    metrics = metrics %>% filter(tool != "msv_lt")
    metrics = metrics %>% filter(tool != "msv_lts")
    metrics = metrics %>% mutate(
        tools=case_when(
            tool=="msv_ltg" ~ "minisv_ltg",
            tool=="msv_ltgs" ~ "minisv_ltg",
            tool=="msv_ltgs_asmonlyname" ~ "minisv",
            tool=="snf" ~ "sniffles2", 
            .default = tool
        )
    )
    metrics = metrics %>% mutate(
        group=case_when(
            tool=="msv_ltg" ~ "caller",
            tool=="msv_ltgs" ~ "caller+asm",
            tool=="msv_ltgs_asmonlyname" ~ "caller+asm_onlyname",
            .default = group
        )
    )
    metrics_cutoff2 = metrics %>% filter(count == 2)
    metrics_cutoff3 = metrics %>% filter(count == 3)
    metrics_cutoff4 = metrics %>% filter(count == 4)
    colnames(metrics)[1] = 'metrics'
    colnames(metrics_cutoff2)[1] = 'metrics'
    colnames(metrics_cutoff3)[1] = 'metrics'
    colnames(metrics_cutoff4)[1] = 'metrics'
    write_tsv(metrics, out_path[['stat']])
    write_tsv(metrics_cutoff2, out_path[['stat2']])

    pdf(out_path[['bar_pdf2']], width=18, height=9)
    p1 = ggplot(data=metrics_cutoff2) +
      geom_bar_pattern(
          aes(
	        x=group, y=X4,
                fill=factor(metrics), 
		pattern=factor(group),
	  ), 
          stat = "identity", position = 'stack',
          colour          = 'black',
	  pattern_fill = "black",
	  pattern_angle = 45,
	  pattern_density = 0.03,
	  pattern_key_scale_factor = 0.6,
          pattern_spacing = 0.05) + 
      scale_fill_manual(values = custom_colors) + 
      scale_pattern_manual(values = c(caller = "none", `caller+asm` = "stripe", `caller+asm_onlyname` = "stripe")) +
      xlab("different tools") +
      facet_grid(cell_line~tools) +
      get_theme(angle=45, size=12) +
      ggtitle("Mosaic SV evaluation") + ylab("RN/RP ratio at\nread count cutoff 2")
    p2 = ggplot(data=metrics_cutoff3) +
      geom_bar_pattern(
          aes(
	        x=group, y=X4,
                fill=factor(metrics), 
		pattern=factor(group),
	  ), 
          stat = "identity", position = 'stack',
          colour          = 'black',
	  pattern_fill = "black",
	  pattern_angle = 45,
	  pattern_density = 0.03,
	  pattern_key_scale_factor = 0.6,
          pattern_spacing = 0.05) + 
      scale_fill_manual(values = custom_colors) + 
      scale_pattern_manual(values = c(caller = "none", `caller+asm` = "stripe", `caller+asm_onlyname` = "stripe")) +
      xlab("different tools") +
      facet_grid(cell_line~tools) +
      get_theme(angle=45, size=12) +
      ggtitle("Mosaic SV evaluation") + ylab("RN/RP ratio at\nread count cutoff 3")
    p3 = ggplot(data=metrics_cutoff4) +
      geom_bar_pattern(
          aes(
	        x=group, y=X4,
                fill=factor(metrics), 
		pattern=factor(group),
	  ), 
          stat = "identity", position = 'stack',
          colour          = 'black',
	  pattern_fill = "black",
	  pattern_angle = 45,
	  pattern_density = 0.03,
	  pattern_key_scale_factor = 0.6,
          pattern_spacing = 0.05) + 
      scale_fill_manual(values = custom_colors) + 
      scale_pattern_manual(values = c(caller = "none", `caller+asm` = "stripe", `caller+asm_onlyname` = "stripe")) +
      xlab("different tools") +
      facet_grid(cell_line~tools) +
      get_theme(angle=45, size=12) +
      ggtitle("Mosaic SV evaluation") + ylab("RN/RP ratio at\nread count cutoff 4")
    print(p1 + p2 + p3)
    dev.off()

    pdf(out_path[['bar_pdf']], width=9, height=9)
    p = ggplot(data=metrics) + geom_bar(aes(x=tools, y=X4, fill=factor(metrics)), stat = "identity", position = 'stack') + xlab("different tools") + ylab("RN/RP ratio") + facet_grid(cell_line~count, scales = "free") + get_theme(angle=30, size=10) + scale_fill_manual(values = custom_colors) + ggtitle("Mosaic SV calling for tumor-normal mixed sample")
    print(p)
    dev.off()
}


do_bar_chart(snakemake@input, snakemake@output, snakemake@threads, snakemake@params)

