options(warn=1)
library(ggplot2)
library(stringr)
library(tidyverse)
library(ggpattern)
library(ggthemes)
library(patchwork)

custom_colors <- c(
  "FP" = rgb(46, 108, 128, maxColorValue = 255),
  "FN" = rgb(233, 127, 90, maxColorValue = 255)
)


get_theme <- function(size=12, angle=0) {
    defined_theme = theme_clean(base_size=size) + theme(legend.title=element_text(size=size), strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05)) #, legend.position="bottom", legend.box = "horizontal") 
    defined_theme
}


do_bar_chart <- function(data_path, out_path, threads, myparam) {
    # R code
    df_list = list()
    for (path in unlist(data_path)) {
        res = read_tsv(path, col_names = FALSE)
        tools = res[, ncol(res)] %>% pull()

        # interesting
        # in total 15 cols, will change as we append more tools
	res = as.data.frame(res)[, c(-1, -15)]

        # sensitivity
        sensitivity = res %>% select(X2) %>% pull()
	tp = sensitivity[2:length(sensitivity)] * sensitivity[1]
	fn = as.integer(sensitivity[1] - tp)

        ## precision
        print(diag(as.matrix(res)))
	total_pred = as.numeric(diag(as.matrix(res)))
	total_pred = total_pred[2:length(total_pred)]
        print(unlist(as.data.frame(res)[1, 2:ncol(res)]))
        specificity = as.numeric(unlist(as.data.frame(res)[1, 2:ncol(res)]))
	tp = as.integer(total_pred * specificity)
	fp = total_pred - tp
        count = str_extract(path, "count(\\d+)", group=1)
	metrics = as_tibble(data.frame(FP=fp, FN=fn, tools=tools[-1], count=count))
	metrics = metrics %>% mutate(
				     cell_line=gsub("_100k", "", gsub("_hifi1_grch38_count\\d_eval.tsv", "", basename(path))))
        df_list[[path]] = metrics
    }

    metrics = as_tibble(do.call(rbind, df_list))
    metrics$count = as.numeric(metrics$count)
    metrics = metrics %>% pivot_longer(cols=starts_with('F'), names_to="metrics", values_to = "SV_num")
    print(head(metrics))

    metrics = metrics %>% mutate(tool = ifelse(
      grepl("severus_lowaf25", tools),
      "severus", 
      ifelse(grepl("sniffles|snf", tools),
             "snf",
	     ifelse(grepl("grch38l_l\\+t_mosaic", tools),
             "msv_lt",
	     ifelse(grepl("grch38l_l\\+t\\+g\\+s_mosaic", tools),
             "msv_ltgs",
	     ifelse(grepl("grch38l_l\\+t\\+g_mosaic", tools),
             "msv_ltg",
	     ifelse(grepl("grch38l_l\\+t\\+s_mosaic", tools),
             "msv_lts",
             ifelse(
             grepl("_msv_", tools),
             "msv_ltgs_asmonlyname",
	     NA
             )
	     ))
	     )
	     )
      )
    ))

    metrics = metrics %>% mutate(group=ifelse(grepl("_asm\\.vcf", tools), "caller+asm", "caller"))
    metrics = metrics %>% mutate(group2=ifelse(grepl("readname_asm", tools), "_onlyname", ""))
    metrics = metrics %>% mutate(group=paste0(group, group2))

    metrics = metrics %>% filter(tool != "msv_lt")
    metrics = metrics %>% filter(tool != "msv_lts")
    write_tsv(metrics, out_path[['stat']])

    metrics = metrics %>% mutate(
        tools=case_when(
            tool=="msv_ltg" ~ "minisv",
            tool=="msv_ltgs" ~ "minisv",
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

    pdf(out_path[['bar_pdf2']], width=18, height=4)
    p1 = ggplot(data=metrics_cutoff2) +
      geom_bar_pattern(
          aes(
                x=group, y=SV_num,
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
                x=group, y=SV_num,
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
    print(p1 + p2)
    dev.off()

    pdf(out_path[['bar_pdf']], width=9, height=3.5)
    p = ggplot(data=metrics) + geom_bar(aes(x=tool, y=SV_num, fill=factor(metrics)), stat = "identity", position = 'stack') + xlab("different tools") + ylab("FN/FP SV calls number") + facet_grid(cell_line~count, scales = "free") + get_theme(angle=30, size=10) + scale_fill_manual(values = custom_colors) + ggtitle("Mosaic SV calling for tumor-normal mixed sample") 
    print(p)
    dev.off()
}

do_bar_chart(snakemake@input, snakemake@output, snakemake@threads, snakemake@params)

