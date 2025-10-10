library(ggplot2)
library(ggpattern)
library(stringr)
library(tidyverse)
library(ggthemes)
library(patchwork)


custom_colors <- c(
  "FP" = rgb(46, 108, 128, maxColorValue = 255),
  "TP" = rgb(233, 127, 90, maxColorValue = 255)
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
	print(res)
	res = as.data.frame(res)[, c(-1, -15)]
	print(res)

        # sensitivity
        sensitivity = res %>% select(X2) %>% pull()
	tp = sensitivity[2:length(sensitivity)] * sensitivity[1]
	fn = round(sensitivity[1] - tp)

        ## precision
	total_pred = as.numeric(diag(as.matrix(res)))
	total_pred = total_pred[2:length(total_pred)]
        specificity = as.numeric(unlist(as.data.frame(res)[1, 2:ncol(res)]))
	tp = round(total_pred * specificity)
	fp = total_pred - tp
        count = str_extract(path, "count(\\d+)", group=1)
	metrics = as_tibble(data.frame(FP=fp, TP=tp, tools=tools[-1], count=count))
	metrics = metrics %>% mutate(cell_line=gsub("_hifi1_grch38_count\\d_eval.tsv", "", basename(path)))
        df_list[[path]] = metrics
    }

    metrics = as_tibble(do.call(rbind, df_list))
    metrics = metrics %>% pivot_longer(cols=c("FP", "TP"), names_to="metrics", values_to = "SV_num")
    metrics$count = as.numeric(metrics$count)

    metrics = metrics %>% mutate(tool = ifelse(
      grepl("severus", tools),
      "severus", 
      ifelse(grepl("snf", tools),
             "sniffles2",
	     ifelse(grepl("msv_ltgs", tools),
             "minisv_ltg",
	     ifelse(grepl("l\\+tg\\.pair", tools),
             "minisv_ltg",
	     ifelse(grepl("msv_lts_", tools),
             "minisv_lts",
	     ifelse(grepl("savana", tools),
             "savana",
	     ifelse(grepl("nanomonsv", tools),
             "nanomonsv",
	     ifelse(grepl("hg38l\\+t\\.pair", tools),
             "msv_lt",
             NA
	     ))
	     )
	     )
      ))
    )))
    metrics = metrics %>% mutate(group=paste0(ifelse(grepl("asm\\.vcf", tools), "asm_keep", "caller")))
    metrics = metrics %>% filter(tool != "msv_lt")
    metrics = metrics %>% filter(tool != "minisv_lts")
    metrics = metrics %>% mutate(
       tool=case_when(
          tool=="minisv_ltg" ~ "minisv",
          .default=tool
       )
    )

    print(head(metrics))
    metrics = metrics[,-1] %>% pivot_wider(names_from=c("group"), values_from='SV_num')

    metrics = metrics %>% mutate(asm_filter=caller-asm_keep)
    write_tsv(metrics, out_path[['stat']])

    metrics = metrics %>% select(count,cell_line,metrics,tool,asm_keep,asm_filter) %>% pivot_longer(cols=c("asm_filter", "asm_keep"), names_to="group", values_to="SV_num")

    pdf(out_path[['bar_pdf']], width=12, height=3)
    p1 = ggplot(data=metrics %>% filter(metrics == 'FP')) +
      geom_bar_pattern(
          aes(
	        x=count, y=SV_num,
                fill=factor(metrics), 
		pattern=factor(group),
	  ), 
          stat = "identity", position = 'stack',
          colour = 'black',
	  pattern_fill = "black",
	  pattern_angle = 45,
	  pattern_density = 0.03,
	  pattern_key_scale_factor = 0.6,
          pattern_spacing = 0.05) + 
      scale_fill_manual(values = custom_colors) + 
      scale_pattern_manual(values = c(asm_keep = "none", `asm_filter` = "stripe")) +
      xlab("different tools") +
      get_theme(angle=0, size=9) +
      facet_wrap(~tool, ncol=5) + 
      ggtitle("COLO829 tumor-normal paired somatic SV") + ylab("FP SV calls number")

    p2 = ggplot(data=metrics %>% filter(metrics == 'TP')) +
      geom_bar_pattern(
          aes(
	        x=count, y=SV_num,
                fill=factor(metrics), 
		pattern=factor(group),
	  ), 
          stat = "identity", position = 'stack',
          colour = 'black',
	  pattern_fill = "black",
	  pattern_angle = 45,
	  pattern_density = 0.03,
	  pattern_key_scale_factor = 0.6,
          pattern_spacing = 0.05) + 
      scale_fill_manual(values = custom_colors) + 
      scale_pattern_manual(values = c(asm_keep = "none", `asm_filter` = "stripe")) +
      xlab("different tools") +
      facet_wrap(~tool, ncol=5) + 
      get_theme(angle=0, size=9) + geom_hline(yintercept=58, color='red')+
      ggtitle("COLO829 tumor-normal paired somatic SV") + ylab("TP SV calls number")
    print(p1+p2)
    dev.off()
}


do_bar_chart(snakemake@input, snakemake@output, snakemake@threads, snakemake@params)


