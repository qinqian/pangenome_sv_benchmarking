options(warn=1)
library(ggplot2)
library(stringr)
library(tidyverse)
library(ggpattern)
library(ggthemes)
library(patchwork)


custom_colors <- c(
  "FP" = rgb(46, 108, 128, maxColorValue = 255),
  "TP" = rgb(233, 127, 90, maxColorValue = 255)
)


get_theme <- function(size=12, angle=0) {
    #defined_theme = theme_clean(base_size=size) + theme(legend.title=element_text(size=size), strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05), legend.box.spacing = unit(0, "mm"))
    defined_theme = theme_classic(base_size=size) + theme(legend.title=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05), legend.box.spacing = unit(0, "mm"), strip.text = element_blank())
    defined_theme
}


parse_evaluation <- function(data_path) {
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
	#fn = round(sensitivity[1] - tp)
        #fn_rate = 1 - sensitivity[2:length(sensitivity)]

        ## precision
	total_pred = as.numeric(diag(as.matrix(res)))
	total_pred = total_pred[2:length(total_pred)]
        specificity = as.numeric(unlist(as.data.frame(res)[1, 2:ncol(res)]))
	tp = round(total_pred * specificity)
	fp = total_pred - tp
        fp_rate = 1 - specificity

        count = str_extract(path, "count(\\d+)", group=1)
	metrics_fp = as_tibble(data.frame(metrics='FP', SV_num=fp, tools=tools[-1], count=count))
	metrics_fn = as_tibble(data.frame(metrics='TP', SV_num=tp, tools=tools[-1], count=count))
        metrics = bind_rows(metrics_fp, metrics_fn)

	metrics = metrics %>% mutate(cell_line=gsub("_100k", "", gsub("_hifi1_grch38_count\\d_eval.tsv", "", basename(path))))
        df_list[[path]] = metrics
    }

    metrics = as_tibble(do.call(rbind, df_list))
    metrics$count = as.numeric(metrics$count)

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

    metrics = metrics %>% mutate(group=ifelse(grepl("_asm\\.vcf", tools), "asm_keep", "caller"))
    metrics = metrics %>% mutate(group2=ifelse(grepl("readname_asm", tools), "_onlyname", ""))

    # remove only read name evaluation
    metrics = metrics %>% filter(group2 != '_onlyname')

    metrics = metrics %>% mutate(group=paste0(group, group2))

    metrics = metrics %>% filter(tool != "msv_lt")
    metrics = metrics %>% filter(tool != "msv_lts")

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
            tool=="msv_ltgs" ~ "asm_keep",
            tool=="msv_ltgs_asmonlyname" ~ "caller+asm_onlyname",
            .default = group
        )
    )

    metrics = metrics %>% select(metrics,SV_num,tools,count,cell_line,group)
    metrics = metrics %>% pivot_wider(names_from=c("group"), values_from='SV_num')
    metrics = metrics %>% mutate(asm_filter=abs(caller-asm_keep))
    metrics = metrics %>% select(count,cell_line,metrics,tools,asm_keep,asm_filter) %>% pivot_longer(cols=c("asm_filter", "asm_keep"), names_to="group", values_to="SV_num")

    metrics_cutoff2 = metrics %>% filter(count == 2)
    metrics_cutoff2
}


do_bar_chart <- function(data_path, out_path, threads, myparam) {
    metrics = parse_evaluation(data_path['stat'])
    write_tsv(metrics, out_path[['stat']])
    
    pdf(out_path[['bar_pdf2']], width=7.5, height=3.5)
    p1 = ggplot(data=metrics %>% filter(metrics=='FP'), 
          aes(x=tools, y=SV_num,
              fill=factor(metrics),
              pattern=factor(group))) + 
      geom_bar_pattern(
          stat = "identity", position="stack",
          colour          = 'black',
          pattern_fill = "black",
          pattern_angle = 45,
          pattern_density = 0.03,
          pattern_key_scale_factor = 0.6,
          pattern_spacing = 0.05) + 
      scale_fill_manual(values = custom_colors) + 
      scale_pattern_manual(values = c(asm_keep = "none", `asm_filter` = "stripe")) +
      geom_text(aes(y=SV_num, label=SV_num), size=3,
                position = 'stack') +
      get_theme(angle=45, size=12) +
      ggtitle("COLO829 Mosaic SV evaluation (all)") + ylab("FP number at\nread count cutoff 2")+ xlab("different tools")

    p2 = ggplot(data=metrics %>% filter(metrics=='TP'), 
          aes(x=tools, y=SV_num,
              fill=factor(metrics),
              pattern=factor(group))) + 
      geom_bar_pattern(
          stat = "identity", position="stack",
          colour          = 'black',
          pattern_fill = "black",
          pattern_angle = 45,
          pattern_density = 0.03,
          pattern_key_scale_factor = 0.6,
          pattern_spacing = 0.05) + 
      scale_fill_manual(values = custom_colors) + 
      scale_pattern_manual(values = c(asm_keep = "none", `asm_filter` = "stripe")) +
      geom_text(aes(y=SV_num, label=SV_num), size=3,
                position = 'stack') +
      get_theme(angle=45, size=12) +
      geom_hline(yintercept = 58, color='red') +
      ggtitle("COLO829 Mosaic SV evaluation (all)") + ylab("TP number at\nread count cutoff 2")+ xlab("different tools")

    print((p1+p2))
    dev.off()
}

do_bar_chart(snakemake@input, snakemake@output, snakemake@threads, snakemake@params)


