library(ggplot2)
library(stringr)
library(tidyverse)
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
	print(res)
	res = as.data.frame(res)[, c(-1, -15)]
	print(res)

        # sensitivity
        sensitivity = res %>% select(X2) %>% pull()
	tp = sensitivity[2:length(sensitivity)] * sensitivity[1]
	fn = as.integer(sensitivity[1] - tp)

        ## precision
	total_pred = as.numeric(diag(as.matrix(res)))
	total_pred = total_pred[2:length(total_pred)]
        specificity = as.numeric(unlist(as.data.frame(res)[1, 2:ncol(res)]))
	tp = as.integer(total_pred * specificity)
	fp = total_pred - tp
        print(total_pred)
        print(tp)
        count = str_extract(path, "count(\\d+)", group=1)
	metrics = as_tibble(data.frame(FP=fp, FN=fn, tools=tools[-1], count=count))
	metrics = metrics %>% mutate(cell_line=gsub("_hifi1_grch38_count\\d_eval.tsv", "", basename(path)))
        df_list[[path]] = metrics
    }

    metrics = as_tibble(do.call(rbind, df_list))
    metrics = metrics %>% pivot_longer(cols=starts_with('F'), names_to="metrics", values_to = "SV_num")

    metrics = metrics %>% mutate(tool = ifelse(
      grepl("severus", tools),
      "severus", 
      ifelse(grepl("snf", tools),
             "snf",
	     ifelse(grepl("msv_ltgs", tools),
             "msv_ltgs",
	     ifelse(grepl("l\\+tg\\.pair", tools),
             "msv_ltg",
	     ifelse(grepl("msv_lts_", tools),
             "msv_lts",
	     ifelse(grepl("l\\+t\\+s_mosaic", tools),
             "msv_lts",
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
    ))))
    metrics = metrics %>% mutate(tool=paste0(tool, ifelse(grepl("asm\\.vcf", tools), "_asm", "")))
    metrics = metrics %>% filter(tool != "msv_lt")

    write_tsv(metrics, out_path[['stat']])
    print(head(metrics))

    pdf(out_path[['bar_pdf']], width=12, height=3.2)
    p = ggplot(data=metrics) + geom_bar(aes(x=tool, y=SV_num, fill=factor(metrics)), stat = "identity", position = 'stack') + xlab("different tools") + facet_grid(cell_line~count, scales = "free") + get_theme(angle=30, size=10) + scale_fill_manual(values = custom_colors) + ggtitle("COLO829 tumor-normal paired somatic SV") + ylab("FN/FP SV calls number")
    print(p)
    dev.off()
}


do_bar_chart(snakemake@input, snakemake@output, snakemake@threads, snakemake@params)

