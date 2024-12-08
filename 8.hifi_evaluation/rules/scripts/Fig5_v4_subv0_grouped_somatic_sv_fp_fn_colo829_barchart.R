library(ggplot2)
library(ggpattern)
library(stringr)
library(tidyverse)
library(ggthemes)
library(patchwork)


custom_colors <- c(
  "asm_filter" = rgb(46, 108, 128, maxColorValue = 255),
  "asm_keep" = rgb(233, 127, 90, maxColorValue = 255)
)


get_theme <- function(size=12, angle=0) {
    defined_theme = theme_clean(base_size=size) + theme(legend.title=element_text(size=size), strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05)) #, legend.position="bottom", legend.box = "horizontal") 
    defined_theme
}


do_bar_chart <- function(data_path, out_path, threads, myparam) {
    # R code
    df_list = list()
    file_num = length(data_path[['union_count']])

    for (index in seq(file_num)) {
        union_count = read.table(data_path[['union_count']][index], colClasses = 'character', sep='\t')
        asm_union_count = read.table(data_path[['asm_union_count']][index], colClasses = 'character', sep='\t')
	metrics = as_tibble(union_count) %>% left_join(as_tibble(asm_union_count), by='V1')
	df_list[[index]] = metrics %>% mutate(cell_line = gsub('_hifi1_somatic_generation3_eval.tsv', '', basename(data_path[['union_count']][index])))
    }
    metrics = bind_rows(df_list)
    print(metrics)
    colnames(metrics)[1:3] = c('comb', 'caller',  'asm_keep')
    #metrics = data.frame(do.call(rbind, strsplit(metrics$comb, "")), metrics$caller %>% pull(), metrics$asm_keep)
    metrics = metrics %>% mutate(caller=as.numeric(caller), asm_keep=as.numeric(asm_keep))
    metrics = metrics %>% mutate(asm_filter=abs(caller-asm_keep)) %>% select(-caller)

    metrics = metrics %>% pivot_longer(cols=c("asm_filter", "asm_keep"), names_to="metrics", values_to = "SV_num")
    write_tsv(metrics, out_path[['stat']])

    pdf(out_path[['bar_pdf']], width=12, height=8)
    print(metrics)
    p1 = ggplot(data=metrics)+
      geom_bar_pattern(
          aes(
	        x=reorder(comb, SV_num), y=SV_num,
                fill=factor(metrics), 
		pattern=factor(metrics),
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
      facet_wrap(~cell_line, ncol=2) +
      get_theme(angle=0, size=9) +
      ggtitle("COLO829 tumor-normal paired somatic SV") + ylab("FP SV calls number")
    print(p1)
    dev.off()
}


do_bar_chart(snakemake@input, snakemake@output, snakemake@threads, snakemake@params)


