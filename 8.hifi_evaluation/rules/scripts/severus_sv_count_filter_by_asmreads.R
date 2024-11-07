library(ggplot2)
library(ggthemes)
library(stringr)
library(tidyverse)
library(patchwork)

custom_colors <- c(
  "asm_filter" = rgb(46, 108, 128, maxColorValue = 255),
  "asm_support" = rgb(233, 127, 90, maxColorValue = 255)
)

get_theme <- function(size=12, angle=0) {
    defined_theme = theme_clean(base_size=size) + theme(legend.title=element_text(size=size), strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05)) #, legend.position="bottom", legend.box = "horizontal") 
    defined_theme
}

do_bar_chart <- function(input, out_path, threads, myparam) {
    data_path1 = unlist(input[['severus_stat']])

    data_list = list()
    for (file in c(data_path1)) {
        data_list[[file]] = read_tsv(file, col_names=T)
        data_list[[file]] = data_list[[file]] %>% mutate(cell_line = gsub("_hifi1_somatic_generation2_filterasm.stat", "", basename(file)))
    }

    metrics = bind_rows(data_list)
    metrics = metrics %>% mutate(tool = ifelse(
      grepl("severus", cell_line),
      "severus", 
      ifelse(grepl("snf", cell_line),
             "snf",
             "msv"
      )
    ))
    write_tsv(metrics, out_path[['stat']])

    metrics = metrics %>% mutate(cell_line = gsub("msv_|snf_|severus_", "", cell_line))
    print(table(metrics$tool))

    metrics_summarize = list()
    metrics_summarize.before = list()
    for (cutoff in seq(3, 10)) {
        metrics_sub_severus = metrics %>% filter(read_name_number>=cutoff) %>% group_by(tool, cell_line) %>% summarize(total=n()) %>% ungroup()
        metrics_sub_severus_merge = metrics %>% filter(is_consensus) %>% filter(read_name_number>=cutoff) %>% group_by(tool, cell_line) %>% summarize(total=n()) %>% ungroup()

        metrics_sub_severus_dsa = metrics %>% filter(read_name_number>=cutoff) %>% filter(asm_support_onlyreadname>=cutoff) %>% group_by(tool, cell_line) %>% summarize(total=n()) %>% ungroup()
        metrics_sub_severus_merge_dsa = metrics %>% filter(is_consensus) %>% filter(read_name_number>=cutoff) %>% filter(asm_support_onlyreadname>=cutoff) %>% group_by(tool, cell_line) %>% summarize(total=n()) %>% ungroup()

	metrics_sub_severus_comb =  metrics_sub_severus %>% left_join(metrics_sub_severus_dsa, by=c("tool", "cell_line"))
	metrics_sub_severus_merge_comb = metrics_sub_severus_merge %>% left_join(metrics_sub_severus_merge_dsa, by=c("tool", "cell_line"))
	metrics_sub2 = bind_rows(metrics_sub_severus_comb %>% mutate(comb='severus+dsa'), metrics_sub_severus_merge_comb %>% mutate(comb='severus+merge+dsa'))
	metrics_sub2 = metrics_sub2 %>% mutate(asm_filter=total.x-total.y) %>% mutate(cutoff = cutoff) 

	metrics_sub_severus_comb = metrics_sub_severus_comb %>% mutate(asm_filter=total.x-total.y) %>% rename(asm_support=total.y) %>% select(tool, cell_line, asm_filter, asm_support) %>% pivot_longer(cols=starts_with("asm"), names_to="process", values_to="count")
	metrics_sub_severus_merge_comb = metrics_sub_severus_merge_comb %>% mutate(asm_filter=total.x-total.y) %>% rename(asm_support=total.y) %>% select(tool, cell_line, asm_filter, asm_support) %>% pivot_longer(cols=starts_with("asm"), names_to="process", values_to="count")
        metrics_sub_severus_comb = metrics_sub_severus_comb %>% mutate(cutoff=cutoff) %>% mutate(comb="severus+dsa")
        metrics_sub_severus_merge_comb = metrics_sub_severus_merge_comb %>% mutate(cutoff=cutoff) %>% mutate(comb="severus+merge+dsa")

	metrics_sub = bind_rows(metrics_sub_severus_comb, metrics_sub_severus_merge_comb)
        metrics_summarize[[cutoff]] = metrics_sub
        metrics_summarize.before[[cutoff]] = metrics_sub2
    }
    metrics_summarize.before = bind_rows(metrics_summarize.before)
    write_tsv(metrics_summarize.before, out_path[['test_stat_sv_num']])
    metrics_summarize = bind_rows(metrics_summarize)
    write_tsv(metrics_summarize, out_path[['stat_sv_num']])
    metrics_summarize = metrics_summarize %>% filter(cutoff > 2 & cutoff <= 6)
    ## tool    cell_line process     count cutoff comb

    pdf(out_path[['bar_pdf']], width=9, height=9)
    p = ggplot(data=metrics_summarize) + geom_bar(aes(x=comb, y=count, fill=factor(process)), stat = "identity", position = 'stack') + xlab("Self-assembly overlapped read cutoff") + ylab("SV numbers") + facet_grid(cell_line~cutoff, scales = "free") + get_theme(angle=30) + scale_fill_manual(values = custom_colors)
    print(p)
    dev.off()

    metrics_summarize = list()
    metrics_summarize.before = list()
    for (cutoff in seq(3, 10)) {
        metrics_sub_severus = metrics %>% filter(read_name_number>=cutoff) %>% group_by(tool, cell_line) %>% summarize(total=n()) %>% ungroup()
        metrics_sub_severus_merge = metrics %>% filter(is_consensus) %>% filter(read_name_number>=cutoff) %>% group_by(tool, cell_line) %>% summarize(total=n()) %>% ungroup()

        metrics_sub_severus_dsa = metrics %>% filter(read_name_number>=cutoff) %>% filter(asm_support>=cutoff) %>% group_by(tool, cell_line) %>% summarize(total=n()) %>% ungroup()
        metrics_sub_severus_merge_dsa = metrics %>% filter(is_consensus) %>% filter(read_name_number>=cutoff) %>% filter(asm_support>=cutoff) %>% group_by(tool, cell_line) %>% summarize(total=n()) %>% ungroup()

	metrics_sub_severus_comb =  metrics_sub_severus %>% left_join(metrics_sub_severus_dsa, by=c("tool", "cell_line"))
	metrics_sub_severus_merge_comb = metrics_sub_severus_merge %>% left_join(metrics_sub_severus_merge_dsa, by=c("tool", "cell_line"))
	metrics_sub2 = bind_rows(metrics_sub_severus_comb %>% mutate(comb='severus+dsa'), metrics_sub_severus_merge_comb %>% mutate(comb='severus+merge+dsa'))
	metrics_sub2 = metrics_sub2 %>% mutate(asm_filter=total.x-total.y) %>% mutate(cutoff = cutoff) 

	metrics_sub_severus_comb = metrics_sub_severus_comb %>% mutate(asm_filter=total.x-total.y) %>% rename(asm_support=total.y) %>% select(tool, cell_line, asm_filter, asm_support) %>% pivot_longer(cols=starts_with("asm"), names_to="process", values_to="count")
	metrics_sub_severus_merge_comb = metrics_sub_severus_merge_comb %>% mutate(asm_filter=total.x-total.y) %>% rename(asm_support=total.y) %>% select(tool, cell_line, asm_filter, asm_support) %>% pivot_longer(cols=starts_with("asm"), names_to="process", values_to="count")
        metrics_sub_severus_comb = metrics_sub_severus_comb %>% mutate(cutoff=cutoff) %>% mutate(comb="severus+dsa")
        metrics_sub_severus_merge_comb = metrics_sub_severus_merge_comb %>% mutate(cutoff=cutoff) %>% mutate(comb="severus+merge+dsa")

	metrics_sub = bind_rows(metrics_sub_severus_comb, metrics_sub_severus_merge_comb)
        metrics_summarize[[cutoff]] = metrics_sub
        metrics_summarize.before[[cutoff]] = metrics_sub2
    }
    metrics_summarize.before = bind_rows(metrics_summarize.before)
    write_tsv(metrics_summarize.before, out_path[['test_stat_sv_num']])
    metrics_summarize = bind_rows(metrics_summarize)
    write_tsv(metrics_summarize, out_path[['stat_sv_num']])
    metrics_summarize = metrics_summarize %>% filter(cutoff > 2 & cutoff <= 6)
    ## tool    cell_line process     count cutoff comb
    print(metrics_summarize)
    pdf(out_path[['bar_pdf2']], width=9, height=9)
    p = ggplot(data=metrics_summarize) + geom_bar(aes(x=comb, y=count, fill=factor(process)), stat = "identity", position = 'stack') + xlab("Self-assembly overlapped read cutoff") + ylab("SV numbers") + facet_grid(cell_line~cutoff, scales = "free") + get_theme(angle=30) + scale_fill_manual(values = custom_colors)
    print(p)
    dev.off()
}

do_bar_chart(snakemake@input, snakemake@output, snakemake@threads, snakemake@params)


