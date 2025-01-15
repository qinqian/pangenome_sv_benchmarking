library(ggplot2)
library(ggthemes)
library(stringr)
library(tidyverse)
library(patchwork)


custom_colors <- c(
  "filtered by asm" = rgb(46, 108, 128, maxColorValue = 255),
  "kept by asm" = rgb(233, 127, 90, maxColorValue = 255)
)

custom_colors2 <- c(
  "filtered by asm" = "black",
  "kept by asm" = "gray"
)


get_theme <- function(size=12, angle=0) {
    defined_theme = theme_clean(base_size=size) + theme(legend.title=element_blank(), strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05), legend.position="bottom", legend.box = "horizontal")
    defined_theme
}

do_bar_chart <- function(input, out_path, threads, myparam) {
    data_path1 = unlist(input[['severus_stat']])
    data_path2 = unlist(input[['nano_stat']])
    data_path3 = unlist(input[['savana_stat']])
    data_path4 = unlist(input[['snf_stat']])
    data_path5 = unlist(input[['msv_stat']])
    data_path6 = unlist(input[['msv_lts_stat']])

    data_list5 = list()
    for (file in c(data_path5)) {
        data_list5[[file]] = read_tsv(file, col_names=T)
        data_list5[[file]] = data_list5[[file]] %>% mutate(cell_line = gsub("_hifi1_somatic_generation2_filterasm.stat", "", basename(file)))
    }
    data_list6 = list()
    for (file in c(data_path6)) {
        data_list6[[file]] = read_tsv(file, col_names=T)
        data_list6[[file]] = data_list6[[file]] %>% mutate(cell_line = gsub("_hifi1_somatic_generation2_filterasm.stat", "", basename(file)))
    }

    data_list1 = list()
    for (file in c(data_path1)) {
        data_list1[[file]] = read_tsv(file, col_names=T)
        data_list1[[file]] = data_list1[[file]] %>% mutate(cell_line = gsub("_hifi1_somatic_generation2_filterasm.stat", "", basename(file)))
    }
    data_list2 = list()
    for (file in c(data_path2)) {
        data_list2[[file]] = read_tsv(file, col_names=T)
        data_list2[[file]] = data_list2[[file]] %>% mutate(cell_line = gsub("_hifi1_somatic_generation2_filterasm.stat", "", basename(file)))
    }
    data_list3 = list()
    for (file in c(data_path3)) {
        data_list3[[file]] = read_tsv(file, col_names=T)
        data_list3[[file]] = data_list3[[file]] %>% mutate(cell_line = gsub("_hifi1_somatic_generation2_filterasm.stat", "", basename(file)))
    }
    data_list4 = list()
    for (file in c(data_path4)) {
        data_list4[[file]] = read_tsv(file, col_names=T)
        data_list4[[file]] = data_list4[[file]] %>% mutate(cell_line = gsub("_hifi1_somatic_generation2_filterasm.stat", "", basename(file)))
    }
    metrics1 = bind_rows(data_list1)
    metrics2 = bind_rows(data_list2)
    metrics3 = bind_rows(data_list3)
    metrics4 = bind_rows(data_list4)
    metrics5 = bind_rows(data_list5)
    metrics6 = bind_rows(data_list6)
    metrics5$svid = paste0('svid', metrics5$svid)
    metrics6$svid = paste0('svid', metrics6$svid)

    metrics = bind_rows(metrics1, metrics2, metrics3, metrics4, metrics5, metrics6)
    metrics = metrics %>% mutate(tool = ifelse(
      grepl("severus", cell_line),
      "severus", 
      ifelse(grepl("snf", cell_line),
             "snf",
             ifelse(grepl("nano", cell_line),
                   "nanomonsv",
                  ifelse(grepl("sava", cell_line), 
                          "savana",
                  ifelse(grepl("msv_ltgs", cell_line),
                         "msv_ltgs",
                         "msv_lts")
      )
    ))))
    
    print(table(metrics$tool))
    metrics$cell_line = gsub('ltgs_', '', metrics$cell_line)
    metrics$cell_line = gsub('lts_', '', metrics$cell_line)
    write_tsv(metrics, out_path[['stat']])

    metrics = metrics %>% mutate(cell_line = gsub("msv_|snf_|severus_|nanomonsv_|savana_", "", cell_line))
    metrics_summarize = list()
    metrics_summarize.before = list()
    for (cutoff in seq(3, 10)) {
        metrics_sub_all = metrics %>% filter(read_name_number>=cutoff) %>% group_by(tool, cell_line) %>% summarize(total=n()) %>% ungroup()
        metrics_sub_all_dsa = metrics %>% filter(read_name_number>=cutoff) %>% filter(asm_support_onlyreadname>=cutoff) %>% group_by(tool, cell_line) %>% summarize(total=n()) %>% ungroup()

	metrics_sub_comb = metrics_sub_all %>% left_join(metrics_sub_all_dsa, by=c("tool", "cell_line"))
	metrics_sub_comb = metrics_sub_comb %>% mutate(asm_filter=total.x-total.y) %>% mutate(cutoff = cutoff) 
	print(head(metrics_sub_comb))
	metrics_sub_comb = metrics_sub_comb %>% mutate(asm_filter=total.x-total.y) %>% rename(asm_support=total.y) %>% select(tool, cell_line, asm_filter, asm_support, cutoff) %>% pivot_longer(cols=starts_with("asm"), names_to="process", values_to="count")
        metrics_summarize[[cutoff]] = metrics_sub_comb
    }
    metrics_summarize = bind_rows(metrics_summarize)
    print(table(metrics_summarize$tool))
    write_tsv(metrics_summarize, out_path[['stat_sv_num']])
    metrics_summarize = metrics_summarize %>% filter(cutoff >= 3 & cutoff <= 7)
    ##metrics_summarize = metrics_summarize %>% filter(tool != "msv_lts")

    metrics_summarize = metrics_summarize %>% mutate(
        tools = case_when(
            tool=="msv_ltgs" ~ "minisv_ltg",
            tool=="msv_lts" ~ "minisv_lt",
            tool=="snf" ~ "sniffles2", 
            .default = tool
        )
    )

    pdf(out_path[['bar_pdf']], width=9, height=9)
    p = ggplot(data=metrics_summarize, aes(x=tools, y=count, fill=factor(process))) + geom_bar(position="stack", stat="identity") + xlab("Self-assembly overlapped read cutoff") + ylab("SV numbers") + facet_grid(cell_line~cutoff, scales = "free") + get_theme(angle=30) + scale_fill_manual(values = custom_colors)
    print(p)
    dev.off()

    metrics = metrics %>% mutate(cell_line = gsub("msv_lts_|msv_ltgs_|snf_|severus_|nanomonsv_|savana_", "", cell_line))
    metrics_summarize = list()
    metrics_summarize.before = list()
    for (cutoff in seq(3, 10)) {
        metrics_sub_all = metrics %>% filter(read_name_number>=cutoff) %>% group_by(tool, cell_line) %>% summarize(total=n()) %>% ungroup()
        metrics_sub_all_dsa = metrics %>% filter(read_name_number>=cutoff) %>% filter(asm_support>=cutoff) %>% group_by(tool, cell_line) %>% summarize(total=n()) %>% ungroup()

	metrics_sub_comb = metrics_sub_all %>% left_join(metrics_sub_all_dsa, by=c("tool", "cell_line"))
	metrics_sub_comb = metrics_sub_comb %>% mutate(asm_filter=total.x-total.y) %>% mutate(cutoff = cutoff) 
	print(head(metrics_sub_comb))
	metrics_sub_comb = metrics_sub_comb %>% mutate(asm_filter=total.x-total.y) %>% rename(asm_support=total.y) %>% select(tool, cell_line, asm_filter, asm_support, cutoff) %>% pivot_longer(cols=starts_with("asm"), names_to="process", values_to="count")
        metrics_summarize[[cutoff]] = metrics_sub_comb
    }

    metrics_summarize = bind_rows(metrics_summarize)
    write_tsv(metrics_summarize, out_path[['stat_sv_num']])

    #metrics_summarize = metrics_summarize %>% filter(tool != "msv_lts")
    metrics_summarize = metrics_summarize %>% filter(cutoff >= 3 & cutoff <= 7)

    metrics_summarize = metrics_summarize %>% mutate(
        tools = case_when(
            tool=="msv_ltgs" ~ "minisv_ltg",
            tool=="snf" ~ "sniffles2", 
            tool=="msv_lts" ~ "minisv_lt",
            .default = tool
        )
    )

    metrics_summarize = metrics_summarize %>% mutate(
          process=case_when(
              process == "asm_support" ~ "kept by asm",
              process == "asm_filter" ~ "filtered by asm",
              .default = process,
          )
    )

    pdf(out_path[['bar_pdf2']], width=12, height=8)
    p = ggplot(data=metrics_summarize, aes(x=tools, y=count, fill=factor(process))) + geom_bar(position="stack", stat="identity") + xlab("Self-assembly overlapped read cutoff") + ylab("SV numbers") + facet_grid(cell_line~cutoff, scales = "free") + get_theme(angle=30) + scale_fill_manual(values = custom_colors)
    print(p)
    dev.off()

    pdf(out_path[['bar_pdf2_bw']], width=12, height=8)
    p = ggplot(data=metrics_summarize, aes(x=tools, y=count, fill=factor(process))) + geom_bar(position="stack", stat="identity") + xlab("Self-assembly overlapped read cutoff") + ylab("SV numbers") + facet_grid(cell_line~cutoff, scales = "free") + get_theme(angle=30) + scale_fill_manual(values = custom_colors2)
    print(p)
    dev.off()
}


do_bar_chart(snakemake@input, snakemake@output, snakemake@threads, snakemake@params)


