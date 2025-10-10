library(ggplot2)
library(stringr)
library(tidyverse)

do_bar_chart <- function(data_path, out_path, threads, myparam) {
    # R code
    df_list = list()
    for (path in unlist(data_path)) {
        res = read_tsv(path, col_names = FALSE)
        print(res)
        # sensitivity
        sensitivity = res %>% select(X2) %>% pull()
        # precision
        specificity = unlist(as.data.frame(res)[1, 2:(ncol(res)-1)])
        count = str_extract(path, "count(\\d+)", group=1)
        if (is.na(count)) {
           count = str_extract(path, "cutoff(\\d+)", group=1)
        }
        tools = res[, ncol(res)] %>% pull()
        metrics = data.frame(sensitivity=sensitivity[-1], specificity=specificity[-1], f1=2*specificity[-1]*sensitivity[-1]/(sensitivity[-1]+specificity[-1]),
                             tool=tools[-1], count=count, path=path)
        print(metrics)
        df_list[[path]] = metrics
    }

    metrics = do.call(rbind, df_list)
    #print(metrics)
    print(dim(metrics))
    #print(as_tibble(metrics) %>% pivot_longer(names_to="tool", values_to = "count"))

    metrics$genome = ifelse(grepl("chm13", metrics$tool), "chm13", "hg38")
    metrics$file = metrics$tool
    metrics$cell_line = str_extract(metrics$tool, "(COLO829)[T_cgo\\/]", group=1)
    metrics$cell_line2 = str_extract(metrics$tool, "(ont1|ont2)", group=1)
    metrics$cell_line2[is.na(metrics$cell_line2)] = ''
    metrics$cell_line = paste0(metrics$cell_line, metrics$cell_line2)

    metrics$tool = str_extract(metrics$tool, "(minisv|gafcall|minisv_pair|minisv_mosaic|nanomonsv|savana|severus|sniffles_mosaic|sniffles|svision|cutesv)", group=1)
    metrics$tool[metrics$tool=='gafcall'] = 'minisv'
    metrics$param = str_extract(metrics$file, "(l_t\\+s|l_l\\+\\g\\+s|l_t\\+g|l_t\\+g\\+s|l\\+s|l_l\\+t\\+g\\+s|l_l\\+x|g_g\\+x|l_l\\+g|l_l\\+t|l\\+tgs)_merge", group=1)
    metrics$param[is.na(metrics$param)] = " "
    metrics$tool = paste0(metrics$tool, metrics$param)

    write_tsv(metrics, out_path[['table1']])
    ggplot(metrics, aes(x=sensitivity, y=specificity, shape=factor(count))) + geom_point(aes(colour = factor(tool)), size = 5) + facet_wrap(genome~cell_line, scales='free')+theme_classic() + xlim(0, 1) + ylim(0, 1) + ylab('Precision') + xlab("Recall")
    ggsave(out_path[['plot1']], width=16, height=10)

    metrics = metrics %>% group_by(tool, cell_line, genome) %>% slice_max(order_by=f1, n=1, with_ties=F)
    write_tsv(metrics, out_path[['table2']])
    ggplot(metrics, aes(y=f1, x=reorder(tool, f1), fill=tool)) + geom_bar(stat='identity') + facet_grid(genome~cell_line, scales='free')+theme_classic() + ylab('F1 score') + theme(axis.text.x = element_text(angle=90, size=12, hjust=1))
    ggsave(out_path[['plot2']], width=9.8, height=5.5)
}

do_bar_chart(snakemake@input, snakemake@output, snakemake@threads, snakemake@params)
