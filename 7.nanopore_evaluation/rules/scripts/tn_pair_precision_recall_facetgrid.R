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
    print(dim(metrics))

    metrics$genome = ifelse(grepl("chm13", metrics$tool), "chm13", "hg38")
    print(metrics$tool)
    metrics$file = metrics$tool
    metrics$cell_line = str_extract(metrics$file, "(HG002rep1|HG002rep2)", group=1)
    metrics$tool = str_extract(metrics$tool, "(gafcall|minisv_pair|minisv_mosaic|nanomonsv|savana|severus|sniffles_mosaic|sniffles|svision|cutesv)", group=1)
    metrics$tool[metrics$tool=='gafcall'] = 'minisv'
    metrics$param = str_extract(metrics$file, "(l_l\\+t\\+g\\+s|l_l\\+x|g_g\\+x|l_l\\+g|l_l\\+t|l\\+tgs)", group=1)
    metrics$param[is.na(metrics$param)] = " "
    metrics$toolparam = paste0(metrics$tool, metrics$param)

    ggplot(metrics, aes(x=sensitivity, y=specificity, shape=factor(count))) + geom_point(aes(colour = factor(toolparam)), size = 3.5) + facet_grid(cell_line~genome, scales='free')+theme_classic() + xlim(0, 1) + ylim(0, 1) + ylab('Precision') + xlab("Recall")
    print(out_path)
    ggsave(out_path[['pr_plot']], width=11, height=9.5)

    write_tsv(metrics, out_path[['table1']])
    metrics = metrics %>% group_by(tool, cell_line, genome) %>% slice_max(order_by=f1,n=1, with_ties=F)
    write_tsv(metrics, out_path[['table2']])
    ggplot(metrics, aes(y=f1, x=reorder(tool, f1), fill=tool)) + geom_bar(stat='identity') + facet_grid(cell_line~genome, scales='free')+theme_classic() + ylab('F1 score')
    ggsave(out_path[['f1_plot']], width=8, height=5.5)
}

do_bar_chart(snakemake@input, snakemake@output, snakemake@threads, snakemake@params)
