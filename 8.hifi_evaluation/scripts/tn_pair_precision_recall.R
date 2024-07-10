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
    #metrics$tool = gsub("msv|msv\\.gz|vcf|vcf\\.gz", "", gsub("output/", "", gsub("chm13|grch38", "", metrics$tool)))
    #metrics$tool = gsub("../1a.alignment_sv_tools/", "", gsub("chm13|grch38", "", metrics$tool))
    print(metrics$tool)
    print(str_extract(metrics$tool, "(gafcall|minisv_pair|minisv_mosaic|nanomonsv|savana|severus|sniffles_mosaic|sniffles|svision)", group=1))
    metrics$file = metrics$tool
    metrics$tool = gsub("gafcall", "minisv", str_extract(metrics$tool, "(gafcall|minisv_pair|minisv_mosaic|nanomonsv|savana|severus|sniffles_mosaic|sniffles|svision)", group=1))
    metrics$param = str_extract(metrics$file, "(l\\+g|l\\+gs|l\\+tgs|l\\+ts|l\\+tg|l\\+tgs)\\.pair", group=1)
    metrics$param[is.na(metrics$param)] = " "
    metrics$tool = paste0(metrics$tool, metrics$param)

    ggplot(metrics, aes(x=sensitivity, y=specificity, shape=factor(count))) + geom_point(aes(colour = factor(tool)), size = 5) + facet_wrap(~genome, scales='free')+theme_classic() + xlim(0, 1) + ylim(0, 1) + ylab('Precision') + xlab("Recall")
    ggsave(out_path[['plot']], width=15, height=5.5)
    write_tsv(metrics, out_path[['table']])

    metrics_hg38 = metrics %>% filter(genome == 'hg38')
    metrics_chm13 = metrics %>% filter(genome == "chm13")
    ggplot(metrics_hg38, aes(x=sensitivity, y=specificity, shape=factor(count))) + geom_point(aes(colour = factor(tool)), size = 5) + facet_wrap(~genome, scales='free')+theme_classic() + xlim(0, 1) + ylim(0, 1) + ylab('Precision') + xlab("Recall")
    ggsave(out_path[['hg38plot']], width=9.5, height=2.5)
    ggplot(metrics_chm13, aes(x=sensitivity, y=specificity, shape=factor(count))) + geom_point(aes(colour = factor(tool)), size = 5) + facet_wrap(~genome, scales='free')+theme_classic() + xlim(0, 1) + ylim(0, 1) + ylab('Precision') + xlab("Recall")
    ggsave(out_path[['chm13plot']], width=9.5, height=2.5)

    metrics = metrics %>% group_by(tool, genome) %>% slice_max(order_by=f1,n=1, with_ties=F)
    ggplot(metrics, aes(y=f1, x=reorder(tool, f1), fill=tool)) + geom_bar(stat='identity') + facet_grid(~genome, scales='free')+theme_classic() + ylab('F1 score')
    ggsave(out_path[['f1plot']], width=15, height=5.5)
}

do_bar_chart(snakemake@input, snakemake@output, snakemake@threads, snakemake@params)
