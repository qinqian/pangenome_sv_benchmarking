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
        metrics = data.frame(sensitivity=sensitivity[-1], specificity=specificity[-1], 
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
    metrics$tool = str_extract(metrics$tool, "(gafcall|minisv_pair|minisv_mosaic|nanomonsv|savana|severus|sniffles_mosaic|sniffles|svision)", group=1)

    ggplot(metrics, aes(x=sensitivity, y=specificity, shape=factor(count))) + geom_point(aes(colour = factor(tool)), size = 5) + facet_wrap(~genome, scales='free')+theme_classic() + xlim(0, 1) + ylim(0, 1) + ylab('Precision') + xlab("Recall")
    print(out_path)
    ggsave(out_path[['plot']], width=15, height=5.5)
    write_tsv(metrics, out_path[['table']])
}

do_bar_chart(snakemake@input, snakemake@output, snakemake@threads, snakemake@params)
