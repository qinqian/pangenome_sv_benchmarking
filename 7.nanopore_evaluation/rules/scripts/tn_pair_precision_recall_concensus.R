library(ggplot2)
library(stringr)
library(tidyverse)
library(tidytext)

do_bar_chart <- function(data_path, out_path, threads, myparam) {
    # R code
    df_list = list()
    for (path in unlist(data_path)) {
        res = read_tsv(path, col_names = FALSE)
        res = res %>% mutate(X4=1-X4) %>% select(X1, X4, X5)
        res = res %>% pivot_wider(names_from = X1, values_from=X4, values_fill=NA)
        colnames(res) = c("file", "sensitivity", "precision")
        count = str_extract(path, "count(\\d+)", group=1)
        res$count = count
        res$eval_file = path
        ##    X1       X2    X3     X4 X5
        ##   <chr> <dbl> <dbl>  <dbl> <chr>
        ## 1 RN      200    56 0.28   output/severus/H1437_ont1/chm13/somatic_SVs/severus…
        ## 2 RP      235    35 0.149  output/severus/H1437_ont1/chm13/somatic_SVs/severus…
        ## 3 RN      196    63 0.321  output/minisv_pair/H1437_ont1_pair_chm13l_l+t+g+s_c…
        ## 4 RP     8746  8553 0.978  output/minisv_pair/H1437_ont1_pair_chm13l_l+t+g+s_c…
        res = res %>% mutate(f1=2*sensitivity*precision/(sensitivity+precision))
        df_list[[path]] = res
    }

    metrics = do.call(rbind, df_list)
    print(dim(metrics))
    metrics$genome = ifelse(grepl("chm13", metrics$eval_file), "chm13", "hg38")
    metrics$cell_line = str_extract(metrics$eval_file, "(HCC1937|H1437|HCC1395|H2009|HCC1954)", group=1)
    print(metrics)

    metrics$tool = str_extract(metrics$file, "(gafcall|minisv_pair|minisv_mosaic|nanomonsv|savana|severus|sniffles_mosaic|sniffles|svision|cutesv)", group=1)
#l_l+t+g+s
    metrics$param = str_extract(metrics$file, "(l_l\\+t\\+g\\+s|l_l\\+t\\+g\\+s|l_l\\+x|g_g\\+x|l_l\\+g|l_l\\+t|l\\+tgs)_merge", group=1)
    metrics$param[is.na(metrics$param)] = " "
    metrics$tool = paste0(metrics$tool, metrics$param)

    write_tsv(metrics, out_path[['table1']])
    ggplot(metrics, aes(x=sensitivity, y=precision, shape=factor(count))) + geom_point(aes(colour = factor(tool)), size = 5) + facet_grid(cell_line~genome, scales='free')+theme_bw() + xlim(0, 1) + ylim(0, 1) + ylab('Precision') + xlab("Recall")
    ggsave(out_path[['plot1']], width=12, height=9.5)

    metrics = metrics %>% group_by(tool, cell_line, genome) %>% slice_max(order_by=f1, n=1, with_ties=F)
    write_tsv(metrics, out_path[['table2']])
 
    metrics %>% 
        #ggplot(aes(y=f1, x=reorder_within(tool, f1, list(genome, cell_line)), fill=tool)) + geom_bar(stat='identity') + scale_x_reordered() + facet_wrap(genome~cell_line, ncol=5)+theme_bw() + ylab('F1 score') + theme(axis.text.x = element_text(angle=90, size=12, hjust=1))
        #ggplot(aes(y=f1, x=reorder_within(tool, f1, list(genome, cell_line)), fill=tool)) + geom_bar(stat='identity') + scale_x_reordered() + facet_wrap(genome~cell_line, scales='free_x', ncol=5)+theme_bw() + ylab('F1 score') + theme(axis.text.x = element_text(angle=90, size=12, hjust=1))
        ggplot(aes(y=f1, x=reorder(tool, f1), fill=tool)) + geom_bar(stat='identity') + facet_grid(genome~cell_line, scales='free_x') + theme_bw() + ylab('F1 score') + theme(axis.text.x = element_text(angle=90, size=12, hjust=1))
        
    ggsave(out_path[['plot2']], width=12, height=8)
}

do_bar_chart(snakemake@input, snakemake@output, snakemake@threads, snakemake@params)
