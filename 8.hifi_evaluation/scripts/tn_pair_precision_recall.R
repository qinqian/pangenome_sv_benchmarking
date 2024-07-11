library(ggplot2)
library(ggthemes)
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
                             tool=tools[-1], count=as.numeric(count), path=path)
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

    ggplot(metrics, aes(x=sensitivity, y=specificity, shape=factor(count))) + geom_point(aes(colour = factor(tool)), size = 5) + facet_wrap(~genome, scales='free')+theme_clean() + xlim(0, 1) + ylim(0, 1) + ylab('Precision') + xlab("Recall")
    ggsave(out_path[['plot']], width=15, height=5.5)
    write_tsv(metrics, out_path[['table']])


    precision <- seq(0.01, 1.0, by = 0.01)
    recall <- seq(0.01, 1.0, by = 0.01)
    # Create a data frame with all combinations of precision and recall
    grid <- expand.grid(Precision = precision, Recall = recall)
    # Compute the F1 score for each combination of precision and recall
    grid <- as_tibble(grid) %>%
        mutate(F1 = 2 * (Precision * Recall) / (Precision + Recall))

    defined_theme = theme_clean() + theme(legend.title=element_text(size=7), strip.text=element_text(size=7), legend.text=element_text(size=7), axis.title.x=element_text(size=7), axis.title.y=element_text(size=7), axis.text.y=element_text(size=7), axis.text.x=element_text(size=7))

    metrics_hg38 = metrics %>% filter(genome == 'hg38' & count >= 4)
    metrics_chm13 = metrics %>% filter(genome == "chm13" & count >= 4)

    # remove minisvl+tg, minisvl+ts
    metrics_hg38 = metrics_hg38 %>% filter(!(tool %in% c('minisvl+tg', 'minisvl+ts')))
    metrics_chm13 = metrics_chm13 %>% filter(!(tool %in% c('minisvl+tg', 'minisvl+ts', 'minisvl+g')))

    # Plot contours of constant F1 score
    ggplot() +
        geom_contour(data=grid, aes(x = Precision, y = Recall, z = F1), linetype="dashed", linewidth=0.35, color='gray', bins = 10) + 
        geom_point(data=metrics_hg38, aes(x=sensitivity, y=specificity, colour = factor(tool), size=count), alpha=0.6) + scale_size(range = c(1, 3)) + xlim(0, 1) + ylim(0, 1) + ylab('Precision') + xlab("Recall") +  
        geom_line(data=metrics_hg38, aes(x=sensitivity, y=specificity, colour = factor(tool))) + 
        defined_theme
    #ggplot(metrics_hg38, aes(x=sensitivity, y=specificity, shape=factor(count))) + geom_point(aes(colour = factor(tool)), size = 5)     
    ggsave(out_path[['hg38plot']], width=5.5, height=4.6)

    ggplot() +
        geom_contour(data=grid, aes(x = Precision, y = Recall, z = F1), linetype="dashed", linewidth=0.35, color='gray', bins = 10) + 
        geom_point(data=metrics_chm13, aes(x=sensitivity, y=specificity, colour = factor(tool), size=count), alpha=0.6) + scale_size(range = c(1, 3)) + xlim(0, 1) + ylim(0, 1) + ylab('Precision') + xlab("Recall") +  
        geom_line(data=metrics_chm13, aes(x=sensitivity, y=specificity, colour = factor(tool))) + 
        defined_theme
    #ggplot(metrics_chm13, aes(x=sensitivity, y=specificity, shape=factor(count))) + geom_point(aes(colour = factor(tool)), size = 5) + facet_wrap(~genome, scales='free')+xlim(0, 1) + ylim(0, 1) + ylab('Precision') + xlab("Recall") + defined_theme
    ggsave(out_path[['chm13plot']], width=5.5, height=4.6)

    metrics = metrics %>% group_by(tool, genome) %>% slice_max(order_by=f1,n=1, with_ties=F)
    ggplot(metrics, aes(y=f1, x=reorder(tool, f1), fill=tool)) + geom_bar(stat='identity') + facet_grid(~genome, scales='free')+theme_clean() + ylab('F1 score')
    ggsave(out_path[['f1plot']], width=12, height=5.5)
}

do_bar_chart(snakemake@input, snakemake@output, snakemake@threads, snakemake@params)
