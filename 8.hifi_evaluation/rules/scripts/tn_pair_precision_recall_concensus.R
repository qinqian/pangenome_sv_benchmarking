library(ggplot2)
library(stringr)
library(tidyverse)
library(tidytext)
library(ggthemes)
library(patchwork)

get_theme <- function(size=7, angle=0) {
    defined_theme = theme_bw(base_size=size) + theme(legend.title=element_text(size=size), strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05)) #, legend.position="bottom", legend.box = "horizontal") 
    defined_theme
}

generate_grid <- function() {
    precision <- seq(0.01, 1.0, by = 0.01)
    recall <- seq(0.01, 1.0, by = 0.01)
    # Create a data frame with all combinations of precision and recall
    grid <- expand.grid(precision = precision, sensitivity = recall)
    # Compute the F1 score for each combination of precision and recall
    grid <- as_tibble(grid) %>%
        mutate(F1 = 2 * (precision * sensitivity) / (precision + sensitivity))
    grid
}

plot_prec_recall <- function(grid, metrics) {
    p = ggplot() +
        geom_contour(data=grid, aes(x = precision, y = sensitivity, z = F1), linetype="dashed", linewidth=0.45, color='gray', bins = 10) + 
        geom_point(data=metrics, aes(x=sensitivity, y=precision, colour = factor(tool), size=count), alpha=0.5) + scale_size_continuous(name = "Count", breaks = c(3, 4, 5, 10), range = c(1, 5)) + xlim(0, 1) + ylim(0, 1) + ylab('Precision') + xlab("Recall") +  
        geom_path(data=metrics[order(metrics$count),], aes(x=sensitivity, y=precision, colour = factor(tool))) + 
        #geom_line(data=metrics, aes(x=sensitivity, y=precision, colour = factor(tool))) + 
        facet_wrap(cell_line~genome, nrow=5, ncol=2, scales = "free") + get_theme()
    p
    #ggplot(metrics, aes(x=sensitivity, y=precision, shape=factor(count))) + geom_point(aes(colour = factor(tool)), size = 5) + facet_grid(cell_line~genome, scales='free')+theme_bw() + xlim(0, 1) + ylim(0, 1) + ylab('Precision') + xlab("Recall") + get_theme()
}


parse_list_metrics_files <- function(data_path) {
    df_list = list()
    for (path in unlist(data_path)) {
        res = read_tsv(path, col_names = FALSE)
        print(res)
        res = res %>% mutate(X4=1-X4) %>% select(X1, X4, X5)
        res = res %>% pivot_wider(names_from = X1, values_from=X4, values_fill=NA)
        colnames(res) = c("file", "sensitivity", "precision")
        count = str_extract(path, "count(\\d+)", group=1)
        res$count = count
        res$eval_file = path
        res = res %>% mutate(f1=2*sensitivity*precision/(sensitivity+precision))
        df_list[[path]] = res
    }

    metrics = do.call(rbind, df_list)
    metrics$count = as.numeric(metrics$count)
 
    # start from 3 read counts
    metrics = metrics %>% filter(count >= 3)

    metrics$genome = ifelse(grepl("chm13", metrics$eval_file), "chm13", "hg38")
    metrics$cell_line = str_extract(metrics$eval_file, "(COLO829|HCC1937|H1437|HCC1395|H2009|HCC1954|NCI1437|NCI2009)", group=1)
  
    # Exclude COLO829 from the analysis
    metrics = metrics %>% filter(!grepl("COLO829", metrics$cell_line))

    #NOTE: is H2009 equal to NCI2009?
    #NOTE: is H1437 equal to NCI1437?
    metrics$cell_line = gsub("H2009", "NCI2009",  metrics$cell_line)
    metrics$cell_line = gsub("H1437", "NCI1437",  metrics$cell_line)

    metrics$tool = str_extract(metrics$file, "(gafcall|minisv_pair|minisv_mosaic|nanomonsv|savana|severus\\/|sniffles_mosaic|sniffles|svision|cutesv)", group=1)
#l_l+t+g+s
    metrics$tool = gsub('gafcall', 'minisv', metrics$tool)
    print(table(metrics$tool))

    metrics$param = str_extract(metrics$file, "(l\\+g|l_l\\+t\\+g\\+s|l_l\\+t\\+g\\+s|l_l\\+x|g_g\\+x|l_l\\+g|l_l\\+t|l\\+tgs|l\\+gs|l\\+ts|l\\+tg)(_c2s0|\\.pair-c2s0)", group=1)
    metrics$param[is.na(metrics$param)] = ""
    metrics$tool = paste0(metrics$tool, metrics$param)

    metrics$tool = gsub('minisv_pairl_l\\+t\\+g\\+s', 'minisvl\\+tgs', metrics$tool)
    metrics$tool = ifelse(metrics$genome == 'chm13' & (grepl("l\\+gs", metrics$tool)), gsub("l\\+gs", "l\\+tgs", metrics$tool), metrics$tool)
    metrics$tool = gsub('minisvl\\+tgs', 'msv:tgs', metrics$tool)

    metrics = metrics %>% filter(!(tool %in% c("minisvl+g", "minisvl+tg", "minisvl+ts")))
    print(table(metrics$tool))
    metrics
}


do_bar_chart <- function(data_path, out_path, threads, myparam) {
    # R code
    metrics = parse_list_metrics_files(data_path[['ont']])
    metrics_hifi = parse_list_metrics_files(data_path[['hifi']])
    write_tsv(metrics, out_path[['table1']])
    write_tsv(metrics_hifi, out_path[['table1_hifi']])

    grid = generate_grid()
    hifi.p = plot_prec_recall(grid, metrics_hifi)
    ggsave(out_path[['plot1_hifi']], width=9, height=13.5)

    ont.p = plot_prec_recall(grid, metrics)
    ggsave(out_path[['plot1']], width=9, height=13.5)

    pdf(out_path[['joint_plot']], width=18, height=13.5)
    print(hifi.p + ont.p)
    dev.off()

    metrics = metrics %>% group_by(tool, cell_line, genome) %>% slice_max(order_by=f1, n=1, with_ties=F)
    write_tsv(metrics, out_path[['table2']])
 
    metrics %>% 
        #ggplot(aes(y=f1, x=reorder_within(tool, f1, list(genome, cell_line)), fill=tool)) + geom_bar(stat='identity') + scale_x_reordered() + facet_wrap(genome~cell_line, ncol=5)+theme_bw() + ylab('F1 score') + theme(axis.text.x = element_text(angle=90, size=12, hjust=1))
        #ggplot(aes(y=f1, x=reorder_within(tool, f1, list(genome, cell_line)), fill=tool)) + geom_bar(stat='identity') + scale_x_reordered() + facet_wrap(genome~cell_line, scales='free_x', ncol=5)+theme_bw() + ylab('F1 score') + theme(axis.text.x = element_text(angle=90, size=12, hjust=1))
        ggplot(aes(y=f1, x=reorder(tool, f1), fill=tool)) + geom_bar(stat='identity') + facet_grid(genome~cell_line, scales='free_x') + theme_bw() + ylab('F1 score') + theme(axis.text.x = element_text(angle=90, size=12, hjust=1))
        
    ggsave(out_path[['plot2']], width=12, height=8)
}

do_bar_chart(snakemake@input, snakemake@output, snakemake@threads, snakemake@params)
