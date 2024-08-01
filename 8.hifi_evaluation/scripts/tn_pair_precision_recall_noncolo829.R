library(ggplot2)
library(ggthemes)
library(stringr)
library(tidyverse)
library(patchwork)

custom_colors <- c(
  "minisvl+tg" = rgb(249, 134, 130, maxColorValue = 255),
  "minisvl+tgs" = rgb(187, 142, 33, maxColorValue = 255),
  "minisvl+gs" = rgb(187, 142, 33, maxColorValue = 255),
  "minisvl+g" = rgb(249, 134, 130, maxColorValue = 255),
  "nanomonsv" = rgb(147, 203, 118, maxColorValue = 255),
  "savana" = rgb(22, 183, 139, maxColorValue = 255),
  "severus" = rgb(16, 174, 228, maxColorValue = 255),
  "sniffles" = rgb(154, 130, 251, maxColorValue = 255),
  "svision" = rgb(248, 89, 206, maxColorValue = 255)
)

get_theme <- function(size=7, angle=0) {
    defined_theme = theme_clean(base_size=size) + theme(legend.title=element_text(size=size), strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05)) #, legend.position="bottom", legend.box = "horizontal") 
    defined_theme
}


load_metrics <- function(data_path) {
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
        metrics = data.frame(sensitivity=sensitivity[-1], precision=specificity[-1], f1=2*specificity[-1]*sensitivity[-1]/(sensitivity[-1]+specificity[-1]),
                             tool=tools[-1], count=as.numeric(count), path=path)
        print(metrics)
        df_list[[path]] = metrics
    }
    metrics = do.call(rbind, df_list)
    metrics$genome = ifelse(grepl("chm13", metrics$tool), "chm13", "hg38")
    metrics$file = metrics$tool
    metrics$tool = gsub("_pair", "", gsub("gafcall", "minisv", str_extract(metrics$tool, "(gafcall|minisv_pair|minisv_mosaic|nanomonsv|savana|severus|sniffles_mosaic|sniffles|svision)", group=1)))
    metrics$cell_line = str_extract(metrics$file, "(HCC1937|H1437|HCC1395|H2009|HCC1954|NCI1437|NCI2009)", group=1)
    print(table(metrics$cell_line))
   
    # previous gafcall output path name
    # different from current minisv 
    # for ultralong reads COLO829

    metrics$param = ifelse(grepl('minisv_pair|minisv_mosaic_asm', metrics$file), str_extract(metrics$file, "(l\\+t\\+g|l\\+g|l\\+t\\+g\\+s|l\\+x|g\\+x|l\\+t|l\\+tg)_(c2s0|mosaic)", group=1), str_extract(metrics$file, "(l\\+tg|l\\+g|l\\+gs|l\\+tgs|l\\+ts|l\\+tg)\\.pair", group=1))

    metrics$param[is.na(metrics$param)] = ""
    metrics$tool = paste0(metrics$tool, metrics$param)
    metrics
}


generate_grid <- function() {
    precision <- seq(0.01, 1.0, by = 0.01)
    recall <- seq(0.01, 1.0, by = 0.01)
    # Create a data frame with all combinations of precision and recall
    grid <- expand.grid(precision = precision, sensitivity=recall)
    # Compute the F1 score for each combination of precision and recall
    grid <- as_tibble(grid) %>%
        mutate(F1 = 2 * (precision * sensitivity) / (precision + sensitivity))
    grid
}


#plot_prec_recall <- function(grid, metrics) {
#    p = ggplot() +
#        geom_contour(data=grid, aes(x = Precision, y = Recall, z = F1), linetype="dashed", linewidth=0.45, color='gray', bins = 10) + 
#        geom_point(data=metrics, aes(x=sensitivity, y=specificity, colour = factor(tool), size=count), alpha=0.5) + scale_size(range = c(1, 4)) + xlim(0, 1) + ylim(0, 1) + ylab('Precision') + xlab("Recall") +  
#        geom_path(data=metrics[order(metrics$count),], aes(x=sensitivity, y=specificity, colour = factor(tool))) + 
#        #geom_line(data=metrics, aes(x=sensitivity, y=specificity, colour = factor(tool))) + 
#        get_theme()
#    p
#}

plot_f1 <- function(metrics) {
    metrics = metrics %>% group_by(tool, genome) %>% slice_max(order_by=f1,n=1, with_ties=F)
    ggplot(metrics, aes(y=f1, x=reorder(tool, f1), fill=tool)) + geom_bar(stat='identity') + facet_grid(~genome, scales='free')+theme_clean() + ylab('F1 score')
    ggsave(out_path[['f1plot']], width=12, height=5.5)
}

plot_prec_recall <- function(grid, metrics) {
    p = ggplot() +
        geom_contour(data=grid, aes(x = precision, y = sensitivity, z = F1), linetype="dashed", linewidth=0.45, color='gray', bins = 10) + 
        geom_point(data=metrics, aes(x=sensitivity, y=precision, colour = factor(tool), size=count), alpha=0.5) + scale_size(range = c(1, 4)) + xlim(0, 1) + ylim(0, 1) + ylab('Precision') + xlab("Recall") +  
        geom_path(data=metrics[order(metrics$count),], aes(x=sensitivity, y=precision, colour = factor(tool))) + 
        #geom_line(data=metrics, aes(x=sensitivity, y=precision, colour = factor(tool))) + 
        facet_wrap(cell_line~genome, nrow=5, ncol=2, scales = "free") + get_theme()
    p
    #ggplot(metrics, aes(x=sensitivity, y=precision, shape=factor(count))) + geom_point(aes(colour = factor(tool)), size = 5) + facet_grid(cell_line~genome, scales='free')+theme_bw() + xlim(0, 1) + ylim(0, 1) + ylab('Precision') + xlab("Recall") + get_theme()
}

do_bar_chart <- function(input, out_path, threads, myparam) {
    data_path = input[['mixed_noncolo829_hifi']]
    print(data_path)
    metrics = load_metrics(data_path)
    print(head(metrics))
    write_tsv(metrics, out_path[['mixed_table']])
    print(metrics)

    grid <- generate_grid()
    pdf(out_path[['mixedhg38plot']], width=5.5, height=5.6)
    mixed_hg38_p = plot_prec_recall(grid, metrics) + ggtitle("non-COLO829 HiFi tumor-normal 1:4 mixed reads")
    print(mixed_hg38_p)
    dev.off()
    #print((hifi.p1 + ont.p2 + mixed_hg38_p) + plot_layout(guides = "collect") +scale_colour_manual(values = custom_colors))
}

do_bar_chart(snakemake@input, snakemake@output, snakemake@threads, snakemake@params)
