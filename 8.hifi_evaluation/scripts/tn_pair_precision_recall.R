library(ggplot2)
library(ggthemes)
library(stringr)
library(tidyverse)
library(patchwork)

custom_colors <- c(
  "minisvl+tg" = rgb(249, 134, 130, maxColorValue = 255),
  "minisvl+tgs" = rgb(187, 142, 33, maxColorValue = 255),
  "msv:tgs" = rgb(187, 142, 33, maxColorValue = 255),
  "minisvl+gs" = rgb(187, 142, 33, maxColorValue = 255),
  "msv:gs" = rgb(187, 142, 33, maxColorValue = 255),
  "minisvl+g" = rgb(249, 134, 130, maxColorValue = 255),
  "msv:g" = rgb(249, 134, 130, maxColorValue = 255),
  "msv:tg" = rgb(249, 134, 130, maxColorValue = 255),
  "nanomonsv" = rgb(147, 203, 118, maxColorValue = 255),
  "savana" = rgb(22, 183, 139, maxColorValue = 255),
  "severus" = rgb(16, 174, 228, maxColorValue = 255),
  "sniffles" = rgb(154, 130, 251, maxColorValue = 255),
  "svision" = rgb(248, 89, 206, maxColorValue = 255)
)

get_theme <- function(size=7, angle=0) {
    defined_theme = theme_bw(base_size=size) + theme(legend.title=element_text(size=size), strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05)) #, legend.position="bottom", legend.box = "horizontal") 
    defined_theme
}


load_metrics <- function(data_path) {
    df_list = list()
    for (path in unlist(data_path)) {
        res = read_tsv(path, col_names = FALSE)
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
        df_list[[path]] = metrics
    }
    metrics = do.call(rbind, df_list)
    metrics$genome = ifelse(grepl("chm13", metrics$tool), "chm13", "hg38")
    metrics$file = metrics$tool
    metrics$tool = gsub("_pair", "", gsub("gafcall", "minisv", str_extract(metrics$tool, "(gafcall|minisv_pair|minisv_mosaic|nanomonsv|savana|severus|sniffles_mosaic|sniffles|svision)", group=1)))
   
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
    grid <- expand.grid(Precision = precision, Recall = recall)
    # Compute the F1 score for each combination of precision and recall
    grid <- as_tibble(grid) %>%
        mutate(F1 = 2 * (Precision * Recall) / (Precision + Recall))
    grid
}


plot_prec_recall <- function(grid, metrics) {
    p = ggplot() +
        #geom_contour(data=grid, aes(x = Precision, y = Recall, z = F1), linetype="dashed", linewidth=0.45, color='gray', bins = 10) + 
        geom_contour(data=grid, aes(x = Recall, y = Precision, z = F1), linetype="dashed", linewidth=0.45, color='gray', bins = 10) + 
        geom_point(data=metrics, aes(x=sensitivity, y=specificity, colour = factor(tool), size=count), alpha=0.5) + ylab('Precision') + xlab("Recall") + scale_size_continuous(name = "Count", breaks = c(3, 4, 5, 10), range = c(1, 5)) +  #scale_size(range = c(1, 4)) + + xlim(0, 1) + ylim(0, 1) 
        geom_path(data=metrics[order(metrics$count),], aes(x=sensitivity, y=specificity, colour = factor(tool))) + 
        get_theme()
    p
}

plot_f1 <- function(metrics) {
    metrics = metrics %>% group_by(tool, genome) %>% slice_max(order_by=f1,n=1, with_ties=F)
    ggplot(metrics, aes(y=f1, x=reorder(tool, f1), fill=tool)) + geom_bar(stat='identity') + facet_grid(~genome, scales='free')+theme_clean() + ylab('F1 score')
    ggsave(out_path[['f1plot']], width=12, height=5.5)
}

do_bar_chart <- function(input, out_path, threads, myparam) {
    data_path = input[['colo829_hifi']]
    data_path_mixed = input[['colo829_hifi_mosaic']]
    data_path_mixed_100kb = input[['colo829_hifi_mosaic_100kb']]
    data_path_ont = input[['colo829_ont']]

    metrics = load_metrics(data_path)
    metrics_hg38 = metrics %>% filter(genome == 'hg38' & count >= 2)
    metrics_chm13 = metrics %>% filter(genome == "chm13" & count >= 2)
    write_tsv(metrics, out_path[['table']])

    metrics.ont = load_metrics(data_path_ont)
    metrics_hg38.ont  = metrics.ont %>% filter(genome == 'hg38' & count >= 2)
    metrics_chm13.ont = metrics.ont %>% filter(genome == "chm13" & count >= 2)
    write_tsv(metrics.ont, out_path[['ont_table']])

    metrics.mixed = load_metrics(data_path_mixed)
    metrics_hg38.mixed  = metrics.mixed %>% filter(genome == 'hg38' & count >= 2)
    #metrics_chm13.mixed = metrics.mixed %>% filter(genome == "chm13" & count >= 4)
    write_tsv(metrics.mixed, out_path[['mixed_table']])

    print(data_path_mixed_100kb)
    plot_100kb = F
    if (length(unlist(data_path_mixed_100kb)) > 0) {
        plot_100kb = T
        metrics.mixed_100kb = load_metrics(data_path_mixed_100kb)
        metrics.hg38_mixed_100kb = metrics.mixed_100kb %>% filter(genome == 'hg38' & count >= 2)
        write_tsv(metrics.hg38_mixed_100kb, out_path[['mixed_table_100kb']])
    }

    # remove minisvl+tg, minisvl+ts
    metrics_hg38 = metrics_hg38 %>% filter(!(tool %in% c('minisvl+ts')))
    metrics_chm13 = metrics_chm13 %>% filter(!(tool %in% c('minisvl+ts')))

    metrics_hg38.ont = metrics_hg38.ont %>% filter(!(tool %in% c('minisvl+ts', 'minisvl+x', 'minisvl+g', 'minisvg+x', 'minisvl+t')))
    metrics_hg38.ont$tool = gsub("minisvl\\+t\\+g\\+s", "minisvl\\+tgs", metrics_hg38.ont$tool)

    metrics_hg38.mixed$tool = gsub("minisvl\\+t\\+g\\+s", "minisvl\\+tgs", gsub("_mosaic", "", gsub("_mosaic ", "", metrics_hg38.mixed$tool)))
    metrics_hg38.mixed$tool = gsub("minisvl\\+t\\+g", "minisvl\\+tg", metrics_hg38.mixed$tool)

    metrics_chm13.ont = metrics_chm13.ont %>% filter(!(tool %in% c('minisvl+ts', 'minisvl+x', 'minisvg+x')))
    metrics_chm13.ont$tool = gsub("minisvl\\+t\\+g\\+s", "minisvl\\+gs", metrics_chm13.ont$tool)
    metrics_chm13$tool = gsub("minisvl\\+", "msv:", metrics_chm13$tool)
    metrics_chm13.ont$tool = gsub("minisvl\\+", "msv:", metrics_chm13.ont$tool)
    metrics_hg38$tool = gsub("minisvl\\+", "msv:", metrics_hg38$tool)
    metrics_hg38.ont$tool = gsub("minisvl\\+", "msv:", metrics_hg38.ont$tool)
    metrics_hg38.mixed$tool = gsub("minisvl\\+", "msv:", metrics_hg38.mixed$tool)
    
    if (plot_100kb) {
        metrics.hg38_mixed_100kb$tool = gsub("minisvl\\+t\\+g\\+s", "msv:tgs", gsub("_mosaic", "", metrics.hg38_mixed_100kb$tool))
        metrics.hg38_mixed_100kb$tool = gsub("minisvl\\+t\\+g", "msv:tg", metrics.hg38_mixed_100kb$tool)
    }

    grid <- generate_grid()
    pdf(out_path[['mixedhg38plot']], width=5.5, height=5.6)
    mixed_hg38_p = plot_prec_recall(grid, metrics_hg38.mixed) + ggtitle("COLO829 HiFi tumor-normal 1:4 mixed reads")+scale_colour_manual(values = custom_colors)
    print(mixed_hg38_p)
    dev.off()

    # Plot contours of constant F1 score
    if (plot_100kb) {
        pdf(out_path[['hg38plot']], width=8.5, height=7.6)
        hifi.p1 = plot_prec_recall(grid, metrics_hg38) + ggtitle("COLO829 HiFi tumor-normal pair")
        ont.p2 = plot_prec_recall(grid, metrics_hg38.ont) + ggtitle("COLO829 ONT tumor-normal pair")
        mixed_hg38_p_100kb = plot_prec_recall(grid, metrics.hg38_mixed_100kb) + ggtitle("COLO829 HiFi tumor-normal 1:4 mixed reads >100kb SV")
        print((hifi.p1 + ont.p2 +scale_colour_manual(values = custom_colors)) / (mixed_hg38_p + mixed_hg38_p_100kb+scale_colour_manual(values = custom_colors)) + plot_layout(guides = "collect"))
    } else {
        pdf(out_path[['hg38plot']], width=13.5, height=4.6)
        hifi.p1 = plot_prec_recall(grid, metrics_hg38) + ggtitle("COLO829 HiFi tumor-normal pair")
        ont.p2 = plot_prec_recall(grid, metrics_hg38.ont) + ggtitle("COLO829 ONT tumor-normal pair")
        print((hifi.p1 + ont.p2 + mixed_hg38_p) + plot_layout(guides = "collect") +scale_colour_manual(values = custom_colors))
    }
    dev.off()

    plot_prec_recall(grid, metrics_hg38.ont)
    ggsave(out_path[['hg38plot_ont']], width=5.5, height=5.6)

    chm13_ont_p = plot_prec_recall(grid, metrics_chm13.ont)
    ggsave(out_path[['chm13plot_ont']], width=5.5, height=5.6)

    pdf(out_path[['chm13plot']], width=10.5, height=5.6)
    print(table(metrics_chm13.ont$tool))
    p1 = plot_prec_recall(grid, metrics_chm13)
    print((p1 + chm13_ont_p) + plot_layout(guides = "collect") +scale_colour_manual(values = custom_colors))
    dev.off()

    #plot_prec_recall(grid, metrics_chm13.mixed)
    #ggsave(out_path[['mixedchm13plot']], width=5.5, height=4.6)
}

do_bar_chart(snakemake@input, snakemake@output, snakemake@threads, snakemake@params)
