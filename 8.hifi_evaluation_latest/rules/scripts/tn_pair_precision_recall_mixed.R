library(ggplot2)
library(ggthemes)
library(stringr)
library(tidyverse)
library(patchwork)

custom_colors <- c(
  "msv:tg" = rgb(249, 134, 130, maxColorValue = 255),
  "minisvl+tg" = rgb(249, 134, 130, maxColorValue = 255),
  "minisvl+tgs" = rgb(187, 142, 33, maxColorValue = 255),
  "msv:tgsnonormal" = rgb(22, 183, 139, maxColorValue = 255),
  "msv:tgs" = rgb(187, 142, 33, maxColorValue = 255),
  "minisvl+gs" = rgb(187, 142, 33, maxColorValue = 255),
  "minisvl+g" = rgb(249, 134, 130, maxColorValue = 255),
  "nanomonsv" = rgb(147, 203, 118, maxColorValue = 255),
  "savana" = rgb(22, 183, 139, maxColorValue = 255),
  "severus_mosaic" = rgb(16, 174, 228, maxColorValue = 255),
  "severus_somatic" = "orange",
  "severus_lowaf" = rgb(16, 174, 228, maxColorValue = 255),
  "sniffles" = rgb(154, 130, 251, maxColorValue = 255),
  "sniffles_somatic" = rgb(248, 89, 206, maxColorValue = 255),
  "snf_mosaic" = rgb(154, 130, 251, maxColorValue = 255),
  "svision" = rgb(248, 89, 206, maxColorValue = 255),
  "severus_mosaic_filterasm" = "black",
  "severus_mosaic_rawasm" = "gray"
)

get_theme <- function(size=12, angle=0) {
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
    metrics$tool = 
	    ifelse(grepl("somatic_generation", metrics$file) & (!grepl("severus", metrics$file)),
		   "sniffles_somatic",
	    ifelse(grepl('grch38l', metrics$file), 
			  str_extract(metrics$file, "(gafcall|minisv_pair|minisv_mosaic|nanomonsv|savana|severus_lowaf|sniffles_mosaic|sniffles|svision|cutesv)", group=1),
                          str_extract(metrics$file, "(severus_lowaf|sniffles_mosaic|sniffles|svision|cutesv|severus_all|severus_somatic|severus_af)", group=1)))

    metrics$tool[grepl('rawasm', metrics$file)] = 'severus_mosaic_rawasm'
    metrics$tool[grepl('filterasm', metrics$file)] = 'severus_mosaic_filterasm'

#    metrics$tool = ifelse(grepl("grch38g|grch38l|hg38l|chm13l|chm13g", metrics$file),  # minisv
#			  gsub("_pair", "", gsub("gafcall", "minisv", str_extract(metrics$tool, "(gafcall|minisv|minisv_pair|minisv_mosaic)", group=1))),
#			  gsub("_pair", "", gsub("gafcall", "minisv", str_extract(metrics$tool, "(nanomonsv|savana|severus|severus_lowaf|sniffles_mosaic|sniffles|svision)", group=1))))
#
    metrics$tool = gsub('gafcall', 'minisv', metrics$tool)
    metrics$tool = gsub('sniffles_mosaic', 'snf_mosaic', metrics$tool)
    metrics$tool = gsub('severus_all', 'severus_mosaic', metrics$tool)
    metrics$tool = gsub('severus_af', 'severus_mosaic', metrics$tool)

    print(table(metrics$tool))

    # previous gafcall output path name
    # different from current minisv 
    # for ultralong reads COLO829
    metrics$param = str_extract(metrics$file, "(l\\+g|l_l\\+t\\+g|l_l\\+t\\+g\\+s|l_l\\+x|g_g\\+x|l_l\\+g|l_l\\+t|l\\+tgs|l\\+gs|l\\+ts|l\\+tg)_mixed", group=1)
    metrics$param[is.na(metrics$param)] = ""

    metrics$normal = str_extract(metrics$file, "(nonormal)", group=1)
    metrics$normal[is.na(metrics$normal)] = ""
    metrics$tool = paste0(metrics$tool, metrics$param, metrics$normal)

    metrics$tool = gsub('minisv_mosaicl_l\\+t\\+g\\+s', 'msv:tgs', metrics$tool)
    metrics$tool = gsub('minisv_mosaicl_l\\+t\\+g', 'msv:tg', metrics$tool)
    print('--------')
    print(table(metrics$tool))
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
        geom_contour(data=grid, aes(x = Recall, y = Precision, z = F1), linetype="dashed", linewidth=0.45, color='gray', bins = 10) + 
        geom_point(data=metrics, aes(x=sensitivity, y=specificity, colour = factor(tool), size=count), alpha=0.5) + ylab('Precision') + xlab("Recall") + scale_size_continuous(name = "Count", breaks = c(2, 3, 4, 5), range = c(1, 5)) + 
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
    data_path_mixed = input[['mixed_colo829_hifi']]
    data_path_mixed_100kb = input[['mixed_colo829_hifi_100kb']]
    print(data_path_mixed_100kb)

    metrics.mixed = load_metrics(data_path_mixed)
    metrics_hg38.mixed  = metrics.mixed %>% filter(genome == 'hg38' & count >= 2 & count <= 5)
    write_tsv(metrics.mixed, out_path[['mixed_table']])

    plot_100kb = F
    if (length(unlist(data_path_mixed_100kb)) > 0) {
        plot_100kb = T
        metrics.mixed_100kb = load_metrics(data_path_mixed_100kb)
        metrics.hg38_mixed_100kb = metrics.mixed_100kb %>% filter(genome == 'hg38' & count >= 2 & count <= 5)
        write_tsv(metrics.hg38_mixed_100kb, out_path[['mixed_table_100kb']])
    }

    #####metrics_hg38.mixed$tool = gsub("minisvl\\+t\\+g\\+s", "minisvl\\+tgs", gsub("_mosaic", "", gsub("_mosaic ", "", metrics_hg38.mixed$tool)))
    metrics_hg38.mixed$tool = gsub("minisvl\\+t\\+g", "minisvl\\+tg", metrics_hg38.mixed$tool)
    print(metrics_hg38.mixed)
    
    if (plot_100kb) {
        ###metrics.hg38_mixed_100kb$tool = gsub("minisvl\\+t\\+g\\+s", "msv:tgs", gsub("_mosaic", "", metrics.hg38_mixed_100kb$tool))
        metrics.hg38_mixed_100kb$tool = gsub("minisvl\\+t\\+g", "msv:tg", metrics.hg38_mixed_100kb$tool)
    }

    grid <- generate_grid()
    pdf(out_path[['mixedhg38plot']], width=5.5, height=5.6)
    mixed_hg38_p = plot_prec_recall(grid, metrics_hg38.mixed) + ggtitle("COLO829 HiFi tumor-normal\nmixed reads all SVs")+scale_colour_manual(values = custom_colors)
    print(mixed_hg38_p)
    dev.off()

    # Plot contours of constant F1 score
    if (plot_100kb) {
        pdf(out_path[['mixedhg38plot']], width=9.5, height=4)
        mixed_hg38_p_100kb = plot_prec_recall(grid, metrics.hg38_mixed_100kb)+scale_colour_manual(values = custom_colors) + ggtitle("COLO829 HiFi tumor-normal\nmixed reads >100kb SVs")
        print((mixed_hg38_p + mixed_hg38_p_100kb) + plot_layout(guides = "collect"))
    }
    dev.off()

    #plot_prec_recall(grid, metrics_chm13.mixed)
    #ggsave(out_path[['mixedchm13plot']], width=5.5, height=4.6)
}

do_bar_chart(snakemake@input, snakemake@output, snakemake@threads, snakemake@params)
