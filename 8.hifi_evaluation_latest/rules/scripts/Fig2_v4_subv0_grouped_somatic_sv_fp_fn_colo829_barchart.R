library(ggplot2)
library(ggpattern)
library(stringr)
library(tidyverse)
library(ggthemes)
library(patchwork)


custom_colors <- c(
  "FP" = rgb(46, 108, 128, maxColorValue = 255),
  "TP" = rgb(233, 127, 90, maxColorValue = 255)
)

custom_colors2 <- c(
  "FP" = "gray",
  "TP" = "gray"
)


get_theme <- function(size=12, angle=0) {
    defined_theme = theme_minimal(base_size=size) + theme(strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05), legend.position="bottom", legend.box = "horizontal", legend.title=element_blank())
    defined_theme
}


load_data <- function(paths, out_path) {
    df_list = list()
    for (path in unlist(paths)) {
        res = read_tsv(path, col_names = FALSE)
        tools = res[, ncol(res)] %>% pull()
        print(dim(res))
         
        if (grepl("ont", path, fixed=T)) {
	    res = as.data.frame(res)[, c(-1, -13)]
        } else {
	    res = as.data.frame(res)[, c(-1, -15)]
        }

        # sensitivity
        sensitivity = res %>% select(X2) %>% pull()
	tp = sensitivity[2:length(sensitivity)] * sensitivity[1]
	fn = round(sensitivity[1] - tp)

        ## precision
	total_pred = as.numeric(diag(as.matrix(res)))
	total_pred = total_pred[2:length(total_pred)]
        specificity = as.numeric(unlist(as.data.frame(res)[1, 2:ncol(res)]))
	tp = round(total_pred * specificity)
	fp = total_pred - tp
        count = str_extract(path, "count(\\d+)", group=1)
	metrics = as_tibble(data.frame(FP=fp, TP=tp, tools=tools[-1], count=count))
	metrics = metrics %>% mutate(cell_line=gsub("_[hifi1|ont1]_grch38_count\\d_eval.tsv", "", basename(path)))
        df_list[[path]] = metrics
    }

    # R code
    metrics = as_tibble(do.call(rbind, df_list))
    metrics = metrics %>% pivot_longer(cols=c("FP", "TP"), names_to="metrics", values_to = "SV_num")
    metrics$count = as.numeric(metrics$count)

    if (grepl("ont", path, fixed=T)) {
        metrics = metrics %>% mutate(tool = ifelse(
          grepl("severus", tools),
          "severus", 
          ifelse(grepl("snf", tools),
                 "sniffles2",
                 ifelse(grepl("msv_ltgs", tools),
                 "minisv_ltg",
                 ifelse(grepl("hg38l_l\\+tg", tools),
                 "minisv_ltg",
                 ifelse(grepl("msv_lts_", tools),
                 "minisv_lts",
                 ifelse(grepl("savana", tools),
                 "savana",
                 ifelse(grepl("nanomonsv", tools),
                 "nanomonsv",
                 NA
                 ))
                 )
                 )
          ))
        ))
    } else {
        metrics = metrics %>% mutate(tool = ifelse(
          grepl("severus", tools),
          "severus", 
          ifelse(grepl("snf", tools),
                 "sniffles2",
                 ifelse(grepl("msv_ltgs", tools),
                 "minisv_ltg",
                 ifelse(grepl("l\\+tg\\.pair", tools),
                 "minisv_ltg",
                 ifelse(grepl("msv_lts_", tools),
                 "minisv_lts",
                 ifelse(grepl("savana", tools),
                 "savana",
                 ifelse(grepl("nanomonsv", tools),
                 "nanomonsv",
                 ifelse(grepl("hg38l\\+t", tools),
                 "msv_lt",
                 NA
                 ))
                 )
                 )
          ))
        ))
     )
    }
    metrics = metrics %>% mutate(group=paste0(ifelse(grepl("asm\\.vcf", tools), "asm_keep", "caller")))
    metrics = metrics %>% filter(tool != "msv_lt")
    metrics = metrics %>% filter(tool != "minisv_lts")
    metrics = metrics %>% mutate(
       tool=case_when(
          tool=="minisv_ltg" ~ "minisv",
          .default=tool
       )
    )

    print(head(metrics))
    metrics = metrics[,-1] %>% pivot_wider(names_from=c("group"), values_from='SV_num')

    metrics = metrics %>% mutate(asm_filter=caller-asm_keep)

    metrics = metrics %>% select(count,cell_line,metrics,tool,asm_keep,asm_filter) %>% pivot_longer(cols=c("asm_filter", "asm_keep"), names_to="group", values_to="SV_num")

    metrics = metrics %>% mutate(
        group=case_when(
            group == "asm_keep" ~ "kept by asm",
            group == "asm_filter" ~ "filtered by asm",
            .default = group,
        )
    )
    write_tsv(metrics, out_path)
}


do_bar_chart <- function(data_path, out_path, threads, myparam) {

    metrics = load_data(data_path[['stat']], out_path[['stat']])
    metrics_ont = load_data(data_path[['stat_ont']], out_path[['stat_ont']])

    metrics = metrics %>% mutate(tool=
       case_when(
           tool=="severus" ~ "Severus",
           tool=="savana" ~ "SAVANA",
           tool=="sniffles2" ~ "Sniffles2",
           tool=="delly" ~ "Delly",
           tool=="cutesv" ~ "cuteSV",
           tool=="svision" ~ "SVision-pro",
           .default = tool
       )
    )

    metrics_ont = metrics_ont %>% mutate(tool=
       case_when(
           tool=="severus" ~ "Severus",
           tool=="savana" ~ "SAVANA",
           tool=="sniffles2" ~ "Sniffles2",
           tool=="delly" ~ "Delly",
           tool=="cutesv" ~ "cuteSV",
           tool=="svision" ~ "SVision-pro",
           .default = tool
       )
    )

    pdf(out_path[['bar_pdf_2_ont_bw']], width=6, height=2.8)
    p3 = ggplot(data=metrics_ont %>% filter(metrics == 'FP')) +
      geom_bar_pattern(
          aes(
	        x=count, y=SV_num,
                fill=factor(metrics), 
		pattern=factor(group),
	  ), 
          stat = "identity", position = 'stack',
          colour = 'black',
	  pattern_fill = "black",
	  pattern_angle = 45,
	  pattern_density = 0.03,
	  pattern_key_scale_factor = 0.6,
          pattern_spacing = 0.05) + 
      scale_fill_manual(values = custom_colors2) + ylim(0, 550) + 
      scale_pattern_manual(values = c(`kept by asm` = "none", `filtered by asm` = "stripe")) +
      xlab("") +
      get_theme(angle=0, size=9) +
      facet_wrap(~tool, ncol=5) + 
      ggtitle("") + ylab("#FP SVs") + theme(legend.position='none')
      #guides(fill = "none", 
      #       pattern = guide_legend(override.aes = list(
      #         pattern = c("stripe", "none"),
      #         fill = c(`kept by asm` = "gray", `filtered by asm` = "gray"))))
    print(p3)
    dev.off()

    pdf(out_path[['bar_pdf_1_ont_bw']], width=6, height=2.8)
    p4 = ggplot(data=metrics_ont %>% filter(metrics == 'TP')) +
      geom_bar_pattern(
          aes(
	        x=count, y=SV_num,
                fill=factor(metrics), 
		pattern=factor(group),
	  ), 
          stat = "identity", position = 'stack',
          colour = 'black',
	  pattern_fill = "black",
	  pattern_angle = 45,
	  pattern_density = 0.03,
	  pattern_key_scale_factor = 0.6,
          pattern_spacing = 0.05) + 
      scale_fill_manual(values = custom_colors2) + 
      scale_pattern_manual(values = c(`kept by asm` = "none", `filtered by asm` = "stripe")) +
      xlab("") +
      get_theme(angle=0, size=9) +
      facet_wrap(~tool, ncol=5) + 
      ggtitle("") + ylab("#TP SVs") + ylim(0, 65) + 
      get_theme(angle=0, size=9) + geom_hline(yintercept=58, color='black')+
      guides(fill = "none", 
             pattern = guide_legend(override.aes = list(
               pattern = c("stripe", "none"),
               fill = c(`kept by asm` = "gray", `filtered by asm` = "gray"))))
    print(p4)
    dev.off()

    pdf(out_path[['bar_pdf_1']], width=6, height=2.8)
    p1 = ggplot(data=metrics %>% filter(metrics == 'FP')) +
      geom_bar_pattern(
          aes(
	        x=count, y=SV_num,
                fill=factor(metrics), 
		pattern=factor(group),
	  ), 
          stat = "identity", position = 'stack',
          colour = 'black',
	  pattern_fill = "black",
	  pattern_angle = 45,
	  pattern_density = 0.03,
	  pattern_key_scale_factor = 0.6,
          pattern_spacing = 0.05) + ylim(0, 550) +
      scale_fill_manual(values = custom_colors) + 
      scale_pattern_manual(values = c(`kept by asm` = "none", `filtered by asm` = "stripe")) +
      xlab("") +
      get_theme(angle=0, size=9) +
      facet_wrap(~tool, ncol=5) + 
      ggtitle("") + ylab("#FP SVs") + theme(legend.position='none')
      guides(fill = "none", 
             pattern = guide_legend(override.aes = list(
               pattern = c("stripe", "none"),
               fill = c(`kept by asm` = rgb(46, 108, 128, maxColorValue = 255), `filtered by asm` = rgb(46, 108, 128, maxColorValue = 255)))))
    print(p1)
    dev.off()

    pdf(out_path[['bar_pdf_1_bw']], width=6, height=2.8)
    p1 = ggplot(data=metrics %>% filter(metrics == 'FP')) +
      geom_bar_pattern(
          aes(
	        x=count, y=SV_num,
                fill=factor(metrics), 
		pattern=factor(group),
	  ), 
          stat = "identity", position = 'stack',
          colour = 'black',
	  pattern_fill = "black",
	  pattern_angle = 45,
	  pattern_density = 0.03,
	  pattern_key_scale_factor = 0.6,
          pattern_spacing = 0.05) + ylim(0, 550) +
      scale_fill_manual(values = custom_colors2) + 
      scale_pattern_manual(values = c(`kept by asm` = "none", `filtered by asm` = "stripe")) +
      xlab("") +
      get_theme(angle=0, size=9) +
      facet_wrap(~tool, ncol=5) + 
      ggtitle("") + ylab("#FP SVs") + theme(legend.position='none') +
      guides(fill = "none", 
             pattern = guide_legend(override.aes = list(
               pattern = c("stripe", "none"),
               fill = c(`kept by asm` = "gray", `filtered by asm` = "gray"))))
    print(p1)
    dev.off()

    pdf(out_path[['bar_pdf_2']], width=6, height=2.8)
    p0 = ggplot(data=metrics %>% filter(metrics == 'TP')) +
      geom_bar_pattern(
          aes(
	        x=count, y=SV_num,
                fill=factor(metrics), 
		pattern=factor(group),
	  ), 
          stat = "identity", position = 'stack',
          colour = 'black',
	  pattern_fill = "black",
	  pattern_angle = 45,
	  pattern_density = 0.03,
	  pattern_key_scale_factor = 0.6,
          pattern_spacing = 0.05) + 
      scale_fill_manual(values = custom_colors) + 
      scale_pattern_manual(values = c(`kept by asm` = "none", `filtered by asm` = "stripe")) +
      xlab("") +
      facet_wrap(~tool, ncol=5) + 
      get_theme(angle=0, size=9) + geom_hline(yintercept=58, color='red')+
      ggtitle("") + ylab("#TP SVs")+
      guides(fill = "none", 
             pattern = guide_legend(override.aes = list(
               pattern = c("stripe", "none"),
               fill = c(`kept by asm` = rgb(233, 127, 90, maxColorValue = 255), `filtered by asm` = rgb(233, 127, 90, maxColorValue = 255)))))
    print(p0)
    dev.off()

    pdf(out_path[['bar_pdf_2_bw']], width=6, height=2.8)
    p2 = ggplot(data=metrics %>% filter(metrics == 'TP')) +
      geom_bar_pattern(
          aes(
	        x=count, y=SV_num,
                fill=factor(metrics), 
		pattern=factor(group),
	  ), 
          stat = "identity", position = 'stack',
          colour = 'black',
	  #pattern_fill = "black",
	  pattern_angle = 45,
	  pattern_density = 0.03,
	  pattern_key_scale_factor = 0.6,
          pattern_spacing = 0.05) +  ylim(0, 65) + 
      scale_fill_manual(values = custom_colors2) + 
      scale_pattern_manual(values = c(`kept by asm` = "none", `filtered by asm` = "stripe")) +
      xlab("") +
      facet_wrap(~tool, ncol=5) + 
      get_theme(angle=0, size=9) + geom_hline(yintercept=58, color='black')+
      ggtitle("") + ylab("#TP SVs")+
      guides(fill = "none", 
             pattern = guide_legend(override.aes = list(
               pattern = c("stripe", "none"),
               fill = c(`kept by asm` = "gray", `filtered by asm` = "gray"))))
    print(p2)
    dev.off()

    pdf(out_path[['bar_pdf_all_bw']], width=9, height=6)
    print(((p1+p3) / (p2+p4))+plot_annotation(tag_levels = list(c('A', 'B', 'C', 'D'), '1'))) # + plot_layout(guides='collect', ncol=2, nrow=2, widths = c(4, 4), heights=c(1, 1)) & theme(legend.position = "bottom"))
    dev.off()
}


do_bar_chart(snakemake@input, snakemake@output, snakemake@threads, snakemake@params)


