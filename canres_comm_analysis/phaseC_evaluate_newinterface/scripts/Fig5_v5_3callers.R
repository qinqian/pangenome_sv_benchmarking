library(ggplot2)
library(ggpattern)
library(stringr)
library(cowplot)
library(tidyverse)
library(ggthemes)
library(patchwork)


custom_colors <- c(
  "filtered by asm" = rgb(46, 108, 128, maxColorValue = 255),
  "kept by asm" = rgb(233, 127, 90, maxColorValue = 255)
)

custom_colors2 <- c(
  "kept by filters" = rgb(127, 0, 202, max = 255),
  "filtered by lg" = rgb(20, 143, 97, max = 255),
  "filtered by lgs" = rgb(220, 141, 13, max = 255)
)

get_theme <- function(size=12, angle=0) {
    #defined_theme = theme_clean(base_size=size) + theme(legend.title=element_text(size=size), strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05)) #, legend.position="bottom", legend.box = "horizontal") 
    defined_theme = theme_clean(base_size=size) + theme(legend.title=element_text(size=size), strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05), legend.position="none")
    defined_theme
}



upset_simple = function(df, cl, xlog10=F, ylog10=F, axis_off=F, label="", theme=F) {
  df = df %>% select(-l_count, -cell_line, -lg_filter, -lgs_filter)
  print('---------')
  print(head(df))
  # total count sort
  df = df %>% arrange(desc(count+count2+count3))
  perf_metrics = data.frame()
  ##perf_metrics = data.frame(
  ##  x=1:n,
  ##  y=sapply(1:n, function(i) sum(df[which(df[,i] == 1), "count"])),
  ##  class='asm_keep'
  ##)
  ##perf_metrics = perf_metrics %>% mutate(tp=y, fn=total_asm-y) %>% mutate(tool=case_when(
  ##    x==1 ~ "Severus",
  ##    x==2 ~ "SAVANA",
  ##    x==3 ~ "nanomonsv",
  ##))
  n = 3 # sv caller number
  print(n)
  sums = data.frame(
    x=1:n,
    y=sapply(1:n, function(i) sum(df[which(df[,i] == 1), "count"])),
    class='asm_keep'
  )
  sums2 = data.frame(
    x=1:n,
    y=sapply(1:n, function(i) sum(df[which(df[,i] == 1), "count2"])),
    class='filter'
  )
  sums3 = data.frame(
    x=1:n,
    y=sapply(1:n, function(i) sum(df[which(df[,i] == 1), "count3"])),
    class='filter2'
  )
  print('---------')
  sums = rbind(sums, sums2, sums3)
  sums = sums %>% mutate(
        class=case_when(
            class == "asm_keep" ~ "kept by filters",
            class == "filter" ~ "filtered by lgs",
            class == "filter2" ~ "filtered by lg",
            .default = class,
        )
  )
  print(sums)

  if (xlog10) {
    sums$y = log10(sums$y)
    ylabel = "Set Size (log 10)"
  } else {
    ylabel = "#SV"
  }
  
  if (ylog10) {
    df$count = log10(df$count)
    xlabel = "Intersection Size (log10)"
  } else {
    xlabel = "#SV"
  }

  if (axis_off) {
      axis_x = element_blank()
      labels = rep("", n)
  } else {
      axis_x = element_text(vjust=1)
      labels = colnames(df)[1:n]
  }
  if (theme) {
     legend="left"
  } else {
     legend="none"
  }
  p1 = ggplot(sums, aes(x=x, y=y, fill=factor(class))) + 
    geom_bar(stat="identity", width=0.3) +
    scale_x_continuous(breaks=1:n, labels=labels, limits=c(0.5, n+0.5), position = "bottom") + 
    theme_minimal() +
    theme(
      text = element_text(size = 20),
      axis.title.y=element_blank(),
      axis.title.x = element_text(vjust=1, size=12),
      axis.text.x  = element_text(vjust=1, size=12),
      panel.grid.major.y=element_blank(),
      panel.grid.minor.y=element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "null"),
      legend.position="none"
    ) +
    scale_fill_manual(values = custom_colors2) + 
    ylab(ylabel) +
    coord_flip() + 
    scale_y_continuous(trans = "reverse")
  
  point_size=3
  line_size=0.75
  text_scale=1
  processed_df = as_tibble(df) %>% mutate(x=1:nrow(df)) %>% pivot_longer(cols=c('count', 'count2', 'count3')) %>% mutate(
        name=case_when(
            name == "count" ~ "kept by filters",
            name == "count2" ~ "filtered by lgs",
            name == "count3" ~ "filtered by lg",
            .default = name,
        )
  )

  p2 = ggplot(processed_df) + 
    geom_bar(aes(x=x, y=value, fill=factor(name)), stat="identity", width=0.5, position='stack') +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(vjust=1, size=12),
      axis.title.y = element_text(vjust=1, size=12),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(size=20, face='bold'),
      plot.margin = unit(c(0, 0, 0, 0), "null"),
      legend.text = element_text(size=20),
      legend.position=legend, 
      legend.title=element_blank(),
      plot.tag = element_text(
                  size = rel(2),
                  hjust = 0.5, vjust = 0.5)) +
    labs(title=NULL, x=NULL, y=NULL, tag=NULL) +
    ylab(xlabel) +
    scale_fill_manual(values = custom_colors2) + 
    scale_x_continuous(limits=c(0, nrow(df)+1), expand=c(0, 0)) + 
    scale_y_continuous(expand=c(0, 0)) + ggtitle(paste0(cl))

  nsets = 3
  nbars = nrow(df)
  Mat_data = make_mat_data(df)
  shading_data = make_shading_data(nsets, nbars)
  all_points = expand.grid(y=1:nsets, x=1:nbars)
  p3 = ggplot() + theme_minimal() + theme(
      plot.margin = unit(c(0, 0, 0, 0), "null"),
      legend.position="none",
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
    ) + xlab("") + ylab("") + scale_x_continuous(limits=c(0, nbars+1), expand=c(0, 0)) + 
    geom_rect(
      data=shading_data, 
      aes(
        xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax
      ), fill=shading_data$color, alpha=0.5
    ) + geom_point(
      data=all_points, aes(x=x, y=y),
      color="darkgrey", size=point_size, shape=16, alpha=0.5
    ) + geom_point(
      data=Mat_data, aes(x=x, y=y),
      size=point_size, shape=16
    ) + geom_line(
      data=Mat_data, 
      aes(group=group, x=x, y=y),
      linewidth=line_size
    ) + scale_color_identity()
  
  guide = guide_area() 
  plot = guide + p2 + p1 + p3 +  plot_layout(ncol = 2,  widths = c(4, 10), heights=c(3.2, 0.8), guides="collect")
  ####plot = guide + p2 + p1 +  plot_layout(ncol = 2,  widths = c(4, 10), heights=c(3.2, 0.8), guides="collect")
  list(perf_metrics, plot)
}


make_mat_data = function(df) {
  n = ncol(df)-3
  df = df[,1:n]
  df2 = NULL
  for (i in 1:nrow(df)) {
    w = which(df[i,] == 1)
    name = paste(colnames(df)[w], collapse="&")
    for (j in w) {
      df2 = rbind(df2, data.frame(x=i, y=j, group=name))
    }
  }
  df2
}

make_shading_data = function(nrows, ncols) {
  data.frame(
    xmin=0.5,
    xmax=ncols+0.5,
    ymin=0.5:(nrows-0.5),
    ymax=1.5:(nrows+0.5),
    color=c("gray", "white")[(0:(nrows-1) %% 2) + 1],
    alpha=0.5
  )
}


make_matrix_plot <- function(df, point_size=3, line_size=0.75, text_scale=1, shading_data=NULL) {
  nsets = ncol(df) - 3
  nbars = nrow(df)
  Mat_data = make_mat_data(df)
  if (is.null(shading_data)) {
    shading_data = make_shading_data(nsets, nbars)
  }
  all_points = expand.grid(y=1:nsets, x=1:nbars)
  
  ggplot_gtable(ggplot_build(
    ggplot()
    + theme_minimal()
    + theme(
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "null")
    )
    + xlab("") 
    + ylab("")
    + scale_y_continuous(
        labels=colnames(df)[1:nsets],
        breaks=c(1:nsets),
        limits=c(0.5,(nsets+0.5)),
        expand=c(0, 0)
    )
    + scale_x_continuous(limits=c(0, nbars+1), expand=c(0, 0))
    + geom_rect(
      data=shading_data, 
      aes_string(
        xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax"
      ), fill=shading_data$color, alpha=0.5
    )
    + geom_point(
      data=all_points, aes_string(x="x", y="y"),
      color="darkgrey", size=point_size, shape=16, alpha=0.5
    )
    + geom_point(
      data=Mat_data, aes_string(x="x", y="y"),
      size=point_size, shape=16
    )
    + geom_line(
      data=Mat_data, 
      aes_string(group="group", x="x", y="y"),
      size=line_size
    )
    + scale_color_identity()
  ))
}


do_bar_chart <- function(data_path, out_path, threads, myparam) {
    # R code
    df_list = list()
    file_num = length(data_path[['union_count']])

    for (index in seq(file_num)) {
        union_count = read.table(data_path[['union_count']][index], colClasses = 'character', sep='\t')
        g_union_count = read.table(data_path[['g_union_count']][index], colClasses = 'character', sep='\t')
        gs_union_count = read.table(data_path[['gs_union_count']][index], colClasses = 'character', sep='\t')
	metrics = as_tibble(union_count) %>% left_join(as_tibble(g_union_count), by='V1')
	metrics = as_tibble(metrics) %>% left_join(as_tibble(gs_union_count), by='V1')
	df_list[[index]] = metrics %>% mutate(cell_line = gsub("_hifi1_new_interface", "", basename(dirname(data_path[['union_count']][index]))))
return
    }

    metrics = bind_rows(df_list)
    colnames(metrics)[1:4] = c('comb', 'caller',  'lg_keep', 'lgs_keep')
 
    metrics = data.frame(do.call(rbind, lapply(strsplit(metrics$comb, ""), as.numeric)), metrics$caller, metrics$lg_keep, metrics$lgs_keep, metrics$cell_line)
    colnames(metrics) = c("Severus", "SAVANA", "nanomonsv", "l_count", "lg_keep", "lgs_keep", "cell_line")

    #metrics = data.frame(metrics[, c(1,2,3,4,5,6,7)])
    print(head(metrics))
    metrics = metrics %>% 
        mutate(l_count=as.numeric(l_count)) %>% 
        mutate(lg_keep=as.numeric(lg_keep)) %>%
        mutate(lgs_keep=as.numeric(lgs_keep)) %>%
        mutate(lg_filter=l_count-lg_keep) %>% 
        mutate(lgs_filter=lg_keep-lgs_keep) %>% 
        mutate(count=lgs_keep) %>% 
        mutate(count2=lgs_filter) %>% 
        mutate(count3=lg_filter) %>% 
        mutate(count2=ifelse(count2>0, count2, 0)) %>% 
        mutate(count3=ifelse(count3>0, count3, 0))

    write_tsv(metrics, out_path[['stat']])
    metrics = metrics %>% dplyr::select(-lg_keep) %>% dplyr::select(-lgs_keep)

    pdf(out_path[['bar_pdf']], width=20, height=11)
    #plot_list <- list() 
    #labels = toupper(letters)
    #n = 0
    #for (cl in unique(metrics$cell_line)) {
    #    n <- n + 1
    #    metrics_cl = metrics %>% filter(cell_line == cl)
    #    print(cl)
    #    if (cl == "COLO829" || cl == "HCC1954") {
    #        if (cl == "COLO829") {
    #            results = upset_simple(metrics_cl, cl, axis_off=F, label=labels[n], theme=T)
    #            plot_list[[cl]] = results[[2]]
    #        } else {
    #            results = upset_simple(metrics_cl, cl, axis_off=F, label=labels[n], theme=F)
    #            plot_list[[cl]] = results[[2]]
    #        }
    #    } else {
    #        results = upset_simple(metrics_cl, cl, axis_off=T, label=labels[n], theme=F)
    #        plot_list[[cl]] = results[[2]]
    #    }
    #}
    #print(wrap_plots(plot_list), ncol=2, guides='collect')
    dev.off()

    pdf(out_path[['bar_pdf_bw']], width=20, height=11)
    plot_list <- list() 
    labels = toupper(letters)
    n = 0
    out.list = list()
    for (cl in unique(metrics$cell_line)) {
        n <- n + 1
        metrics_cl = metrics %>% filter(cell_line == cl)
        if (cl == "COLO829" || cl == "HCC1954") {
            if (cl == "COLO829") {
                results = upset_simple(metrics_cl, cl, axis_off=F, label=labels[n], theme=T)
                plot_list[[cl]] = results[[2]]
            } else {
                results = upset_simple(metrics_cl, cl, axis_off=F, label=labels[n], theme=F)
                plot_list[[cl]] = results[[2]]
            }
        } else {
            results = upset_simple(metrics_cl, cl, axis_off=T, label=labels[n], theme=F)
            plot_list[[cl]] = results[[2]]
        }
        #out.list[[cl]] = cbind(results[[1]], cl)
    }
    print(wrap_plots(plot_list), ncol=2, guides='collect')
    dev.off()

    #write_tsv(do.call(rbind, out.list), out_path[['metrics']])
}


do_bar_chart(snakemake@input, snakemake@output, snakemake@threads, snakemake@params)


