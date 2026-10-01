library(ggplot2)
library(ggpattern)
library(stringr)
library(tidyverse)
library(ggthemes)
library(patchwork)

custom_colors_fnfp <- c(
  "FP" = rgb(46, 108, 128, maxColorValue = 255),
  "FN" = rgb(233, 127, 90, maxColorValue = 255)
)

custom_colors_fnfp_bw <- c(
  "FP" = "gray70",
  "FN" = "gray70"
)

get_theme <- function(size = 12, angle = 0) {
  defined_theme <- theme_minimal(base_size = size) +
    theme(
      strip.text = element_text(size = size),
      legend.text = element_text(size = size),
      axis.title.x = element_text(size = size),
      axis.title.y = element_text(size = size),
      axis.text.y = element_text(size = size),
      axis.text.x = element_text(
        size = size,
        angle = angle,
        hjust = 1,
        vjust = 1.05
      ),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_blank()
    )

  defined_theme
}


load_minisv_gnomad_fnfp <- function(stat_path, out_path = NULL) {
  res <- read_tsv(
    stat_path,
    col_names = FALSE,
    show_col_types = FALSE
  )
  callset_paths <- res[[ncol(res)]]

  mat <- res[, 2:(ncol(res) - 1)] %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix()

  rownames(mat) <- callset_paths
  colnames(mat) <- callset_paths

  sv_counts <- diag(mat)

  # Because you reordered eval output, truth should be row 1 / column 1.
  # Still keep this check to make the script safer.
  truth_idx <- 1

  if (!str_detect(callset_paths[truth_idx], "truth")) {
    stop("The first row is not recognized as truth. Please check input order.")
  }

  truth_count <- sv_counts[truth_idx]
  print(truth_count)

  metrics <- tibble(
    path = callset_paths,
    row_idx = seq_along(callset_paths),
    SV_total = sv_counts,

    # This is precision-like and is used to estimate FP.
    frac_caller_in_truth = mat[truth_idx, ],

    # This is sensitivity-like and is used to estimate FN.
    frac_truth_in_caller = mat[, truth_idx]
  ) %>%
    mutate(
      TP_for_FP = round(SV_total * frac_caller_in_truth),
      FP = SV_total - TP_for_FP,

      TP = round(truth_count * frac_truth_in_caller),
      #FN = truth_count - TP_for_FN,
      tool = case_when(
        str_detect(path, "nanomonsv") ~ "nanomonsv",
        str_detect(path, "savana13|savana") ~ "SAVANA",
        str_detect(path, "severus16|severus") ~ "Severus",
        str_detect(path, "sniffles2|snf") ~ "Sniffles2",
        TRUE ~ NA_character_
      ),
      gnomad_filter = case_when(
        str_detect(path, "filtered\\.gnomad001\\.vcf") ~ "gnomad001",
        str_detect(path, "filtered\\.gnomad01\\.vcf") ~ "gnomad01",
        #str_detect(path, "filtered\\.gnomad0\\.vcf") ~ "gnomad0",
        TRUE ~ "original"
      )
    ) 
  metrics = metrics %>%
    filter(row_idx != 1) %>%
    filter(!is.na(tool)) %>%
    filter(tool != "Sniffles2") %>%
    filter(!str_detect(path, "truth")) %>%
    mutate(
	   truth_count=truth_count,
      tool = factor(
        tool,
        #levels = c("nanomonsv", "SAVANA", "Severus", "Sniffles2")
        levels = c("nanomonsv", "SAVANA", "Severus")
      ) # ,
      #gnomad_filter = factor(
      #  gnomad_filter,
      #  levels = c("original", "gnomad01", "gnomad001", "gnomad0")
      #)
    ) %>%
    select(
      tool,
      gnomad_filter,
      path,
      SV_total,
      FP,
      #FN,
      TP,
      truth_count,
      #TP_for_FP, TP_for_FN,
      frac_caller_in_truth,
      frac_truth_in_caller
    )


 write_tsv(metrics, file='metrics_test.tsv')

#tool    TP_for_FP       TP_for_FN       gnomad_filter   path    SV_total        FP      FN      frac_caller_in_truth    frac_truth_in_caller
#nanomonsv       46      46      gnomad001       COLO829_hifi1_hifi1_nanomonsv_filtered.gnomad001.vcf    52      6       12      0.8846  0.7931
#nanomonsv       46      46      gnomad01        COLO829_hifi1_hifi1_nanomonsv_filtered.gnomad01.vcf     54      8       12      0.8519  0.7931

  metrics_long1 <- metrics %>%
    #select(tool, gnomad_filter, FP, TP_for_FP) %>%
	 mutate(SV_num=FP) %>% 
    select(tool, gnomad_filter, SV_num) %>%
    mutate(metrics="FP") %>%
    #pivot_longer(
    #  cols = c("FP"),
    #  names_to = "metrics",
    #  values_to = "SV_num"
    #) %>%
    mutate(
      #metrics = factor(metrics, levels = c("FP", "TP_for_FP")),
      #gnomad_filter = factor(gnomad_filter, levels=c("original", "gnomad01", "gnomad001", "gnomad0"))
      gnomad_filter = factor(gnomad_filter, levels=c("original", "gnomad01", "gnomad001"))
    )

  metrics_long2 <- metrics %>%
    #select(tool, gnomad_filter, FN, TP_for_FN) %>%
	 mutate(SV_num=TP) %>% 
    select(tool, gnomad_filter, SV_num) %>%
    #pivot_longer(
    #  cols = c("FN", "TP_for_FN"),
    #  names_to = "metrics",
    #  values_to = "SV_num"
    #) %>%
    mutate(metrics="TP") %>%
    mutate(
      #metrics = factor(metrics, levels = c("FN", "TP_for_FN")),
      #gnomad_filter = factor(gnomad_filter, levels=c("original", "gnomad01", "gnomad001", "gnomad0"))
      gnomad_filter = factor(gnomad_filter, levels=c("original", "gnomad01", "gnomad001"))
    )

  metrics_long = bind_rows(metrics_long1, metrics_long2)

  if (!is.null(out_path)) {
    write_tsv(metrics_long, out_path)
  }

  return(metrics_long)
}


plot_minisv_gnomad_fnfp <- function(metrics, bw = FALSE) {
  fill_colors <- if (bw) {
    custom_colors_fnfp_bw
  } else {
    custom_colors_fnfp
  }
  p <- ggplot(data = metrics) +
    #geom_bar_pattern(
    geom_bar(
      aes(
        x = gnomad_filter,
        y = SV_num,
        #fill = metrics,
        #pattern = metrics
      ),
      stat = "identity",
      position = "stack",
      colour = "black",
      #pattern_fill = "black",
      #pattern_angle = 45,
      #pattern_density = 0.03,
      #pattern_key_scale_factor = 0.6,
      #pattern_spacing = 0.05
    ) +
    #geom_text(
    #  aes(
    #    x = gnomad_filter,
    #    y = SV_num,
    #    label = ifelse(SV_num > 0, SV_num, "")
    #  ),
    #  position = position_stack(vjust = 0.5),
    #  size = 3
    #) +
    #scale_fill_manual(values = fill_colors) +
    #scale_pattern_manual(values = c(
    #  "FP" = "none",
    #  "FN" = "stripe"
    #)) +
    facet_wrap(~tool, ncol = 4) +
    xlab("") +
    ylab("#SVs") +
    ggtitle("False positives and false negatives after gnomAD filtering") +
    get_theme(angle = 35, size = 10) # +
    #guides(
    #  fill = guide_legend(
    #    override.aes = list(
    #      pattern = c("none", "stripe")
    #    )
    #  ),
    #  pattern = "none"
    #)
  return(p)
}

plot_fn_fp <- function(x, title="COLO829 HiFi") {
  metrics <- load_minisv_gnomad_fnfp(
    stat_path = x,
    out_path = "test.csv"
  )
  p1 <- plot_minisv_gnomad_fnfp(
    metrics %>% filter(metrics=='FP'),
      #filter(as.character(.data$metrics) %in% c("TP_for_FP", "FP")),
      #select(c("TP_for_FP", "FP")),
    bw = FALSE
  ) +
    ylab("#FP SVs") +
    ggtitle(title)
  p2 <- plot_minisv_gnomad_fnfp(
    #metrics = metrics %>%
    #  filter(as.character(.data$metrics) %in% c("TP_for_FN", "FN")),
    metrics %>% filter(metrics=='TP'),
    bw = FALSE
  ) +
    ylab("#TP SVs") + ylim(0, 60) + geom_hline(yintercept=58, color='black') +
    ggtitle("")
  p1 / p2
}

do_gnomad_fnfp_bar_chart <- function(out_path, threads = 1, myparam = NULL) {
  data_path <- "COLO829_hifi1_eval.tsv"
  pp1 = plot_fn_fp(data_path)
  data_path <- "COLO829_ont1_eval.tsv"
  pp2 = plot_fn_fp(data_path, title="COLO829 ONT")

  pdf(out_path[["bar_pdf"]], width = 8, height = 6)
  print(pp1 | pp2)
  dev.off()
}

do_gnomad_fnfp_bar_chart(list(tsv='gnomad_outfig.tsv', bar_pdf='gnomad_outfig.pdf'))

