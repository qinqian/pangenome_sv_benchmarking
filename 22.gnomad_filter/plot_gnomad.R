library(tidyverse)
library(ggplot2)
library(ggpattern)
library(ggthemes)
stat_dir <- "."   # change to your folder path if needed

df <- list.files(stat_dir, pattern = "*\\.stat$", full.names = TRUE) %>%
  map_dfr(function(f) {
    bn <- basename(f)

    # Example filename:
    # HCC1954_hifi1_hifi1_savana13_gnomadaf001.stat
    meta <- str_match(
      bn,
      "^(.*?)_(hifi1|ont1|ont2)_hifi1_(savana13|severus16|nanomonsv|sniffles2)_(gnomad(?:af)?[0-9]+)\\.stat$"
    )

    sample <- meta[, 2]
    platform <- meta[,3]
    tool <- meta[, 4]
    cutoff_tag <- meta[, 5]

    x <- readLines(f, warn = FALSE)
    x <- x[nzchar(trimws(x))]

    # Treat empty files as missing rather than zero
    if (length(x) == 0) {
      return(tibble(
        sample = sample,
        tool = tool,
        cutoff = cutoff_tag,
        total = NA_integer_,
        filtered = NA_integer_
      ))
    }

    nums <- scan(text = x[1], quiet = TRUE)

    tibble(
      sample = sample,
      cutoff = cutoff_tag,
      tool = tool,
      platform=platform,
      total = nums[1],
      filtered = nums[2]
    )
  }) %>%
  mutate(
    cutoff = recode(
      cutoff,
      gnomad01 = "0.01",
      gnomadaf001 = "0.001",
      gnomadaf0 = "0"
    )
  )

df = df %>% filter(!platform%in%c('ont1', 'ont2'))
df = df %>% filter(cutoff != "0")
# Check combined table

df =  df %>% mutate(tool=case_when(
tool == 'savana13' ~ 'SAVANA',
tool == 'severus16' ~ 'Severus',
tool == 'sniffles2' ~ 'Sniffles2',
.default = tool
))
write_tsv(df, file='write_summary.tsv')

plot_df <- df %>%
  filter(!is.na(total), !is.na(filtered)) %>%
  mutate(kept=total-filtered) %>%  select(-total) %>% 
  pivot_longer(
    cols = c(kept, filtered),
    names_to = "status",
    values_to = "count"
  ) %>% 
  mutate(
    status = recode(status, kept = "Kept SVs", filtered = "Filtered SVs"),
    status = factor(status, levels = c("Filtered SVs", "Kept SVs")),
    cutoff = factor(cutoff)
  )

write_tsv(plot_df, file='write_summary2.tsv')

p <- ggplot(plot_df, aes(x = tool, y = count, fill = status)) +
  geom_col(width = 0.75) +
  facet_wrap(sample~ cutoff, nrow = 6, scales='free') +
  labs(
    x = "SV Tools", 
    y = "Number of SVs",
    fill = NULL,
    title = "SVs under different gnomAD AF cutoffs"
  ) +
  theme_clean(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )
print(p)

ggsave(
  "gnomad_filter_stacked_barplot.pdf",
  p,
  width = 8,
  height = 13
)

