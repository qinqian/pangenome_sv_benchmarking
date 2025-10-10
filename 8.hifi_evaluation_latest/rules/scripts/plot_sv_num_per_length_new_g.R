library(ggplot2)
library(ggbreak)
library(tidyverse)
library(ggthemes)
library(patchwork)

custom_colors <- c(
  ">100k" = rgb(245, 106, 100, maxColorValue = 255),
  ">1M" = rgb(24, 176, 61, maxColorValue = 255),
  "translocation" = rgb(85, 148, 249, maxColorValue = 255),
  "all_except_inv" = rgb(21, 183, 187, maxColorValue = 255),
  ">20k" = rgb(221, 180, 104, maxColorValue = 255)
)


get_theme <- function(size=12, angle=45) {
    defined_theme = theme_clean(base_size=size) + theme(legend.title=element_text(size=size), strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05), legend.position="bottom", legend.box = "horizontal") 
    defined_theme
}


preprocess_data <- function(x, count_thresh=2) {
    # columns: translocation   >1M     >100k   >20k    all_except_inv  inv     file
    res = x
    res.small =res[, c(4,5,8,9)]
    # 4: >20k, 5, all except inv
    # 5 -> 100bp - 20kb
    res.small[, 2] = res.small[, 2] - res[, 1]
    # 4 -> 20k - 100k
    res.small[, 1] = res.small[, 1] - res[, 3]
    res.small = res.small %>% pivot_longer(cols=c('all_except_inv', '>20k'), names_to='Size')

    # for different thresholds exploration on minisv
    res = res[, -c(6, 7)]
    # 100bp - 20kb
    res[,5] = res[,5] - res[,4]
    # 20kb - 100kb
    res[,4] = res[,4] - res[,3]
    # 100kb - 1Mb
    res[,3] = res[,3] - res[,2]
    # > 1Mb
    res[,2] = res[,2] - res[,1] 
    res = res[,-4]
    res = res %>% pivot_longer(cols=c('translocation', `>1M`, `>100k`), names_to='Size')
    return(list(res, res.small))
}

plot_bar <- function(x, is_germline=F, add_break=F, facet_wrap=F, size=12) {
    if (is_germline) {
        if (facet_wrap) {
            p=ggplot(x, aes(x=reorder(file, value), y=value, fill=Size)) + geom_bar(position='stack', stat='identity') + facet_wrap(~cell_line, ncol=5)+ylab(expression("#germline SV")) + xlab("") + get_theme(size, 45)
        } else {
            #p=ggplot(x, aes(x=reorder(file, value), y=value, fill=Size)) + geom_bar(position='stack', stat='identity') + facet_grid(genome~cell_line, scales='free')+ylab(expression("#germline SV")) + xlab("") + get_theme(size, 45)
            p=ggplot(x, aes(x=reorder(file, value), y=value, fill=Size)) + geom_bar(position='stack', stat='identity') + ylab(expression("#germline SV")) + xlab("") + get_theme(size, 45)
        }
    } else {
        if (facet_wrap) {
            p=ggplot(x, aes(x=reorder(file, value), y=value, fill=Size)) + geom_bar(position='stack', stat='identity') + facet_wrap(~cell_line, ncol=5)+ylab(expression("#mosaic SV")) + xlab("") + get_theme(size, 45)
            #p=ggplot(x, aes(x=reorder(file, value), y=value, fill=Size)) + geom_bar(position='stack', stat='identity') + ylab(expression("#mosaic SV")) + xlab("") + get_theme(size, 45)
        } else {
            p=ggplot(x, aes(x=reorder(file, value), y=value, fill=Size)) + geom_bar(position='stack', stat='identity') + facet_grid(genome~cell_line, scales='free')+ylab(expression("#mosaic SV")) + xlab("") + get_theme(size, 45) #+ scale_y_cut(breaks=c(50, 500), which=c(1, 3), scales=c(3, 0.5))
        }
    }
    p
}


plot_bar_thresholds <- function(x, size=12, add_break=F) {
    #p=ggplot(x, aes(x=count, y=value, fill=Size)) + geom_bar(position='stack', stat='identity') + facet_grid(genome~cell_line, scales='free') + ylab(expression("#mosaic SV")) + xlab("") + get_theme(size=size, angle=0)
    p=ggplot(x, aes(x=count, y=value, fill=Size)) + geom_bar(position='stack', stat='identity') + facet_wrap(~cell_line, ncol=5) + ylab(expression("#mosaic SV")) + xlab("") + get_theme(size=size, angle=0)
    #if (add_break) {
    #    p = p + scale_y_break(c(100, 1000), scales="free")
    #}
    p
}

do_bar_chart <- function(data_path, out_path, threads, myparam) {
    # R code
    df = read_tsv(data_path)
    df = preprocess_data(df)
    
    print(table(df[[1]]$Size))
    print('-----')
    print(table(df[[2]]$Size))

    pdf(out_path[['pdf']], width=9.6, height=9)
    p1=ggplot(df[[1]], aes(x=refs, y=value, fill=Size)) + geom_bar(position='stack', stat='identity') + facet_wrap(~cell_line, ncol=5) + ylab(expression("#mosaic SV")) + xlab("") + get_theme(size=15, angle=45)
    p2=ggplot(df[[2]], aes(x=refs, y=value, fill=Size)) + geom_bar(position='stack', stat='identity') + facet_wrap(~cell_line, ncol=5) + ylab(expression("#mosaic SV")) + xlab("") + get_theme(size=15, angle=45)
    print(p1/p2)
    dev.off()
}


do_bar_chart(snakemake@input[[1]], snakemake@output, snakemake@threads, snakemake@params)



