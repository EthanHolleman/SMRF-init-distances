# Visualize distances to TSS to SMRF-seq peaks
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

read.bed.like <- function(file.path){

    df <- as.data.frame(read.table(file = file.path, sep = '\t', header = FALSE))
    print(head(df))
    colnames(df) <- c(
        'name', 'start', 'end', 'read_name', 'score', 'strand', 'TSS_dist', 'condition'
        )
    # remove pFC19Fixed names? Need to ask stella about that one
    df <- subset(df, name != 'pFC19FIXED')

}

boxplot.tss.distances <- function(df){
    colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(df$condition)))
    ggplot(df, aes(x=condition, y=TSS_dist, fill=condition)) + 
    geom_boxplot(color='black') + theme_pubr() +
    facet_wrap(~name) + theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
        ) + scale_fill_manual(values=colors)
}

peaks.per.condition <- function(df){

    colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(df$condition)))
    ggplot(df, aes(x=condition, fill=condition)) +
    geom_bar(color='black') + theme_pubr() +
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
        ) + scale_fill_manual(values=colors) + facet_wrap(~name)

}

product.plot <- function(df){

    df$ymin <- seq(from = 0.5, to = nrow(df)*2, length.out = nrow(df))
    df$ymax <- df$ymin + 0.75

    x_max = unique(df$template_length)
    template_name = unique(df$template)

    ggplot(df, 
        aes(
            ymin=ymin,
            ymax=ymax,
            xmin=product_start,
            xmax=product_end,
            fill=anneal_limit)
        ) + geom_rect(color='black') + xlim(0, x_max) + theme_pubr() +
            theme(
                axis.text.y = element_blank(),
                axis.ticks = element_blank()
            ) +
        scale_fill_viridis() + 
        labs(
            x='Template', 
            y='PCR products', 
            title=template_name, 
            fill='Min primer anneal length') +
        theme(text = element_text(size=20, face='bold'))

}


main <- function(){

    file.path <- snakemake@input$distances
    df <- read.bed.like(file.path)
    box <- boxplot.tss.distances(df)
    bar <- peaks.per.condition(df)
    main.plot <- ggarrange(
        box, bar, nrow=1, ncol=2
    )
    ggsave(snakemake@output$plot, main.plot, dpi=600, width=16, height=10)



}

if (!interactive()){

    main()

}