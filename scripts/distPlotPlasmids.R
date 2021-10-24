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
        ) + scale_fill_manual(values=colors) +
        theme(text = element_text(size = 15, face='bold')) +
        labs(
            title='Distance TSS to peak start',
            x='Condition', y='Distance to T3 promoter'
        ) 
}

peaks.per.condition <- function(df){

    colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(df$condition)))
    ggplot(df, aes(x=condition, fill=condition)) +
    geom_bar(color='black') + theme_pubr() +
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
        ) + scale_fill_manual(values=colors) + facet_wrap(~name) +
        theme(text = element_text(face='bold', size = 15)) +
        labs(
            title='SMRF-seq peaks per condition', x='Condition', y='Peak count'
        ) 

}


main <- function(){

    file.path <- snakemake@input$distances
    df <- read.bed.like(file.path)
    box <- boxplot.tss.distances(df)
    bar <- peaks.per.condition(df)
    main.plot <- ggarrange(
        box, bar, nrow=1, ncol=2
    )
    ggsave(snakemake@output$plot, main.plot, dpi=600, width=18, height=12)



}

if (!interactive()){

    main()

}