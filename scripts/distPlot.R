# Visualize distances to TSS to SMRF-seq peaks
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

read.bed.like <- function(file.path){

    df <- as.data.frame(read.table(file = file.path, sep = '\t', header = FALSE))
    print(head(df))
    colnames(df) <- c(
        'name', 'start', 'end', 'read_name', 'score', 'strand', 'TSS_dist', 'condition',
        'tss_start', 'tss_end'
        )
    # remove pFC19Fixed names? Need to ask stella about that one
    df <- subset(df, name != 'pFC19FIXED')

}

gene.length.mean.dist <- function(df){

    df.agg <- aggregate(df[, c('TSS_dist', 'tss_start', 'tss_end')], list(df$name), mean)
    colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(df$condition)))

    ggplot(
        df.agg, 
        aes(x=(tss_end - tss_start), y=TSS_dist)
        ) + stat_summary(fun.data= mean_cl_normal) + 
        geom_smooth(method='lm', color='black') +
        geom_point(aes(fill=Group.1), pch=21, size=5, color='black', alpha=0.7) + 
        theme_pubr() + 
        labs(
            x='Gene length', y='Mean initiation distance', 
            fill='Gene', title='Mean initiation distance vs gene length') +
        theme(text = element_text(size = 15, face='bold')) +
        scale_fill_manual(values=colors) +
        stat_cor()

}

boxplot.tss.percent.distances <- function(df){

    df$gene.length <- df$tss_end - df$tss_start
    df$percent.init <- df$TSS_dist / df$gene.length
    colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(df$condition)))
    print('NUMBER VALUES LESS THAN 0')
    write.table(subset(df, df$percent.init < 0), file='test.tsv', quote=FALSE, sep='\t', col.names=NA)
    reg.genes <- ggplot(subset(df, !grepl('SNO', name, fixed=TRUE)), aes(x=condition, y=percent.init, fill=condition)) + 
        geom_boxplot(color='black') + theme_pubr() +
        scale_fill_manual(values=colors) +
        theme(text = element_text(size = 15, face='bold')) +
        theme(legend.position='None') +
        labs(
            title='Distance TSS to peak start as proportion of \ngene length',
            y='Proportion gene length', x=''
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust=1))
    
    sno.genes <- ggplot(subset(df, grepl('SNO', name, fixed=TRUE)), aes(x=condition, y=percent.init, fill=condition)) + 
        geom_boxplot(color='black') + theme_pubr() +
        scale_fill_manual(values=colors) +
        theme(text = element_text(size = 15, face='bold')) +
        theme(legend.position='None') +
        theme(axis.text.x = element_text(angle = 45, hjust=1)) +
        labs(y='', x='')
    
    ggarrange(reg.genes, sno.genes, nrow=1, ncol=2, widths = c(2, 0.7))
}

boxplot.tss.distances <- function(df){
    colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(df$condition)))

    tt <- df$condition
    # only include conditions with greater than 20 reads
    df <- df[df$name %in% names(tt[tt > 20]), ]


    ggplot(df, aes(x=condition, y=TSS_dist, fill=condition)) + 
        geom_boxplot(color='black') + theme_pubr() +
        theme(legend.position='None') + scale_fill_manual(values=colors) +
        theme(text = element_text(size = 15, face='bold')) +
        labs(
            title='Distance TSS to peak start',
            ,y='Distance (bp)', x=''
        )  +
        theme(axis.text.x = element_text(angle = 45, hjust=1))
}

peaks.per.condition <- function(df){

    colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(df$condition)))
    ggplot(df, aes(x=condition, fill=condition)) +
        geom_bar(color='black', size=1) + theme_pubr() +
        theme(legend.position='None') + 
        scale_fill_manual(values=colors) +
        theme(text = element_text(face='bold', size = 15)) +
        labs(
            title='Number of SMRF-seq peaks per gene', y='Peak count', x=''
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust=1))

}

# product.plot <- function(df){

#     df$ymin <- seq(from = 0.5, to = nrow(df)*2, length.out = nrow(df))
#     df$ymax <- df$ymin + 0.75

#     x_max = unique(df$template_length)
#     template_name = unique(df$template)

#     ggplot(df, 
#         aes(
#             ymin=ymin,
#             ymax=ymax,
#             xmin=product_start,
#             xmax=product_end,
#             fill=anneal_limit)
#         ) + geom_rect(color='black') + xlim(0, x_max) + theme_pubr() +
#             theme(
#                 axis.text.y = element_blank(),
#                 axis.ticks = element_blank()
#             ) +
#         scale_fill_viridis() + 
#         labs(
#             x='Template', 
#             y='PCR products', 
#             title=template_name, 
#             fill='Min primer anneal length') +
#         theme(text = element_text(size=20, face='bold'))

# }


main <- function(){

    file.path <- snakemake@input$distances
    df <- read.bed.like(file.path)
    #box <- boxplot.tss.distances(df)
    box <- boxplot.tss.distances (df)
    box.percent <- boxplot.tss.percent.distances(df)
    scatter <- gene.length.mean.dist(df)
    bar <- peaks.per.condition(df)

    main.plot <- ggarrange(
       bar, scatter,
       box, box.percent, nrow=2, ncol=2
    )
    ggsave(snakemake@output$plot, main.plot, dpi=600, width=18, height=12)



}

if (!interactive()){

    main()

}