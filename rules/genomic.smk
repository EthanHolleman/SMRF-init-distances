

rule simplify_hg19_genes:
    input:
        'data/ucsc/ucsc.hg19.genes.tsv'
    output:
        'data/ucsc/ucsc.hg19.genes.simp.tsv'
    threads: 8
    script:'../scripts/cleanUCSCgenes.py'


rule remove_hg19_header:
    input:
        'data/ucsc/ucsc.hg19.genes.simp.tsv'
    output:
        temp('data/ucsc/ucsc.hg19.genes.noheader.tsv')
    shell:'''
    tail -n +2 {input} > {output}
    '''


rule hg19_genes_to_bed:
    # hg19 genes downloaded from UCSC are not in strict bed file format
    # use awk to convert to Bed-6 format. Start and end position are the
    # transcription start and end sites.
    # output format
    # chr   txnStart    txnEnd  GeneSymbol  GeneID  strand
    input:
        'data/ucsc/ucsc.hg19.genes.noheader.tsv'
    output:
        'output/ucsc/ucsc.hg19.genes.clean.bed'
    shell:'''
    awk -F '\t' \'{{print $2 "\t" $4 "\t" $5 "\t" $17 "\t" $1 "\t" $3}}\' {input} > {output}
    '''


rule intersect_hg19_gene:
    conda:
        '../envs/bedtools.yml'
    # insertect SMRF genomic reads with Hg19 genes to get
    # TSS sites for each read
    input:
        genomic_smrf='data/PEAKS_GENOME/{genomic_smrf}.bed',
        hg19_genes='output/ucsc/ucsc.hg19.genes.clean.bed'
    output:
        'output/genomic/hg19-intersect/{genomic_smrf}-hg19-intersect.bed'
    shell:'''
    bedtools intersect -wa -wb -a {input.genomic_smrf} -b {input.hg19_genes} > {output} 
    '''


rule calculate_d_init:
    conda:
        '../envs/Py.yml'
    input:
        bed_intersect='output/genomic/hg19-intersect/{genomic_smrf}-hg19-intersect.bed'
    output:
        'output/genomic/hg19-d_init/{genomic_smrf}-hg19-intersect-d_init.tsv'
    script:'../scripts/peak2TSSGenomic.py'


rule cat_d_init_file:
    input:
        expand(
            'output/genomic/hg19-d_init/{genomic_smrf}-hg19-intersect-d_init.tsv',
            genomic_smrf=list(df_genomic['stem'])
        )
    output:
        'output/genomic/genomic-smrf-peaks-d_init-concat.tsv'
    shell:'''
    cat {input} > {output}
    '''

rule plot_genomic_init_distances:
    conda:
        '../envs/R.yml'
    input:
        distances='output/genomic/genomic-smrf-peaks-d_init-concat.tsv'
    output:
        plot='output/genomic/plots/genomic-tss-dists.pdf'
    script:'../scripts/distPlot.R'
        



