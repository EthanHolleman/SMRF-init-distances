

rule remove_hg19_header:
    input:
        'data/ucsc/ucsc.hg19.genes.tsv'
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
    awk \'{{print $2 "\t" $4 "\t" $5 "\t" $17 "\t" $1 "\t" $3}}\' {input} > {output}
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

