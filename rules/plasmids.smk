rule find_promoters:
    conda:
        '../envs/Py.yml'
    output:
        'output/promoter-locations/SS-bounds-{sample}.bed'
    params:
        fasta=lambda wildcards: df.loc[
            df_plasmids['sample_name'] == wildcards.sample]['filepath'].values[0],
        promoter_seq=lambda wildcards: df.loc[
            df_plasmids['sample_name'] == wildcards.sample]['TSS_seq'].values[0],
        promoter_direction=lambda wildcards: df.loc[
            df_plasmids['sample_name'] == wildcards.sample]['promoter_direction'].values[0]
    script:'scripts/find_promoter.py'


rule calculate_distances:
    conda:
        '../envs/Py.yml'
    input:
        promoter_bed='output/promoter-locations/SS-bounds-{sample}.bed',
        peaks=lambda wildcards: Path(
            df.loc[df_plasmids['sample_name'] == wildcards.sample]['peak_dir'].values[0]
            ).joinpath('{SMRF_sample}.bed')
    params:
        condition=lambda wildcards: wildcards.SMRF_sample
    output:
        'output/TSS-distances/{sample}/{SMRF_sample}/{sample}-{SMRF_sample}.TSS.dist.tsv'
    script:'scripts/peak2TSS.py'


rule cat_distance_files:
    input:
        distances
    output:
        'output/cat-TSS-distances.tsv'
    shell:'''
    cat {input} > {output}
    '''
