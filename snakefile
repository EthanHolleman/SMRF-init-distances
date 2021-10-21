import pandas as pd

SAMPLES = 'samples/samples.tsv'
df = pd.read_csv(SAMPLES, sep='\t').set_index('sample_name', drop=False)


def get_SMRF_samples(*args, **kwargs):
    SMRF = {}
    for each_sample in list(df['sample_name']):
        SMRF_dir = Path(df.loc[
            df['sample_name'] == each_sample]['peak_dir'].values[0]
        )
        SMRF[each_sample] = [fp.stem for fp in SMRF_dir.iterdir()]
    
    return SMRF

SMRF = get_SMRF_samples()

distances = []
for each_sample in SMRF:
    e = expand(
        expand(
        'output/TSS-distances/{sample}/{SMRF_sample}/{sample}-{SMRF_sample}.TSS.dist.tsv',
        sample=each_sample, allow_missing=True
        ), SMRF_sample=SMRF[each_sample]
    )
    distances += e


rule all:
    input:
        'output/plots/TSS-distance-plot.pdf'


rule find_promoters:
    conda:
        '../envs/Py.yml'
    output:
        'output/promoter-locations/SS-bounds-{sample}.bed'
    params:
        fasta=lambda wildcards: df.loc[
            df['sample_name'] == wildcards.sample]['filepath'].values[0],
        promoter_seq=lambda wildcards: df.loc[
            df['sample_name'] == wildcards.sample]['TSS_seq'].values[0],
        promoter_direction=lambda wildcards: df.loc[
            df['sample_name'] == wildcards.sample]['promoter_direction'].values[0]
    script:'scripts/find_promoter.py'


rule calculate_distances:
    conda:
        '../envs/Py.yml'
    input:
        promoter_bed='output/promoter-locations/SS-bounds-{sample}.bed',
        peaks=lambda wildcards: Path(
            df.loc[df['sample_name'] == wildcards.sample]['peak_dir'].values[0]
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


rule plot_distances:
    conda:
        'envs/R.yml'
    input:
        distances='output/cat-TSS-distances.tsv'
    output:
        plot='output/plots/TSS-distance-plot.pdf'
    script:'scripts/distPlot.R'




