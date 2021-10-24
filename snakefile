import pandas as pd

PLASMIDS = 'samples/plasmids.tsv'
GENOMIC = 'samples/genomic.tsv'

df_plasmids = pd.read_csv(PLASMIDS, sep='\t').set_index('sample_name', drop=False)
df_genomic = pd.read_csv(GENOMIC, sep='\t').set_index('stem', drop=False)


GENOMIC = expand(
    'output/genomic/hg19-d_init/{genomic_smrf}-hg19-intersect-d_init.tsv',
    genomic_smrf=list(df_genomic['stem'])
)


def get_SMRF_samples(*args, **kwargs):
    SMRF = {}
    for each_sample in list(df_plasmids['sample_name']):
        SMRF_dir = Path(df_plasmids.loc[
            df_plasmids['sample_name'] == each_sample]['peak_dir'].values[0]
        )
        SMRF[each_sample] = [fp.stem for fp in SMRF_dir.iterdir()]
    
    return SMRF

SMRF = get_SMRF_samples()

distances = []
for each_sample in SMRF:
    e = expand(
        expand(
        'output/plasmids/TSS-distances/{sample}/{SMRF_sample}/{sample}-{SMRF_sample}.TSS.dist.tsv',
        sample=each_sample, allow_missing=True
        ), SMRF_sample=SMRF[each_sample]
    )
    distances += e


include: 'rules/genomic.smk'
include: 'rules/plasmids.smk'

rule all:
    input:
        'output.tar.gz'

rule gzip_output:
    input:
        genomic='output/plasmids/plots/plasmids-tss-dists.pdf',
        plasmids='output/genomic/plots/genomic-tss-dists.pdf'
    output:
        'output.tar.gz'
    shell:'''
    tar -zcvf output.tar.gz output/
    '''





