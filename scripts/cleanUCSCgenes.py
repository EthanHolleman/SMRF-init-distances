# UCSC gene list has lots of duplicates, here we want simplified version
# one entry per gene with a human readable gene symbol
from multiprocessing.pool import ThreadPool
import multiprocessing
import pandas as pd

pool = ThreadPool(multiprocessing.cpu_count())


SYMBOL_COL = "hg19.kgXref.geneSymbol"


def get_earliest_TSS(gene_symbol, genes):

    entries = genes.loc[genes[SYMBOL_COL] == gene_symbol]
    entries = entries.sort_values("hg19.knownGene.txStart")
    return list(entries.iloc[0].values)  # return first row only


def main():
    print("START")
    genes = pd.read_csv(str(snakemake.input), sep="\t")
    gene_symbols = set(genes[SYMBOL_COL])  # unique gene symbols
    args = [(symbol, genes) for symbol in gene_symbols]
    gene_reps = pool.starmap(get_earliest_TSS, args)

    gene_reps_df = pd.DataFrame(gene_reps)

    gene_reps_df.to_csv(str(snakemake.output), sep="\t", index=False)


if __name__ == "__main__":
    main()
