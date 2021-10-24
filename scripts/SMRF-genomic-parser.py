# Create sample table from SMRF-seq genomic reads. Not used in the
# the workflow itself

import pandas as pd
from pathlib import Path

DATA_DIR = "data/PEAKS_GENOME"


def filename_to_fields(filepath):
    file_split = filepath.name.split("_")
    return {
        "filepath": str(filepath),
        "stem": Path(filepath).stem,
        "gene": file_split[1].replace("gene", ""),
        "strand": file_split[2].lower(),
    }


def main():

    files = Path(DATA_DIR).iterdir()
    df = pd.DataFrame([filename_to_fields(f) for f in files])
    df.to_csv("samples/genomic.tsv", index=False, sep="\t")


if __name__ == "__main__":
    main()
