from peak2TSS import read_bed_like
import pandas as pd


def calculate_d_init(bed_row):

    peak_start, peak_end = bed_row[1], bed_row[2]
    tss_start, tss_end = bed_row[8], bed_row[9]
    peak_orr, gene_orr = bed_row[5], bed_row[12]
    gene_symbol = bed_row[10]

    if peak_orr == "+":
        dist = peak_start - tss_start
    else:
        dist = tss_end - peak_end

    return (
        gene_symbol,
        peak_start,
        peak_end,
        bed_row[3],
        0,
        bed_row[5],  # strand
        dist,
        gene_symbol,
        tss_start,
        tss_end
    )
    # name  start   end readname    score=0 strand  distance_to_tss gene_symbol


def d_init(bed_list):
    return pd.DataFrame([calculate_d_init(each_row) for each_row in bed_list])


def main():
    bed_like = snakemake.input["bed_intersect"]
    output = str(snakemake.output)
    bed_rows = read_bed_like(bed_like)
    if bed_rows:
        distances = d_init(bed_rows)
        distances.to_csv(output, sep="\t", index=False, header=None)
    else:
        with open(output, 'w') as output:
            output.write('')


if __name__ == "__main__":
    main()
