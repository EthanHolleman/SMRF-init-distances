from Bio import SeqIO


def write_promter_bed(*args, bed_path):
    with open(str(bed_path), "w") as handle:
        handle.write("\t".join([str(a) for a in args]))
    return bed_path


def find_promoter(promoter_seq, record):
    start = str(record.seq).lower().find(promoter_seq.lower())
    end = start + len(promoter_seq)
    return start, end


def main():

    # get inputs from snakemake interface
    fasta = snakemake.params["fasta"]
    promoter_seq = snakemake.params["promoter_seq"]
    output_path = str(snakemake.output)

    record = SeqIO.read(str(fasta), "fasta")
    promoter = find_promoter(promoter_seq, record)
    write_promter_bed(
        record.id,
        promoter[0],
        promoter[1],
        snakemake.params["promoter_direction"],
        bed_path=output_path,
    )


if __name__ == "__main__":
    main()
