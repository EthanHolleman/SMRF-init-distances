# Calculate bed file of distances of SMRF-seq peaks to TSS sites.
import pandas as pd


def read_bed_like(peak_bed_path):
    try:
        return pd.read_csv(peak_bed_path, sep="\t", header=None).values.tolist()
    except Exception:
        return None


def calculate_distance_from_TSS_start(peak, TSS_orientation, TSS_start, TSS_end):
    if TSS_orientation == "-":
        # if negative orientation end of peak is where transcription would
        # have initiated
        peak.append(TSS_end - peak[2])
    else:
        peak.append(peak[1] - TSS_start)

    return peak


def add_TSS_distance(peak_list, *args):
    peaks = []
    for each_peak in peak_list:
        dist_peak = calculate_distance_from_TSS_start(each_peak, *args)
        peaks.append(dist_peak)
    return pd.DataFrame(peaks)


def main():

    promoter_coords = snakemake.input["promoter_bed"]
    peaks = snakemake.input["peaks"]
    output_path = str(snakemake.output)

    peak_list = read_bed_like(peaks)
    promoter = read_bed_like(promoter_coords)[0]  # should only be one item

    TSS_dist_peaks = add_TSS_distance(peak_list, promoter[-1], promoter[1], promoter[2])
    # append condition to each row
    condition = "_".join(snakemake.params["condition"].split("_")[1:])
    assert len(TSS_dist_peaks.columns) == 7

    TSS_dist_peaks[len(TSS_dist_peaks.columns)] = condition

    TSS_dist_peaks.to_csv(output_path, sep="\t", index=False, header=None)


if __name__ == "__main__":
    main()
