from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO, Phylo
import matplotlib.pyplot as plt
from Bio.Seq import reverse_complement
from collections import defaultdict
import pandas as pd
import numpy as np
import random

# general part


def translate_records(file_path: str) -> list:
    records = list(SeqIO.parse(file_path, "fasta"))

    translated_records = []
    for record in records:
        normalized_seq = (
            record.seq.upper()
        )  # .replace('N', '') #.replace('X', '') # @TODO: ?
        translated_seq = normalized_seq.translate(to_stop=True)
        if len(translated_seq) == 0:
            continue
        translated_records.append(
            SeqRecord(normalized_seq, id=record.id, description=record.description)
        )

    return translated_records


def extract_experiments(file_path: str) -> tuple[pd.DataFrame]:
    experiments_df = pd.read_csv(file_path, sep="\t")
    experiments_df = experiments_df[["curation_status", "coordinate_hg38", "seq_hg38"]]
    experiments_df.dropna(inplace=True)
    experiments_df["seq_hg38"] = experiments_df["seq_hg38"].str.upper()

    positive_df = experiments_df[experiments_df["curation_status"] == "positive"]
    negative_df_v1 = experiments_df[experiments_df["curation_status"] == "negative"][
        ["curation_status", "coordinate_hg38", "seq_hg38"]
    ]

    return experiments_df, positive_df, negative_df_v1


# Get random negatives


def parse_coordinates(coord: str) -> tuple[str]:
    chrom, positions = coord.split(":")
    start, end = map(int, positions.split("-"))
    return chrom, start, end


def generate_random_regions(
    num_regions: int, lengths: int, chrom_sizes: int, positive_regions: int
) -> list:
    negative_regions = []
    positive_intervals = {(chrom, start, end) for chrom, start, end in positive_regions}

    for length in lengths:
        while True:
            chrom = random.choice(list(chrom_sizes.keys()))
            max_start = chrom_sizes[chrom] - length
            if max_start <= 0:
                continue  # Skip chromosomes too small for this length
            start = random.randint(0, max_start)
            end = start + length

            # Check for overlap with positive regions
            overlap = any(
                chrom == pos_chrom and not (end <= pos_start or start >= pos_end)
                for pos_chrom, pos_start, pos_end in positive_intervals
            )
            if not overlap:
                negative_regions.append((chrom, start, end))
                break

    return negative_regions


def extract_sequences(parsed_genome, regions):
    chrom_dict = {record.id: record.seq for record in parsed_genome}
    sequences = []

    for chrom, start, end in regions:
        seq = chrom_dict[chrom][start:end]
        seq_set = set(seq)
        if seq_set == set("X"):
            continue
        sequences.append((chrom, start, end, str(seq)))
    return sequences


def generate_random_negatives(parsed_genome, positive_df):
    chrom_sizes = {record.id: len(record.seq) for record in parsed_genome}

    positive_regions = []
    for coord in positive_df["coordinate_hg38"]:
        if coord in set(["na", "nan", np.nan]):
            continue
        chrom, start, end = parse_coordinates(coord)
        positive_regions.append((chrom, start, end))

    positive_lengths = [end - start for _, start, end in positive_regions]
    negative_regions = generate_random_regions(
        len(positive_lengths), positive_lengths, chrom_sizes, positive_regions
    )

    negative_sequences = extract_sequences(parsed_genome, negative_regions)

    filtered_negatives = [seq for seq in negative_sequences if "N" not in seq[3]]
    print(f"Generated {len(filtered_negatives)} valid negative samples.")

    to_df = []
    for chrom, start, end, seq in filtered_negatives:
        to_df.append(("negative", f"{chrom}:{start}-{end}", seq))

    negative_df_v2 = pd.DataFrame(
        to_df, columns=["curation_status", "coordinate_hg38", "seq_hg38"]
    )
    return negative_df_v2
