from Bio.Seq import reverse_complement
from collections import defaultdict
import pandas as pd


def count_kmers(seq: str, k: int) -> dict:
    seq = seq.upper()
    length = len(seq)

    kmer_counts = defaultdict(int)

    for i in range(length - k + 1):
        kmer = seq[i : i + k]
        rev_kmer = reverse_complement(kmer)
        canonical_kmer = min(kmer, rev_kmer)  # choose smaller
        kmer_counts[canonical_kmer] += 1

    kmer_frequencies = {kmer: count / length for kmer, count in kmer_counts.items()}

    return kmer_frequencies


def process_dataframe_with_kmers(df: pd.DataFrame, k: int) -> pd.DataFrame:
    all_kmers = set()
    results = []

    for _, row in df.iterrows():
        seq = row["seq_hg38"]
        curation_status = row["curation_status"]
        kmer_freqs = count_kmers(seq, k)
        all_kmers.update(kmer_freqs.keys())  #

        result_row = {"curation_status": curation_status, **kmer_freqs}
        results.append(result_row)

    result_df = pd.DataFrame(results).fillna(0)

    for kmer in all_kmers:
        if kmer not in result_df.columns:
            result_df[kmer] = 0

    return result_df
