import pandas as pd
import re
import numpy as np


def parse_fasta_to_dataframe(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    headers = []
    sequences = []
    current_sequence = []

    for line in lines:
        if line.startswith(">"):
            if current_sequence:
                sequences.append("".join(current_sequence))
                current_sequence = []
            headers.append(line.strip())
        else:
            current_sequence.append(line.strip())

    if current_sequence:
        sequences.append("".join(current_sequence))

    ids = []
    proteins = []
    species = []

    for header in headers:
        parts = header[1:].split(" ", 1)
        ids.append(parts[0])
        
        description = parts[1] if len(parts) > 1 else ""
        species_match = re.search(r"\[(.*?)\]", description)
        species_name = species_match.group(1) if species_match else ""
        protein_name = description.split("[")[0].strip()
        
        proteins.append(protein_name)
        species.append(species_name)

    df = pd.DataFrame({
        "ID": ids,
        "Protein": proteins,
        "Species": species,
        "Sequence": sequences
    })

    return df


def save_to_fasta(df, output_file):
    with open(output_file, 'w') as fasta_file:
        for index, row in df.iterrows():
            protein_name = row['Protein'].replace(' ', '_')
            species_name = row['Species'].replace(' ', '_')
            
            readable_header = f"{protein_name}_{species_name}"
            fasta_file.write(f">{readable_header} [ID:{row['ID']}]\n")
            fasta_file.write(row['Sequence'] + "\n")


def parse_blast_results_to_matrix(blast_file):
    columns = ["QueryID", "SubjectID", "PercentageIdentity", "AlignmentLength", 
               "Mismatches", "GapOpens", "QueryStart", "QueryEnd", 
               "SubjectStart", "SubjectEnd", "EValue", "BitScore"]
    blast_results = pd.read_csv(blast_file, sep="\t", header=None, names=columns)

    unique_ids = list(set(blast_results["QueryID"]).union(set(blast_results["SubjectID"])))

    similarity_matrix = pd.DataFrame(np.zeros((len(unique_ids), len(unique_ids))), 
                                     index=unique_ids, columns=unique_ids)

    for _, row in blast_results.iterrows():
        similarity_matrix.loc[row["QueryID"], row["SubjectID"]] = row["PercentageIdentity"]

    np.fill_diagonal(similarity_matrix.values, 100.0)

    return similarity_matrix

