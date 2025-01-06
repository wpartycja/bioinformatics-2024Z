from itertools import combinations
import os


def parse_fasta(file_path):
    organisms = set()
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith(">"):
                organism = line.split("[")[-1].rstrip("]\n")
                organisms.add(organism)
    return organisms


def get_data(dir_):
    file_organisms = {}
    for file_name in os.listdir(dir_):
        if file_name.endswith(".txt"):
            organisms = parse_fasta(os.path.join(dir_, file_name))
            file_organisms[file_name] = organisms
    return file_organisms


def find_common_organisms_subset(files, min_common=6, files_wanted=8):
    best_subset_ = []
    best_common_organisms = set()

    for sub_length in range(files_wanted, 1, -1):
        for subset in combinations(list(files.keys()), sub_length):
            common_organisms_ = set.intersection(*(files[file] for file in subset))

            if len(common_organisms_) >= min_common:
                if len(common_organisms_) > len(best_common_organisms):
                    best_subset_ = subset
                    best_common_organisms = common_organisms_
                    break
        if best_subset_:
            break

    return len(best_subset_), len(best_common_organisms), best_subset_, best_common_organisms


def combine_txt_files(input_folder, output_file):
    try:
        with open(output_file, 'w') as outfile:
            for filename in os.listdir(input_folder):
                if filename.endswith(".txt"):
                    file_path = os.path.join(input_folder, filename)
                    with open(file_path, 'r') as infile:
                        outfile.write(infile.read())
        print(f"All text files have been combined successfully into: {output_file}")
    except Exception as e:
        print(f"An error occurred: {e}")
