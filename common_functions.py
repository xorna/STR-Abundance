# This file will house all the common operations used by the 2 repeat density scripts
import copy
from typing import Any

import pandas
from tqdm import tqdm, trange


def import_tsv(file: str) -> pandas.DataFrame:
    """

    Dataframe object initialization using pandas from source Perf TSV file.

    """

    data_frame = pandas.read_csv(filepath_or_buffer=file,
                                 delimiter="\t",
                                 names=("Chromosome", "Repeat Start", "Repeat Stop", "Repeat Class", "Repeat Length",
                                        "Repeat Strand", "Motif Number", "Actual Repeat", "Gene Name", "Gene Start",
                                        "Gene Stop", "Strand", "Genomic Annotation", "Promoter annotation",
                                        "Distance from TSS"))

    return data_frame


def valid_rows(data_frame: pandas.DataFrame) -> int:
    """

    Determining number of valid rows in said dataframe.

    """

    row_number = data_frame[data_frame.columns[1]].count()
    return row_number


def genome_size(file_path: str) -> int:
    """

    Calculating the number of bases in relevant raw genomic data (ref/cell line).
    This function omits all the 'N' nucleotides.

    """

    print("\n\n")

    accepted_bases = ['A', 'T', 'G', 'C']
    total_base_count = 0
    buffer_size = 1000 * 1024 * 1024  # 26.104015735909343

    with open(file=file_path, mode='r') as human_genome:

        while True:

            lines = human_genome.readlines(buffer_size)

            if not lines:
                break

            for line in tqdm(lines):
                line = line.strip()
                line = line.strip('N')

                if line.startswith('>'):
                    continue

                if 'N' in line:
                    for char in line:
                        if char.upper() in accepted_bases:
                            total_base_count += 1
                else:
                    total_base_count += len(line)

    return total_base_count


def unique_elements(data_frame: pandas.DataFrame, valid_row_num: int, element: str) -> dict:
    """

    Lists unique elements (chromosomes, repeats) from the dataframe object.
    Takes pandas DataFrame, valid rows and type of element as the arguments
    Returns list containing unique elements found in the pandas DataFrame columns

    """

    element_set = set()
    repeat_dict = dict()

    # repeat_list = [data_frame.at[i, element] for i in range(valid_row_num) if
    #                (selection := data_frame.at[i, element]) not in element_set and not element_set.add(selection)]
    #
    # return repeat_list

    # for i in trange(valid_row_num):
    #
    #     if data_frame.at[i, element] not in element_set:
    #         element_set.add(data_frame.at[i, element])
    #         repeat_dict.update({data_frame.at[i, element]: i})
    #
    # return repeat_dict

    i = 0
    while i < valid_row_num:
    # for x in range(i, valid_row_num):

        x = copy.deepcopy(i)

        while data_frame.at[i, element] == data_frame.at[i + 1, element]:
            i += 1
        repeat_dict.update({data_frame.at[i, "Chromosome"]: [x, i+1]})
        i += 1

    return repeat_dict


def repeat_counter(valid_row_num: int, data_frame: pandas.DataFrame, curr_chr: str, repeat_data: dict) -> list:
    """

    This function will take the valid row count / Chromosome limit and count repeat number and repeat length motif
    length wise. For abundance measurement purposes, Chromosome limit will be used.
    It will return a list of dictionaries holding repeat number and length.

    """

    repeat_length_dict = {
        "monomer": 0,
        "dimer": 0,
        "trimer": 0,
        "tetramer": 0,
        "pentamer": 0,
        "hexamer": 0
    }

    repeat_num_dict = {
        "monomer": 0,
        "dimer": 0,
        "trimer": 0,
        "tetramer": 0,
        "pentamer": 0,
        "hexamer": 0
    }

    # for i in trange(curr_chr_start, valid_row_num):
    start = repeat_data.get(curr_chr)[0]
    end = repeat_data.get(curr_chr)[1]
    print(start, end)
    print(repeat_data.get(curr_chr))
    for i in trange(start, end):
        # print(data_frame.at[i, "Chromosome"], curr_chr, data_frame.at[i, "Actual Repeat"])

        # if data_frame.at[i, "Chromosome"] == curr_chr:

        current_repeat = data_frame.at[i, "Repeat Class"]
        current_repeat_len = data_frame.at[i, "Repeat Length"]

        # ********* convert this to switch case *********

        if len(current_repeat) == 1:
            repeat_length_dict['monomer'] += current_repeat_len
            repeat_num_dict['monomer'] += 1

        elif len(current_repeat) == 2:
            repeat_length_dict['dimer'] += current_repeat_len
            repeat_num_dict['dimer'] += 1

        elif len(current_repeat) == 3:
            repeat_length_dict['trimer'] += current_repeat_len
            repeat_num_dict['trimer'] += 1

        elif len(current_repeat) == 4:
            repeat_length_dict['tetramer'] += current_repeat_len
            repeat_num_dict['tetramer'] += 1

        elif len(current_repeat) == 5:
            repeat_length_dict['pentamer'] += current_repeat_len
            repeat_num_dict['pentamer'] += 1

        elif len(current_repeat) == 6:
            repeat_length_dict['hexamer'] += current_repeat_len
            repeat_num_dict['hexamer'] += 1

        else:
            pass

    return [repeat_length_dict, repeat_num_dict]


def repeat_abundance(repeat_count_list: list, raw_base_counts: int, tsv_file: str, chr_name: str):
    """

    This function will take the list of dictionaries holding repeat count (length & number) data, raw base count and
    output json file as its arguments.
    It will return repeat abundance data in the json file specified in the arguments

    """
    repeat_densities_length_dict = dict()

    for motif_len, value in tqdm(repeat_count_list[0].items()):
        repeat_density = (value / raw_base_counts) * 100

        repeat_densities_length_dict.update({f"{motif_len} repeat density": repeat_density})

    repeat_count_list[0] = repeat_densities_length_dict

    output = {'Repeat Density in Length': list(repeat_count_list[0].values()),
              'Repeat Number': list(repeat_count_list[1].values())}

    output_df = pandas.DataFrame(data=output, index=list(repeat_count_list[1].keys()))
    # output_df.to_csv(path_or_buf="D:/ACTREC/Poster/PrimaryAnalysis/Results/RepAbn.tsv",
    #                  sep='\t')

    output_df.to_string(buf=f"{tsv_file}/RepAbn{chr_name}.txt")


def chromosome_parser(chromosome: str, fasta: str) -> tuple[list[Any], int]:
    buffer_size = 1000 * 1024 * 1024  # 26.104015735909343
    chr_data = list()

    with open(file=fasta, mode='r') as fasta_file:

        in_chr = False
        new_chr_reached = False
        base_count = int()

        while True:

            lines = fasta_file.readlines(buffer_size)

            if not lines or new_chr_reached:
                break

            for line in tqdm(lines):

                line = line.strip()
                # line = line.strip('N')

                if line.startswith(f'>{chromosome}') and not in_chr:
                    in_chr = True
                    print(line)
                    print(f"Set in_chr to {in_chr}")
                elif line.startswith('>') and in_chr:
                    in_chr = False
                    new_chr_reached = True
                    print(line)
                    print(f"Set in_chr to {in_chr}")
                    print(f"Set new_chr_reached to {new_chr_reached}")
                    break

                if in_chr and line != '':
                    # print(f"Appending seq. line: {line}")
                    # chr_data.append(line)
                    base_count += len(line)

    return chr_data, base_count
