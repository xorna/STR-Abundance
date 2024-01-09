# script for calculating repeat abundance chromosome/ scaffold wise
# %%#

# import json
import sys
import multiprocessing

sys.path.append('D:/ACTREC/Poster')
import common.common_functions as common_functions


# used to take PERF file and GENOME file as command-line arguments
# try:
#
#     # sys.argv is global, so any changes to it will reflect elsewhere in the code.
#     # take sys.argv into a local-scoped variable for use
#     args = sys.argv
#     perf_file = args[1]
#     genome_file = args[2]
#     output_file = args[3]
#
#     print(f"Command-line options provided: \n"
#           f"\t {perf_file} as PERF file \n"
#           f"\t {genome_file} as Genome file \n"
#           f"\t {output_file} as output JSON \n")
#
# except IndexError:
#
#     # in case of missing argument, this will be printed to STDOUT
#     raise SystemExit(f"Usage: {sys.argv[0]} <PERF file> <Genome file> <Output json file>")

# reading perf tsv file perf_data = common_functions.import_tsv(
# file="/home/gourab/Public/Tejas_2023/Perf_Runs/Human_ref/GRCh38.p14_perf_annotation.tsv") perf_data =
# common_functions.import_tsv(file="/home/gourab/Public/Tejas_2023/Perf_Runs/Human_ref/T2T-CHM13v2.0_perf_annotation
# .tsv")

# for chromosome, chr_start in chromosomes_list.items():
#
#     chr_list, chr_base_count = common_functions.chromosome_parser(chromosome=chromosome, fasta=file_path)
#
#     print("Running Repeat Counter \n")
#     rep_count_list = common_functions.repeat_counter(valid_row_num=row_number,
#                                                      data_frame=perf_data,
#                                                      curr_chr=chromosome,
#                                                      repeat_data=chromosomes_list)
#
#     print("Saving results")
#     common_functions.repeat_abundance(repeat_count_list=rep_count_list,
#                                       raw_base_counts=chr_base_count,
#                                       tsv_file="D:/ACTREC/Poster/PrimaryAnalysis/Results/RefGenome/GRCh38.p14",
#                                       chr_name=chromosome)

def run(chromosome):

    file_path = "D:/ACTREC/Review Paper/HumanRefGenome/GRCh38.p14/GenomeFasta/GCF_000001405.40_GRCh38.p14_genomic.fna"

    chr_list, chr_base_count = common_functions.chromosome_parser(chromosome=chromosome, fasta=file_path)

    print("Running Repeat Counter \n")
    rep_count_list = common_functions.repeat_counter(valid_row_num=row_number,
                                                     data_frame=perf_data,
                                                     curr_chr=chromosome,
                                                     repeat_data=chromosomes_list)

    print("Saving results")
    common_functions.repeat_abundance(repeat_count_list=rep_count_list,
                                      raw_base_counts=chr_base_count,
                                      tsv_file="D:/ACTREC/Poster/PrimaryAnalysis/Results/RefGenome/GRCh38.p14",
                                      chr_name=chromosome)


if __name__ == '__main__':

    print("\n Reading PERF TSV file")
    # perf_data = common_functions.import_tsv(file=perf_file)
    perf_data = common_functions.import_tsv(file="D:/ACTREC/Review Paper/Scripts/Repeat Density/files/T2T-CHM13v2.0_perf_annotation.tsv")

    # determining number of valid data rows
    print("\n Determining number of valid rows")
    row_number = common_functions.valid_rows(data_frame=perf_data)
    print(f"... Valid rows = {row_number}")

    chromosomes_list = common_functions.unique_elements(data_frame=perf_data,
                                                        valid_row_num=row_number,
                                                        element="Chromosome")

    print(len(chromosomes_list), list(chromosomes_list.items())[:5])
    print(len(chromosomes_list), list(chromosomes_list.items())[515:])

    pool = multiprocessing.Pool(processes=2)

    pool.map(run, list(chromosomes_list.keys()))

# %%
