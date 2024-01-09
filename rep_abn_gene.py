# %%

# this script is written to calculate gene wise repeat abundance

# %%#

import sys
import pandas

sys.path.append('D:/ACTREC/Poster')
import common.common_functions as common_functions

# %%

# used to take PERF folder and output folder as command line argument
try:

    # sys.argv is global, so any changes to it will reflect elsewhere in the code.
    # take sys.argv into a local-scoped variable for use
    args = sys.argv
    genome_file = args[1]
    perf_file = args[2]
    out_folder = args[3]

except IndexError:

    # in case of missing argument, this will be printed to STDOUT
    raise SystemExit(f"Usage: {sys.argv[0]} <genome file> <perf file> <Output folder>")


# %%

def run(chromosome, file_path, out_path):
    

    out_path_csv = out_path + '\\' + f"gene_RepAbn{chromosome}.csv"
    out_path_csv = out_path_csv.replace('\\', '/')

    chr_list, chr_base_count = common_functions.chromosome_parser(chromosome=chromosome, fasta=file_path)

    print(chr_base_count)
    
    repeat_length_dict = {}

    start = chromosomes_list.get(chromosome)[0]
    end = chromosomes_list.get(chromosome)[1]
    print(start, end)
    print(chromosomes_list.get(chromosome))
    for i in range(start, end):

        gene_name = perf_data.at[i, "Gene Name"]
        repeat_length = perf_data.at[i, "Repeat Length"]

        if gene_name not in repeat_length_dict:
            repeat_length_dict.update({
                gene_name: repeat_length
            })
        else:
            repeat_length_dict[gene_name] += repeat_length
        
    print(repeat_length_dict)

    repeat_abn_chr_dict = {}
    repeat_abn_chr_dict.update({chromosome: dict()})

    for gene, rep_len in repeat_length_dict.items():

        repeat_abn_chr_dict[chromosome].update({
            gene: (rep_len / chr_base_count) / 100
        })
    
    print(repeat_abn_chr_dict)

    print("Saving results")
    
    out_df = pandas.DataFrame(
        data=repeat_abn_chr_dict
    )

    out_df = out_df.rename(columns = {chromosome: 'RDL'})

    print(out_df)

    out_df.to_csv(
        path_or_buf=out_path_csv
    )



# %%

if __name__ == '__main__':

    print("\n Reading PERF TSV file")
    perf_data = common_functions.import_tsv(file=perf_file)

    # determining number of valid data rows
    print("\n Determining number of valid rows")
    row_number = common_functions.valid_rows(data_frame=perf_data)
    print(f"... Valid rows = {row_number}")

    chromosomes_list = common_functions.unique_elements(data_frame=perf_data,
                                                        valid_row_num=row_number,
                                                        element="Chromosome")

    print(len(chromosomes_list), list(chromosomes_list.items())[:5])
    print(len(chromosomes_list), list(chromosomes_list.items())[515:])
    

    for chr in chromosomes_list.keys():
        run(chromosome=chr, file_path=genome_file, out_path=out_folder)


# %%
