# %%

import sys
import pandas
from glob import glob
from tqdm import tqdm
from pathlib import Path

# %%

# used to take PERF folder and output folder as command line argument
try:

    # sys.argv is global, so any changes to it will reflect elsewhere in the code.
    # take sys.argv into a local-scoped variable for use
    args = sys.argv
    perf_folder = args[1]
    output_folder = args[2]

    print(f"Command-line options provided: \n"
          f"\t {perf_folder} as PERF input folder \n"
          f"\t {output_folder} as output folder \n")


    perf_path = Path(perf_folder)

    if perf_path.is_dir():
        print(f"The directory {perf_folder} exists.")
    else:
        raise FileNotFoundError(f"The directory {perf_folder} does not exist.")

    output_path = Path(output_folder)

    if output_path.is_dir():
        print(f"The directory {output_folder} exists.")
    else:
        print(f"The directory {output_folder} does not exist; attempting to create")
        output_path.mkdir(parents=True, exist_ok=False)

        if output_path.is_dir():
            print(f"The directory {output_folder} successfully created.")

    

except IndexError:

    # in case of missing argument, this will be printed to STDOUT
    raise SystemExit(f"Usage: {sys.argv[0]} <PERF folder> <Output folder>")

except:
    
    print("Something went wrong")


# %%

def calculate_RDL(perf_file: str, result_file: str) -> None:

    perf_data = pandas.read_csv(
        filepath_or_buffer=perf_file,
        sep='\t',
        skiprows=6, # first 6 rows contain metadata which messes up read_csv
        # header=None
    )

    RDL_genomic = {
        'monomer': 0,
        'dimer': 0,
        'trimer': 0,
        'tetramer': 0,
        'pentamer': 0,
        'hexamer': 0
    }

    for row in perf_data.iterrows():

        data = row[1]

        if len(data['repeatClass']) == 1:

            RDL_genomic.update(
                {
                    'monomer': RDL_genomic['monomer'] + (data['bases_norm'] / 10_000)
                }
            )
        
        elif len(data['repeatClass']) == 2:

            RDL_genomic.update(
                {
                    'dimer': RDL_genomic['dimer'] + (data['bases_norm'] / 10_000)
                }
            )
        
        elif len(data['repeatClass']) == 3:

            RDL_genomic.update(
                {
                    'trimer': RDL_genomic['trimer'] + (data['bases_norm'] / 10_000)
                }
            )
        
        elif len(data['repeatClass']) == 4:

            RDL_genomic.update(
                {
                    'tetramer': RDL_genomic['tetramer'] + (data['bases_norm'] / 10_000)
                }
            )
        
        elif len(data['repeatClass']) == 5:

            RDL_genomic.update(
                {
                    'pentamer': RDL_genomic['pentamer'] + (data['bases_norm'] / 10_000)
                }
            )
        
        elif len(data['repeatClass']) == 6:

            RDL_genomic.update(
                {
                    'hexamer': RDL_genomic['hexamer'] + (data['bases_norm'] / 10_000)
                }
            )

    # dividing RDL value by 10,000 as the og calculation uses MBP and doesnt multiply by 100. Thus dividing bases_norm by 1,000,000 (million) and then multiplying by 100 (in effect dividing by 10,000) will solve the problem.

    # for motif, bases_norm in RDL_genomic.items():


    # writing to file

    out_df = pandas.DataFrame(
        data=list(RDL_genomic.values()),
        columns=['RDL'],
        index=list(RDL_genomic.keys())
    )

    out_df.to_csv(path_or_buf=result_file)


# %%


perf_folder = perf_folder + '\*.tsv'

perf_file_paths = [f.replace("\\", "/") for f in glob(perf_folder)]

output_folder = output_folder.replace("\\", "/")

# %%

for path in tqdm(perf_file_paths):

    out_path = output_folder + '/' + path.split('/')[-1].split('.tsv')[0] + '.csv'
    
    calculate_RDL(perf_file=path, result_file=out_path)

# %%

# # calculation behind bases_norm
# # bases_norm is calculated by dividing the number of bases/ total repeat length by genome size in MBP
# print((perf_data.at[0, 'bases']/54826556998)*1_000_000)

# x = 54826556998/1000000
# print((perf_data.at[0, 'bases']/x))

# %%

# proving bases value = total repeat length

# sum = int()

# length_dist = perf_data.at[0, 'length_distribution'].split(';')
# print(length_dist)

# for item in length_dist:

#     sum += int(item.split('-')[0]) * int(item.split('-')[1])

# print(sum)

# if sum == perf_data.at[0, 'bases']:
#     print("hence proved")

# %%

# calculating repeat density in length, genomic view

# bases column holds total repeat lengths for the repeat motif (it will also include repeat lengths from all cyclical variations of the repeat motif; A rep motif has A and T variations).

# totalling instances per million for monomers will give Repeat density in Number (genomeic view)

# totalling bases_norm will give Repeat density in Length (genomic view)

# rep_len = 0

# for row in perf_data.iterrows():

#     data = row[1]
#     # print(data)
    
#     if len(data['repeatClass']) == 1:
#         rep_len += data['bases']
#     else:
#         pass
# # 189925741 ecoli base num
# # 54_826_556_998 HCC1187 base num
# rdl_genomic = (rep_len / 189925741) * 1_000_000

# print(rdl_genomic) # we get this value even if we total entries of bases_norm for this particular k-mer
    
# using the above concept, calculating RDL genomic is as simple as totalling bases_norm values for k-mers.

# %%
