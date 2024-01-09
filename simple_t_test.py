# this is the updated t_test file.
# Previously I was comparing one Reference Sequence RDL data with one Cell line RDL data that too chromosome wise.
# This is not a biologically sound comparision as I am effectively taking 23 chromosomes as 23 replicates which is not right.
# I must take 2 groups (RefSeq/Normal) and (CellLine/Tumor) and compare chromosome-wise and repeat-length wise.

# basically:


#               RefSeq        CellLine
# Chr1_mono     x1, x2        y1, y2, y3, y4, y5, y6
# Chr1_di       x1, x2        y1, y2, y3, y4, y5, y6
# .
# .
# .
# Chr23_hexa    x1, x2        y1, y2, y3, y4, y5, y6


# this will be my experimental setup for t-test

# Null-Hypothesis - There is no significant difference between Normal and Tumor.
# Alt-Hypothesis - There is statistically significant difference between Normal and Tumor -> accept if 'p <= 0.05'


# %%

import scipy.stats as stats
import pandas
from xlsxwriter import Workbook
import sys
import os
from glob import glob
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.cluster import hierarchy

# %%

# used to take RDL csv folder(s) cell-line and RefSeq and output folder as command line argument
try:

    # sys.argv is global, so any changes to it will reflect elsewhere in the code.
    # take sys.argv into a local-scoped variable for use
    args = sys.argv
    f1 = args[1]
    f2 = args[2]
    output_folder = args[3]
    num_files = args[4]

    print(f"Command-line options provided: \n"
          f"\t {f1} as folder holding Reference Sequence RDL data \n",
          f"\t {f2} as folder holding Cell-Line RDL data \n",
          f"\t {output_folder} as folder holding t-test output \n",
          f"\t {num_files} as how many files to process, any number between 1-25. \n\t Use 23 to process first 23 files which is equivalent to processing 23 chromosomes (leaving Y and unmapped, if any)",
          )
    

except IndexError:

    # in case of missing argument, this will be printed to STDOUT
    raise SystemExit(f"Usage: {sys.argv[0]} <RDL folder 1> <RDL folder 2> <output folder> <number of chromosomes>")


output_file = output_folder + '\\' + f1.split('\\')[-1] + '_' + f2.split('\\')[-1] + '_' + 't-test' + '.csv'

output_file = output_file.replace('\\', '/')


f1 = f1 + '\*.csv'
f2 = f2 + '\*.csv'

f1_paths = [f.replace("\\", "/") for f in glob(f1)][:int(num_files)]
f2_paths = [f.replace("\\", "/") for f in glob(f2)][:int(num_files)]

print(len(f1_paths), len(f2_paths))


# %%

# getting folder paths from input folders and saving them according to group (normal/tumor)

f1 = "D:\ACTREC\Poster\PrimaryAnalysis\Results\RefGenome"
f2 = "D:\ACTREC\Poster\PrimaryAnalysis\Results\CellLine"
output_folder = "D:\ACTREC\Poster\SecondaryAnalysis\\new_Result\Files"
image_folder = "D:\ACTREC\Poster\SecondaryAnalysis\\new_Result\Graphs"
num_files = 23

output_file = output_folder + '\\' + 't-test' + '.csv'
# output_file = output_folder + '\\' + 't-test' + '.xlsx'

# RS_files = [os.path.join(f1, entry) for entry in os.listdir(f1) if os.path.isdir(os.path.join(f1, entry))]

# CL_files = [os.path.join(f2, entry) for entry in os.listdir(f2) if os.path.isdir(os.path.join(f2, entry))]

# print(RS_files, '\n', CL_files)



# %%

def get_files_in_folder(folder_path):
    file_paths = []
    # Check if the folder exists
    if os.path.exists(folder_path):
        # Iterate through each entry in the folder
        count = 0
        for entry in os.listdir(folder_path):
            entry_path = os.path.join(folder_path, entry)
            # Check if the entry is a file
            if os.path.isfile(entry_path):
                file_paths.append(entry_path)
            
            count += 1
            if count == num_files:
                break

    else:
        print(f"The folder '{folder_path}' does not exist.")

    return file_paths

def get_files_in_folders(root_directory):
    folder_files_dict = {}
    # Check if the root directory exists
    if os.path.exists(root_directory):
        # Iterate through each entry in the root directory
        for entry in os.listdir(root_directory):
            entry_path = os.path.join(root_directory, entry)
            # Check if the entry is a directory
            if os.path.isdir(entry_path):
                folder_files_dict[entry] = get_files_in_folder(entry_path)

    else:
        print(f"The directory '{root_directory}' does not exist.")

    return folder_files_dict

# Get file paths within folders in the specified root directory

ref_seq_files = get_files_in_folders(f1)
cell_line_files = get_files_in_folders(f2)


# # Display the dictionary containing folder names as keys and lists of file paths as values
# for folder, files in cell_line_files.items():
#     print(f"Folder: {folder}")
#     print(f"Files:")
#     for file_path in files:
#         print(file_path)
#     print("---------------")
    


# %%

ref_seq_data = {}

for folder, files in ref_seq_files.items():

    ref_seq_data.update({
            folder.split('_')[0]: dict()
        })

    i = 1
    for file in files:

        file_df = pandas.read_csv(filepath_or_buffer=file)

        # ref_seq_data[folder.split('_')[0]].update(
        #     {
        #     file.split("\\")[-1].split('.')[0].split('_')[-1]: file_df['RDL'].to_list()
        #     }
        #     )
        
        ref_seq_data[folder.split('_')[0]].update({
            f'chr{i}': file_df['RDL'].to_list()
        })

        i += 1

print(ref_seq_data['GRch38.p14'])
print(ref_seq_data['T2T-CHM13'])
print(ref_seq_data['ASH1'])
print(ref_seq_data['hg01243'])
print(ref_seq_data['ASM2283312'])
print(ref_seq_data['CHM1'])


# %%

cell_line_data = {}

for folder, files in cell_line_files.items():

    cell_line_data.update({
            folder: dict()
        })

    i = 1
    for file in files:

        file_df = pandas.read_csv(filepath_or_buffer=file)

        # cell_line_data[folder].update(
        #     {
        #     file.split("\\")[-1].split('.csv')[0].split('_')[-1].split('.')[0]: file_df['RDL'].to_list()
        #     }
        #     )

        cell_line_data[folder].update({
            f'chr{i}': file_df['RDL'].to_list()
        })
    
        i += 1
        

for cl in cell_line_data:
    print(cell_line_data[cl])
    

# %%

# converting the 2 dictionaries into 2 dataframe

grp_1 = pandas.DataFrame(data=ref_seq_data)
print(grp_1)
print(grp_1.at['chr1', 'GRch38.p14'])


# %%

grp_2 = pandas.DataFrame(data=cell_line_data)
print(grp_2)
print(grp_2.at['chr1', 'DU4475'])
# %%

combined_df = pandas.concat([grp_1, grp_2], axis=1)
print(combined_df)

# %%

# list of repeat motif lengths to append to output for readability
rep_motifs = ['monomer', 'dimer', 'trimer', 'tetramer', 'pentamer', 'hexamer']

# just to check how many times the loop is running
run_num_out = 0
run_num_j = 0

# creating an output nested distionary
out_dict = dict()

for index, row in combined_df.iterrows():

    run_num_out += 1

    # print(index)

    # stats.ttest_ind(a=row)

    # print(row)

    # going 1 row deep means going into a chromosome
    out_dict.update({
        index: dict()
    })

    for j in range(6):

        run_num_j += 1

        temp_list = list()


        for i in range(len(combined_df.columns)):

            # run_num_j += 1

            # making a list to hold values of RDLs of all the samples
            temp_list.append(row.iloc[i][j])
        

        # {rep_motifs[j]} j value signifies which motif length RDL value is being taken for t_test

        print(f"\n Analyzing chromosome {index}, {rep_motifs[j]}s \n values are ref_seq:{temp_list[0:6]} cell-line:{temp_list[6:]}")

        # as stated above going once into j means taking a k-mer set to analyze
        out_dict[index].update({
            rep_motifs[j]: any
        })

        # I will plot dist_plot to check for normality, will also calculate shapiro-wilks.
        # dist-plots will have 6 plots in 1 figure

        # img_path = image_folder + '\\' + 'Normality' + '\\' + index

        # if Path(img_path).is_dir():
        #     print(f"The directory {image_folder} exists.")
        # else:
        #     print(f"The directory {image_folder} does not exist; attempting to create")
        #     Path(img_path).mkdir(parents=True, exist_ok=False)


        # fig, axes = plt.subplots(1, 2, figsize=(12, 6))

        # fig.suptitle(f'Histogram plot of {index} {rep_motifs[j]}s', fontsize=16)

        # sns.distplot(temp_list[:6], hist=True, kde=True, bins=5, color="blue", ax=axes[0])

        # axes[0].set_xlabel('Repeat Density')  # Set x-axis label for the first subplot
        # axes[0].set_ylabel('PDF')       # Set y-axis label for the first subplot
        # # axes[0].set_title(f'Histogram plot of {index} {rep_motifs}[j]s')

        # sns.distplot(temp_list[6:], hist=True, kde=True, bins=5, color="blue", ax=axes[1])

        # axes[1].set_xlabel('Repeat Density')  # Set x-axis label for the first subplot
        # axes[1].set_ylabel('PDF')       # Set y-axis label for the first subplot
        # # axes[1].set_title(f'Histogram plot of {index} {rep_motifs}[j]s')

        # image_file = img_path + '\\' + f'{rep_motifs[j]}.png'

        # plt.savefig(image_file)
        # plt.show()


        # plotting boxplots to check for equal variances between refseq and Cellline

        # Creating a boxplot using pandas

        # df = pandas.DataFrame(data=[temp_list[:2], temp_list[2:]],
        #                       index=['RefSeq', 'CellLine'])
        
        # print(df)

        # df.boxplot()

        # # Setting title and labels
        # plt.title(f'Boxplot of RefSeq and Cellline data {index}')
        # plt.xlabel('Groups')
        # plt.ylabel('Values')

        # # Show the plot
        # plt.show()

        # Performing Levene's test for variance homogeneity
        l_statistic, l_p_value = stats.levene(temp_list[:6], temp_list[6:])

        print(f"\n Levene's test statistic: {l_statistic}")
        print(f"\n P-value: {l_p_value}")

        alpha = 0.05  # Set your desired significance level

        if l_p_value < alpha:
            print("\n Variances are significantly different.")

            statstic, p_value = stats.ttest_ind(a=temp_list[:6], b=temp_list[6:], equal_var=False)

        else:
            print("\n Variances are not significantly different.")

            statstic, p_value = stats.ttest_ind(a=temp_list[:6], b=temp_list[6:], equal_var=True)
            

        # this loop works, the problem is this test doesn't, setting equal_vars flag to False works, only sometimes.
        # statstic, p_value = stats.ttest_ind(a=temp_list[0:2], b=temp_list[3:], equal_var=False)

        print(f"\n Conducting T-Test \n Test statistic: {statstic}, p-value: {p_value}")

        if p_value <= 0.05:

            # rejecting null-hypothesis

            print(f" \n There is statistically significant difference between Normal and Tumor")

            # out_dict[index][rep_motifs[j]].append(statstic)
            # out_dict[index][rep_motifs[j]].append(p_value)
            # out_dict[index][rep_motifs[j]].append(1)
            out_dict[index][rep_motifs[j]] = float(format(p_value, '.3g'))

        else:

            # accepting null-hypothesis

            print(f" \n There is no statistically significant difference between Normal and Tumor")

            # out_dict[index][rep_motifs[j]].append(statstic)
            # out_dict[index][rep_motifs[j]].append(p_value)
            # out_dict[index][rep_motifs[j]].append(0)
            out_dict[index][rep_motifs[j]] = float(format(p_value, '.3g'))
        
        print("-------------------------------")

    # break


# %%

# adding genomic view to the existing data


out_dict.update({
    "genomic": dict()
})


for i in range(6):
    # Initializing an empty list to store the sums
    sums_list = []

    # Iterating through each column
    for col in combined_df.columns:
        # Iterating through each element in the column lists
        col_sum = sum(combined_df[col].apply(lambda x: x[i]))
        sums_list.append(col_sum)

    # print(f"Sum of Column Element {i}: {sums_list}")

        

    # {rep_motifs[j]} j value signifies which motif length RDL value is being taken for t_test

    print(f"\n Analyzing genomic view of {rep_motifs[i]}s \n values are ref_seq:{sums_list[0:6]} cell-line:{sums_list[6:]}")

    out_dict['genomic'].update({
            rep_motifs[i]: any
        })

    # I will plot dist_plot to check for normality, will also calculate shapiro-wilks.
    # dist-plots will have 6 plots in 1 figure

    # img_path = image_folder + '\\' + 'Normality' + '\\' + index

    # if Path(img_path).is_dir():
    #     print(f"The directory {image_folder} exists.")
    # else:
    #     print(f"The directory {image_folder} does not exist; attempting to create")
    #     Path(img_path).mkdir(parents=True, exist_ok=False)


    # fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    # fig.suptitle(f'Histogram plot of {index} {rep_motifs[j]}s', fontsize=16)

    # sns.distplot(temp_list[:6], hist=True, kde=True, bins=5, color="blue", ax=axes[0])

    # axes[0].set_xlabel('Repeat Density')  # Set x-axis label for the first subplot
    # axes[0].set_ylabel('PDF')       # Set y-axis label for the first subplot
    # # axes[0].set_title(f'Histogram plot of {index} {rep_motifs}[j]s')

    # sns.distplot(temp_list[6:], hist=True, kde=True, bins=5, color="blue", ax=axes[1])

    # axes[1].set_xlabel('Repeat Density')  # Set x-axis label for the first subplot
    # axes[1].set_ylabel('PDF')       # Set y-axis label for the first subplot
    # # axes[1].set_title(f'Histogram plot of {index} {rep_motifs}[j]s')

    # image_file = img_path + '\\' + f'{rep_motifs[j]}.png'

    # plt.savefig(image_file)
    # plt.show()


    # Performing Levene's test for variance homogeneity
    l_statistic, l_p_value = stats.levene(sums_list[:6], sums_list[6:])

    print(f"\n Levene's test statistic: {l_statistic}")
    print(f"\n P-value: {l_p_value}")

    alpha = 0.05  # Set your desired significance level

    if l_p_value < alpha:
        print("\n Variances are significantly different.")

        statstic, p_value = stats.ttest_ind(a=sums_list[:6], b=sums_list[6:], equal_var=False)

    else:
        print("\n Variances are not significantly different.")

        statstic, p_value = stats.ttest_ind(a=sums_list[:6], b=sums_list[6:], equal_var=True)
        

    # this loop works, the problem is this test doesn't, setting equal_vars flag to False works, only sometimes.
    # statstic, p_value = stats.ttest_ind(a=temp_list[0:2], b=temp_list[3:], equal_var=False)

    print(f"\n Conducting T-Test \n Test statistic: {statstic}, p-value: {p_value}")

    if p_value <= 0.05:

        # rejecting null-hypothesis

        print(f" \n There is statistically significant difference between Normal and Tumor")

        # out_dict[index][rep_motifs[j]].append(statstic)
        # out_dict[index][rep_motifs[j]].append(p_value)
        # out_dict[index][rep_motifs[j]].append(1)
        out_dict['genomic'][rep_motifs[i]] = float(format(p_value, '.3g'))

    else:

        # accepting null-hypothesis

        print(f" \n There is no statistically significant difference between Normal and Tumor")

        # out_dict[index][rep_motifs[j]].append(statstic)
        # out_dict[index][rep_motifs[j]].append(p_value)
        # out_dict[index][rep_motifs[j]].append(0)
        out_dict['genomic'][rep_motifs[i]] = float(format(p_value, '.3g'))
    
    print("-------------------------------")

# break


# %%

# print(run_num_j)
print(out_dict)
print(out_dict.keys())

# %%

# creating output dataframe

out_df = pandas.DataFrame(data=out_dict)

out_df = out_df.T

print(out_df)

out_df.to_csv(
    path_or_buf=output_file
)



# %%

# converting csv to excel and applying conditional formatting

# Read the CSV file
input_csv = output_file
output_excel = output_folder + '\\' + 't-test' + '.xlsx'
df = pandas.read_csv(input_csv)

# Create a Pandas Excel writer using XlsxWriter as the engine
writer = pandas.ExcelWriter(output_excel, engine='xlsxwriter')

# Write the DataFrame to the Excel file
df.to_excel(writer, sheet_name='Sheet1', index=False)

# Get the xlsxwriter workbook and worksheet objects
workbook = writer.book
worksheet = writer.sheets['Sheet1']

# Create a cell format for light green background
light_green_format = workbook.add_format({'bg_color': '#C6EFCE'})


# Apply conditional formatting based on the specified logic
for idx, col in enumerate(df.columns):
    for row_idx, value in enumerate(df[col]):
        # Check if the value is numeric before applying conditional formatting
        if pandas.notnull(value) and pandas.to_numeric(value, errors='coerce') <= 0.05:
            worksheet.conditional_format(row_idx + 1, idx, row_idx + 1, idx, {'type': 'cell',
                                                                              'criteria': '<=',
                                                                              'value': 0.05,
                                                                              'format': light_green_format})

# Save the workbook
writer.close()

# %%

# creating a heatmap

# Read the CSV file
file_path = output_file
df = pandas.read_csv(file_path, index_col=0)  # Assuming the first column contains row names (chr1, chr2, ...)

# Create custom colormap
# colors = ['green', 'yellow', 'orange']
colors = ['skyblue', 'bisque', 'indianred']
colormap = mcolors.LinearSegmentedColormap.from_list('custom_cmap', colors, N=256)

# Create a heatmap using Seaborn
plt.figure(figsize=(10, 8))  # Adjust the figure size as needed
p_heatmap = sns.heatmap(df, annot=True, cmap=colormap, fmt=".2e")  # Adjust the colormap ('coolwarm') and format ('.2f') as needed

# Add a label to the rightmost color bar
p_colorbar = p_heatmap.collections[0].colorbar  # Get the colorbar object
p_colorbar.set_label('p-values \n\n lower signifies significant difference', fontsize=12, labelpad=15)

# Set labels and title
plt.xlabel('Motif Lengths', labelpad=15, fontsize=13)
plt.ylabel('Chromosomes', labelpad=15, fontsize=13)
plt.title('T-Test Values Heatmap')

# Show the plot
plt.tight_layout()

image_file = image_folder + '\\' + 't_test_heatmap' + '.png'
plt.savefig(image_file)
plt.show()


# %%

# clustering

# Read the CSV file into a DataFrame
data = pandas.read_csv(output_file, index_col=0)  # Assuming row names as index

# Perform hierarchical clustering for rows (chromosomes)
row_clustering = hierarchy.linkage(data.values, method='ward', metric='euclidean')

# Perform hierarchical clustering for columns (categories)
column_clustering = hierarchy.linkage(data.values.T, method='ward', metric='euclidean')

# Reorder rows and columns based on clustering
row_order = hierarchy.dendrogram(row_clustering, no_plot=True)['leaves']
col_order = hierarchy.dendrogram(column_clustering, no_plot=True)['leaves']

data_clustered = data.iloc[row_order, col_order]

# Create custom colormap
colors = ['steelblue', 'paleturquoise', 'indianred']
colormap = mcolors.LinearSegmentedColormap.from_list('custom_cmap', colors, N=256)

# Plot clustered heatmap
plt.figure(figsize=(10, 8))
heatmap = sns.heatmap(data_clustered, cmap=colormap)

# Add a label to the rightmost color bar
colorbar = heatmap.collections[0].colorbar  # Get the colorbar object
colorbar.set_label('p-values \n\n lower signifies significant difference', fontsize=12, labelpad=15)

plt.title('Clustered Heatmap of T-Test Values')
plt.xlabel('Motif lengths', labelpad=15, fontsize=13)
plt.ylabel('Chromosomes', labelpad=15, fontsize=13)

image_file = image_folder + '\\' + 't_test_clustered_heatmap' + '.png'
# plt.savefig(image_file)
plt.show()

 # %%

