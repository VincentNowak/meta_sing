#!/usr/bin/python3

import argparse
from Bio import SeqIO
import os
import fnmatch
import numpy as np
import pandas as pd
import matplotlib
import seaborn

parser = argparse.ArgumentParser(description="This script creates a summary of which BGCs are found in which taxonomies. It uses the contig header/name to link the BGC to the contig, which is then linked to a bin (provided that contig is contained by one of the associated bins) and consequently that bins taxonomy", epilog="References:"+"\n"+"")
parser.add_argument("-r", "--antismash_results", help="gbk summary file (e.g. contigs.gbk) containing the antiSMASH5.1 results", required=True)
parser.add_argument("-b", "--bin_directory", help="directory containing all the bin files in fasta format. Assumes files also end in .fasta", required=True)
parser.add_argument("-t", "--taxonomy_file_bac", help="GTDB-tk Bacteria taxonomy file to be used for graph construction", required=True)
parser.add_argument("-a", "--taxonomy_file_arc", help="GTDB-tk Archaea taxonomy file to be used for graph construction")
parser.add_argument("-o", "--output_dir", help="directory where all the output from this script will go", required=True)
args = parser.parse_args()

# Check if archaeal taxonomy file is provided/if it exists
if args.taxonomy_file_arc is None:
    print("No Archaea taxonomy file provided")
    print("Continuing...\n")

# Check output directory exists and otherwise create it
if os.path.isdir(args.output_dir) == False:
    print("Output directory " + str(args.output_dir) + " does not exist. Creating it now.\n")
    setup_output_dir = "mkdir " + str(args.output_dir)
    os.system(setup_output_dir)

# Count the number of regions (BGCs) and candidate clusters
F = open(args.antismash_results, 'r')
record_count = 0
region_count = 0
cand_cluster_count = 0
for record in SeqIO.parse(F, 'genbank'):
    record_count += 1
    for feature in record.features:
        if feature.type == 'region':
            region_count += 1
        if feature.type == 'cand_cluster':
            cand_cluster_count += 1
# Create a dataframe of the right size
df = pd.DataFrame(np.nan, index=range(cand_cluster_count), columns=['contig', 'contig_length', 'candidate_cluster_length', 'on_contig_edge', 'number_of_protoclusters', 'candidate_cluster_kind', 'product'], dtype=object)
print("The number of records in gene_clusters.gbk is", record_count)
print("The number of regions in gene_clusters.gbk is", region_count)
print("The number of candidate gene cluster in gene_clusters.gbk is", cand_cluster_count)

# Need to open file again to repopulate SeqIO.parse
F2 = open(args.antismash_results, 'r')
# Assign data from Genbank file based on an index count to the df created above
# Output by default is interpreted as a list and here the first, and assumed to be the only, object in the list is thus converted to str or int type
ind_count = 0
for record in SeqIO.parse(F2, 'genbank'):
    for feature in record.features:
        if feature.type == 'cand_cluster':
            df.at[ind_count, 'product'] = str(feature.qualifiers['product'][0])
            df.at[ind_count, 'on_contig_edge'] = str(feature.qualifiers['contig_edge'][0])
            df.at[ind_count, 'number_of_protoclusters'] = int(feature.qualifiers['protoclusters'][0])
            df.at[ind_count, 'contig'] = str(record.annotations['structured_comment']['antiSMASH-Data']['Original ID'])
            df.at[ind_count, 'candidate_cluster_kind'] = str(feature.qualifiers['kind'][0])
            loc = str(feature.location)
            df.at[ind_count, 'candidate_cluster_length'] = int(loc[1:-4].split(':')[1])-int(loc[1:-4].split(':')[0])
            length = record.annotations['structured_comment']['antiSMASH-Data']['Original ID']
            df.at[ind_count, 'contig_length'] = int(length.split('_')[3])
            ind_count += 1

print("The number of non-Other BGCs is", len(df[df['product'] != 'other']))
print("The number of candidate clusters not on a contig edge is", len(df[df['on_contig_edge'] == 'False']))
print("The total length of all BGCs is", df['candidate_cluster_length'].sum(), "bp")

# Count the number of bins and create a dataframe based on the number of contigs in all bins
# Note that this includes contigs < 5000 bp unlike the antiSMASH output
bin_count = 0
df_count = 0
for filename in os.listdir(args.bin_directory):
    if fnmatch.fnmatch(filename, '*.fasta'):
        bin_count += 1
        for sequence in SeqIO.parse(str(args.bin_directory)+"/"+str(filename), 'fasta'):
            df_count += 1
bin_df = pd.DataFrame(np.nan, index=range(df_count), columns=['bin', 'contig'], dtype=object)
bin_df.astype("string")
index_count = 0
for filename in os.listdir(args.bin_directory):
    if fnmatch.fnmatch(filename, '*.fasta'):    
        for sequence in SeqIO.parse(str(args.bin_directory)+"/"+str(filename), 'fasta'):
            bin_df.at[index_count, 'bin'] = filename[:-6]
            bin_df.at[index_count, 'contig'] = sequence.id
            index_count += 1

print("There are", bin_count, "bins/MAGs")
print("There are a total of", len(bin_df), "contigs in all the bins (bin_df),", len(bin_df['contig'].unique()), "of which are unique")

# Merges dataframes using left join, thus maintaining all bins and adding contigs with BGCs (excludes contigs without BGC)
bin_bgc_df = pd.merge(bin_df, df, on='contig', how='left')
print("There are a total of", len(bin_bgc_df['bin'].unique()), "bins in the bin_bgc_df")
print("There are a total of", len(bin_bgc_df['contig'].unique()), "unique contigs in the bin_bgc_df\n")
print("Outputting a tsv file outlining contigs occurring in more than one bin.\n")
doubles = bin_bgc_df[bin_bgc_df.duplicated(subset='contig', keep=False)]
doubles.to_csv(str(args.output_dir)+"/contigs_in_more_than_one_bin.tsv", sep='\t')
print("Outputting a tsv file with the final merged dataframe called 'bin_bgc_df'.\n")

# Making simple bar graph of BGCs per secondary metabolite class
print("Outputting a simple bar graph displaying the number of BGCs per secondary metabolite class.\n")
type_count_df = df['product'].value_counts()
matplotlib.pyplot.bar(type_count_df.index, type_count_df.values)
matplotlib.pyplot.xticks(rotation=90)
matplotlib.pyplot.subplots_adjust(bottom=0.25)
matplotlib.pyplot.savefig(str(args.output_dir)+"/BGCs_per_sec_met_class_.pdf", format='pdf')

# Read in taxonomy file (results.tsv from checkM as part of dRep)
if args.taxonomy_file_arc is not None:
    bin_taxonomy_bac = pd.read_csv(args.taxonomy_file_bac, sep = "\t", index_col=0)
    bin_taxonomy_arc = pd.read_csv(args.taxonomy_file_arc, sep = "\t", index_col=0)
    bin_taxonomy = pd.concat([bin_taxonomy_bac,bin_taxonomy_arc])
else:
    bin_taxonomy = pd.read_csv(args.taxonomy_file_bac, sep = "\t", index_col=0)
print("Taxonomy for "+str(len(bin_taxonomy))+" bins was identified from GTDGtk and will be incorporated into the output")

# Change taxonomic classification into usable chunk, i.e. class level here
for i in bin_taxonomy.index:
    new = bin_taxonomy.at[i,"classification"].split(";")[2]
    bin_taxonomy.at[i,"class"] = new

# Merge the dataframes by bin to add classification to all contigs and thus BGCs
final_bin_bgc_df = pd.merge(bin_bgc_df, bin_taxonomy, how='left', left_on=["bin"], right_on=bin_taxonomy.index)
final_bin_bgc_df.to_csv(str(args.output_dir)+"/final_bin_bgc_df.tsv", sep='\t')

# Count the number of BGCs per taxonomic class per secondary metabolite class and ouput graph
summary_frame = final_bin_bgc_df.groupby(["class", "product"], as_index=False)["candidate_cluster_kind"].count()
g = seaborn.catplot(x="product", y="candidate_cluster_kind", col="class", data=summary_frame, kind="bar", col_wrap=4)
for ax in g.axes.flat:
    ax.set_ylabel("# BGCs")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')
matplotlib.pyplot.subplots_adjust(bottom=0.15)
matplotlib.pyplot.savefig(str(args.output_dir)+"/BGCs_per_sec_met_class_per_tax_class.pdf", format='pdf')
