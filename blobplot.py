#!/usr/bin/python3

import os
import sys
import argparse
import glob
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from Bio.SeqUtils import GC
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="This script constructs a scatterplot of GC content against coverage for all contigs in the final dRep quality-filtered and dereplicated bins. Additional information in this plot is displayed by colour of datapoints indicating taxonomy and size or the datapoint circles indicating length of the contig. The principle of this graph is derived from http://madsalbertsen.github.io/multi-metagenome/. This script also outputs the final dataframe containing all the information accumulated for graph construction", epilog="References:"+"\n"+"Albertsen, M., Hugenholtz, P., Skarshewski, A., Nielsen, K. L., Tyson, G. W., & Nielsen, P. H. (2013). Genome sequences of rare, uncultured bacteria obtained by differential coverage binning of multiple metagenomes. Nature biotechnology, 31(6), 533.")
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
    print("Output directory " + str(args.output_dir) + " does not exist. Creating it now.")
    setup_output_dir = "mkdir " + str(args.output_dir)
    os.system(setup_output_dir)

# Make a datframe detailing info for each binned contig including the bin it is associated with
contigs_in_bins = 0
for mag in glob.glob(str(args.bin_directory)+"/*.fasta"):
    for record in SeqIO.parse(mag, "fasta"):
        contigs_in_bins += 1

contigs_in_bins = 0
blob_df = pd.DataFrame(index = range(contigs_in_bins), columns = ["taxonomy", "mag", "header","length","coverage","gc_content"])
blob_df[:] = 0.0
blob_df = blob_df[["length","coverage","gc_content"]].astype("float")
for mag in glob.glob(str(args.bin_directory)+"/*.fasta"):
    mag_name = mag.split("/")[-1]
    for record in SeqIO.parse(mag, "fasta"):
        blob_df.at[contigs_in_bins, "mag"] = mag_name[:-6]
        blob_df.at[contigs_in_bins, "header"] = record.id
        blob_df.at[contigs_in_bins, "coverage"] = record.id.split('_')[5]
        blob_df.at[contigs_in_bins, "length"] = record.id.split('_')[3]
        blob_df.at[contigs_in_bins, "gc_content"] = GC(record.seq)
        contigs_in_bins += 1

print("The total number of all contigs in all the bins is "+str(contigs_in_bins)+" and a dataframe of length "+str(len(blob_df))+" was constructed\n")

# Read in taxonomy file (results.tsv from checkM as part of dRep)
if args.taxonomy_file_arc is not None:
    bin_taxonomy_bac = pd.read_csv(args.taxonomy_file_bac, sep = "\t", index_col=0)
    bin_taxonomy_arc = pd.read_csv(args.taxonomy_file_arc, sep = "\t", index_col=0)
    bin_taxonomy = pd.concat([bin_taxonomy_bac,bin_taxonomy_arc])
else:
    bin_taxonomy = pd.read_csv(args.taxonomy_file_bac, sep = "\t", index_col=0)
print("Taxonomy for "+str(len(bin_taxonomy))+" bins was identified from GTDGtk and will be incorporated into the output while "+str(len(blob_df["mag"].unique()))+" bins were read in from the bin_directory. These numbers should be equal.\n")

# Change taxonomic classification into usable chunk, e.g. class level or lowest ID
for i in bin_taxonomy.index:
    # Identifies at class_level
    new = bin_taxonomy.at[i,"classification"].split(';')[2]
    bin_taxonomy.at[i,"class"] = new
    # reverse list and if length is bigger than no ID, i.e. len()>f__ (3), then that is lowest ID
    tax_list = bin_taxonomy.at[i,"classification"].split(';')
    rev_tax_list = tax_list[::-1]
    for j in rev_tax_list:
        if len(j) <= 3:
            continue
        elif len(j) > 3:
            bin_taxonomy.at[i,"lowest_ID"] = j
            break
        else:
            print("weird taxonomic classification from GTDBtk")

# Merge all dataframes
blob_df_final = pd.merge(blob_df, bin_taxonomy, how='left', left_on=["mag"], right_on=bin_taxonomy.index)
print("The final dataframe contains a total of "+str(len(blob_df_final[blob_df_final["classification"].notna()==True]))+" contigs, which have taxonomy associated with them, i.e. will end up in the final graph. This should be equal to the initial dataframe length of "+str(len(blob_df))+".\n")
blob_df_final.to_csv(str(args.output_dir)+"/blob_df_final.tsv", sep='\t')

# Make a scatterplot using lowest ID and save to file
m = len(blob_df_final["lowest_ID"].unique())+1
sns.set()
plt.figure(figsize=(20, 20))
ax = sns.scatterplot(x=blob_df_final["gc_content"], y=blob_df_final["coverage"],
                     hue=blob_df_final["lowest_ID"], size=blob_df_final["length"],
                     sizes=(10, 2000), alpha=0.5, data=blob_df_final)
h,l = ax.get_legend_handles_labels()
l1 = ax.legend(h[0:m], l[0:m], loc='upper left', frameon=False)
l2 = ax.legend(h[m:], l[m:], loc='upper right', frameon=False, labelspacing=3)
ax.add_artist(l1)
plt.yscale("log")
print("An error involving backend_gtk3.py will show and tell you that some source ID could not be found when attempting to remove it. This is a known issue with matplotlib (https://github.com/matplotlib/matplotlib/issues/5561/), which was not fixed here since output is still produced as expected.\n")
plt.savefig(str(args.output_dir)+"/dereplicated_genomes_gtdbtk_lowest_ID.pdf", format='pdf')

# Make a scatterplot using class and save to file
d = len(blob_df_final["class"].unique())+1
sns.set()
plt.figure(figsize=(20, 20))
bx = sns.scatterplot(x=blob_df_final["gc_content"], y=blob_df_final["coverage"],
                     hue=blob_df_final["class"], size=blob_df_final["length"],
                     sizes=(10, 2000), alpha=0.5, data=blob_df_final)
e,f = bx.get_legend_handles_labels()
l1 = bx.legend(e[0:d], f[0:d], loc='upper left', frameon=False)
l2 = bx.legend(e[d:], f[d:], loc='upper right', frameon=False, labelspacing=3)
bx.add_artist(l1)
plt.yscale("log")
print("An error involving backend_gtk3.py will show and tell you that some source ID could not be found when attempting to remove it. This is a known issue with matplotlib (https://github.com/matplotlib/matplotlib/issues/5561/), which was not fixed here since output is still produced as expected.\n")
plt.savefig(str(args.output_dir)+"/dereplicated_genomes_gtdbtk_class.pdf", format='pdf')

