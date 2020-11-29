# Metagenomic singularity workflow (meta_sing)
Initial relase: 15/05/2020
This is a WIP

## Requirements
Definitely needed:
- singularity 3.4.2-1 or higher (3.5.0 or higher is recommended due to some security-related updates of singularity)
- wget
- git

Desirable:
- Ubuntu 18.04 LTS

## General Notes
- container_5.sif is still under development but antismash_5_1_ver3 runs antiSMASH v.5.1.0 (Am looking to integrate the 16S_rRNA_microbiome_analysis.py and BGC_to_bin_mapper.py into this container but some of the antiSMASH dependencies have  updated sincer I built antismash_5_1_ver3 and I'm currrently unable to build a container that runs antismash, long story short... WIP)
- All containers were built on a local Ubuntu 18.04.4 LTS machine (Intel CORE i5 vPro 7th Gen) running singularity 3.4.2-1 and were executed on raapoi
- Containers autometa_latest.sif and antismash_5_1_ver3 (/container_5) are currently only able to run on Intel nodes due to the some issue with SSSE3, which is some sort of binary encoding library that is related to Intel archhitecture and incompatible with AMD (not sure how to resolve or whether it even is resolvable)
- Descriptions of each container are accessible by running singularity run-help path/to/container
- Descriptions of each python script are accessible by running python3 python_script.py --help
- meta_sing.sh does not contain help-text but checks order and type of parameters
- Ultimately the databases and .sif files will hopefully end up in a shared folder on raapoi to save everyone downloading containers and databases separately but when that will happen is currently unknown

## Workflow
Before running the workflow, one will need to download two databases, one for autometa and one for GTDB-tk. These are quite large and would be impractical to include in the .sif filed. It should be noted that the SILVA database and antiSMASH database are however included in container_5 since these are smaller. In order to download the databases run the following:
```
# The autometa git needs to be cloned and referenced when running the autometa container due to an issue with the run_autometa.py script within the container. In short, this container was converted from a docker container rather than being a manually installed autometa (tried this but not yet successfully) and the run_autometa.py script makes hmmpress look for .hmm models in the wrong directory/outside the container (I think). 
git clone https://github.com/KwanLab/Autometa.git

# The command doesn't actually check if fake.fasta exist, so it should just download the database (~163 GB) and then crash
singularity exec autometa_latest.sif path/to/autometa/pipeline/make_taxonomy_table.py -db path/to/autometa_db_directory -a fake.fasta

# This will download the GTDB-tk database (~53 GB)
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/gtdbtk_r89_data.tar.gz
tar xvzf gtdbtk_r89_data.tar.gz
```
Running ```bash meta_sing.sh``` will run all containers in sequence with checkpoints after container 1, 2 and 3. Another checkpoint after container 4 will likely also be implemented. The required arguments (in order) are as below. A couple of things to note are:
1. The script currently only checks if the files and directories exist and that expected integers are actually integers but ultimately the algorithms are relied on to report any errors in terms of file format, integrity, etc.
2. Some of the "unnecessary" files produced by algorithms are autmatically deleted throughout the script (see rm commands for details)
```
$1=container_directory $2=threads $3=memory $4=forwards_reads $5=reverse_reads $6=output_directory $7=autometa_db_directory $8=gtdbtk_directory
# Full command
bash meta_sing.sh $1 $2 $3 $4 $5 $6 $7 $8
```

Steps in the workflow:
1. container_1 uses adap_ID.sh (Matt Storey) to identify adapters and add these to the standard Tru-Seq3-PE2.fa adapters provided by Trimmomatic, then trims the reads using these adapters and assembles them using metaSPAdes
2. container_2 runs maxbin2 and metaBAT2 as part of metaWRAP on the assembly to produce two sets of bins/MAGs
3. autometa_latest.sif (container_3) runs the standard Autometa workflow wihout machine-learning to produce another set of bins/MAGs
4. container_4 dereplicates (and qaulity-filters) the three bin/MAG sets using dRep to result in one high-quality set of MAGs, then identifies the taxonomy of these high-quality MAGs using GTDB-tk and then summarises results using the blobplot.py script producing a tsv file and two blobplots (Albertsen et al. 2013)
5. container_5 runs antiSMASH5.1 on the assembly and then associates the BGCs with the bin_taxonomy as well as running a 16S rRNA-based analysis on the assembly producing a bar graph with the relative abundances of taxonomic groups based on k-mer coverage (calculated by metaSPAdes during assembly)
