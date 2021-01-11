#!/bin/bash
# Arguments are as follows
#	$1=container_directory $2=threads $3=memory $4=forwards_reads $5=reverse_reads $6=output_directory $7=autometa_database_directory $8=gtdbtk_directory $9=autometa_binaries_folder

WORKDIR=$(pwd)


# should probably change this to check that all containers are there
# should ultimately incorporate stdout and stderr redirection
if [ ! -d "$1" ]; then
    printf "Argument 1 is missing or does not designate a directory\n" >&2
    exit 1
fi
if [[ ! $2 =~ ^[+-]?[0-9]+$ ]]; then
    printf "Argument 2 is missing or not an integer\n" >&2
    exit 1
fi
if [[ ! $3 =~ ^[+-]?[0-9]+$ ]]; then
    printf "Argument 3 is missing or not an integer\n" >&2
    exit 1
fi
if [ ! -f "$4" ]; then
    printf "Argument 4 is missing or does not designate a file\n" >&2
    exit 1
fi
if [ ! -f "$5" ]; then
    printf "Argument 5 is missing or does not designate a file\n" >&2
    exit 1
fi
if [ -e "$6" ]; then
    printf "The directory $6 already exists.\nEXITING\n"
    exit 1
fi
if [ ! -e "$6" ]; then
    printf "Creating output directory $6\nAll results will be stored here\n"
    mkdir $6
fi
if [ ! -d "$7" ]; then
    printf "Argument 7 is missing or does not designate a directory\n" >&2
    exit 1
fi
printf "Autometa database assumed to be at $7\n"
if [ ! -d "$7" ]; then
    printf "Argument 8 is missing or does not designate a directory\n" >&2
    exit 1
fi
printf "GTDBTK database assumed to be at $8\n"


cd $6

printf "\n#######################################\n"
printf "\nRunning container_1\n"
printf "\n#######################################\n"
singularity run --containall $1/container_1.sif $2 $3 $4 $5 $6/container_1_output

if [ ! -e "$6/container_1_output/assembly/contigs.fasta" ]; then
    printf "No contigs detected.\nAssembly failed.\nEXITING\n"
    exit 1
fi

printf "\nRemoving unpaired reads from Trimmomatic output as well as /K*, /tmp, /misc and /corrected folder from SPAdes output\n"
rm -r $6/container_1_output/assembly/K*
rm -r $6/container_1_output/assembly/corrected/
rm -r $6/container_1_output/assembly/tmp
rm -r $6/container_1_output/assembly/misc
rm $6/container_1_output/forward_trimmed_1U.fq.gz
rm $6/container_1_output/reverse_trimmed_2U.fq.gz 
printf "Also removing trimmomatic.log, since this takes up more space than expected\n"
rm $6/container_1_output/trimmomatic.log

printf "Gunzipping trimmed paired reads and renaming them to allow binning with metaWRAP\n"
gunzip $6/container_1_output/forward_trimmed_1P.fq.gz
mv $6/container_1_output/forward_trimmed_1P.fq $6/container_1_output/trimmed_P_1.fastq 
gunzip $6/container_1_output/reverse_trimmed_2P.fq.gz
mv $6/container_1_output/reverse_trimmed_2P.fq $6/container_1_output/trimmed_P_2.fastq

printf "\n#######################################\n"
printf "\nRunning container_2\n"
printf "\n#######################################\n"
singularity run --containall $1/container_2.sif $2 $3 $6/container_1_output/assembly/contigs.fasta $6/container_1_output/trimmed_P_1.fastq $6/container_1_output/trimmed_P_2.fastq $6/container_2_output

if [ ! -e "$6/container_2_output/metabat2_bins/bin.1.fa" ]; then
    printf "No bins detected.\nmetaWRAP failed.\nEXITING\n"
    exit 1
fi

printf "Removing /work_files directory"
rm -r $6/container_2_output/work_files

printf "\n#######################################\n"
printf "\nRunning container_3/autometa\n"
printf "\n#######################################\n"
singularity exec --containall $1/autometa_latest_2.sif make_taxonomy_table.py -p $2 -db $7 -a $6/container_1_output/assembly/contigs.fasta -l 1000 -o $6/container_3_output 

singularity exec --containall --bind $9:/mnt $1/autometa_latest_2.sif /mnt/pipeline/run_autometa.py --assembly $6/container_3_output/Bacteria.fasta --processors $2 --length_cutoff 1000 --taxonomy_table $6/container_3_output/taxonomy.tab -o $6/container_3_output/autometa_cluster -db $7

singularity exec --containall --bind $9:/mnt $1/autometa_latest_2.sif /mnt/pipeline/cluster_process.py --bin_table $6/container_3_output/autometa_cluster/recursive_dbscan_output.tab -f $6/container_3_output/Bacteria.fasta --output_dir $6/container_3_output/autometa_bins_final

if [ ! -e "$6/container_3_output/autometa_bins_final/cluster_summary_table" ]; then
    printf "No bins detected.\nautometa failed.\nEXITING\n"
    exit 1
fi


printf "\n#######################################\n"
printf "\nRunning container_4\n"
printf "\n#######################################\n"
printf "Copying and renaming all bins into one folder and renaming them to incorporate the binning algorithm\n"
mkdir $6/container_4_output
mkdir $6/container_4_output/bins_to_dereplicate
mkdir $6/container_4_output/checkm_tmp
cd $6/container_4_output/bins_to_dereplicate
cp $6/container_2_output/maxbin2_bins/*.fa $6/container_4_output/bins_to_dereplicate
for i in *.fa; do mv -v $i $i.maxbin2.fasta; done
cp $6/container_2_output/metabat2_bins/*.fa $6/container_4_output/bins_to_dereplicate
rm $6/container_4_output/bins_to_dereplicate/bin.unbinned.fa
for i in *.fa; do mv -v $i $i.metabat2.fasta; done
cp $6/container_3_output/autometa_bins_final/cluster_DBSCAN*.fasta $6/container_4_output/bins_to_dereplicate
cd $6

export SINGULARITY_BIND="$8:/gtdbtk_db_dir"
singularity run --containall --bind $6/container_4_output/checkm_tmp:/tmp $1/container_4.sif $2 $3 $6/container_4_output/bins_to_dereplicate $6/container_4_output/output $8 
if [ ! -e "$6/container_4_output/output/gtdbtk_output/gtdbtk.bac120.summary.tsv" ]; then
    printf "No bacterial bin taxonomy detected.\nautometa failed.\nEXITING\n"
    exit 1
fi

printf "\n#######################################\n"
printf "\nRunning container_5\n"
printf "\n#######################################\n"
mkdir $6/container_5_tmp
singularity run --containall --bind $6/container_5_tmp:/tmp $1/container_5.sif $2 $6/container_1_output/assembly/contigs.fasta $6/container_5_output $6/container_4_output/output/dRep_output/dereplicated_genomes/ $6/container_4_output/output/gtdbtk_output/gtdbtk.bac120.summary.tsv $6/container_4_output/output/gtdbtk_output/gtdbtk.ar122.summary.tsv

