# Metagenomic singularity workflow (meta_sing)
Initial relase: 15/05/2020
Final release: 15/12/2020

## Requirements
Definitely needed:
- singularity 3.4.2-1 or higher (3.5.0 or higher is recommended due to some security-related updates of singularity)
- wget
- git

Desirable:
- Ubuntu 18.04 LTS

## General Notes
- This workflow is very much designed for reproducibility, transferability and ease-of-access rather than computational efficiency. Suggestions are as always welcome but may not be included in the git, as it is intended to allow reproducibility of work yet to be published.
- All containers were built on a local Ubuntu 18.04.4 LTS machine (Intel CORE i5 vPro 7th Gen) running singularity 3.4.2-1 and were executed on Intel nodes in a slurm-based HPC environmnent
- Containers autometa_latest_2.sif and container_5 are currently only able to run on Intel nodes due to an SSSE3 issue, which to the best of my knowledge is a binary encoding library that is related to Intel archhitecture and incompatible with AMD
- Descriptions of each container are accessible by running singularity run-help path/to/container
- Descriptions of each python script are accessible by running python3 python_script.py --help
- meta_sing.sh does not contain help-text but checks order and type of parameters
- Ultimately .sif files will hopefully be uploaded to SingularityHub or similar

## Workflow
Before running the workflow, two databases need to be downloaded, one for autometa and one for GTDB-tk. These are quite large and would be impractical to include in the .sif files. It should be noted that the SILVA database and antiSMASH database are included in container_5 since these are smaller. In order to download the databases run the following:
```
# Firstly, the autometa binaries from the git need to be cloned and fed into meta_sing as a parameter, which will bind them into the singularity container. This container was converted from a docker container rather than being built using the manual install instructions of autometa as a .def file (tried this but not successfull). When trying to run the run_autometa.py script within autometa_latest_2.sif, hmmpress throws an error which I (for the benefit of time) did not resolve. 
git clone https://github.com/KwanLab/Autometa.git

# The command doesn't actually check if fake.fasta (arbitrary name/non-existent file) exist, so it should just download the database (~163 GB) and then crash
singularity exec autometa_latest.sif path/to/autometa/pipeline/make_taxonomy_table.py -db path/to/autometa_db_directory -a fake.fasta

# This will download the GTDB-tk database (~53 GB)
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/gtdbtk_r89_data.tar.gz
tar xvzf gtdbtk_r89_data.tar.gz
```
Running ```bash meta_sing.sh``` will run all containers in sequence with checkpoints after container 1, 2 and 3. The required arguments (in order) are as below. A couple of things to note are:
1. The meta_sing.sh script currently only checks if the files and directories exist and that expected integers are actually integers but ultimately the algorithms themselves should report any errors related to file format, integrity, etc.
2. Some of the "unnecessary" files produced by algorithms are autmatically deleted throughout the script (see rm commands for details) but ultimately users should be able to delete more themselves
```
$1=container_directory $2=threads $3=memory $4=forwards_reads $5=reverse_reads $6=output_directory $7=autometa_db_directory $8=gtdbtk_directory
# Full command
bash meta_sing.sh $1 $2 $3 $4 $5 $6 $7 $8
```

Steps in the workflow:
1. container_1 uses adap_ID.sh (Matt Storey) to identify adapters and adds these to the standard Tru-Seq3-PE2.fa adapters provided by Trimmomatic (Bolger et al. 2014), then uses Trimmomatic v.0.36 to trim reads and assembles them using metaSPAdes v.3.14.0 (Nurk et al. 2017)
2. container_2 runs maxbin2 (Wu et al. 2015) and metaBAT2 (Kang et al. 2019) as part of metaWRAP v.1.2.1 (Uritskiy et al. 2018) on the assembly to produce two sets of bins/MAGs
3. autometa_latest_2.sif (commit d69eb62 from https://hub.docker.com/r/jasonkwan/autometa/builds) runs the standard Autometa (Miller et al. 2019) workflow wihout machine-learning to produce another set of bins/MAGs
4. container_4 dereplicates (and qaulity-filters) the three bin/MAG sets using dRep v.2.6.2 (Olm et al. 2017) to result in one high-quality set of MAGs, then identifies the taxonomy of these high-quality MAGs using GTDB-tk v.1.1.1 (Chaumeil et al. 2020) and then summarises results using the blobplot.py script producing a tsv file and two blobplots (inspired by Albertsen et al. 2013)
5. container_5 runs antiSMASH v.5.1.2 (Blin et al. 2019) on the assembly and then associates the BGCs with the bin_taxonomy producing a few summary bar graphs most notably with the number of BGCs found per taxonomy at the class level and a .csv file containing the information from associating BGCs with MAGs. Note that this container includes the meme-suite, so by using it you agree to the MEME license (http://meme-suite.org/doc/copyright.html)

References:
1. Albertsen, M., Hugenholtz, P., Skarshewski, A., Nielsen, K. L., Tyson, G. W., & Nielsen, P. H. (2013). Genome sequences of rare, uncultured bacteria obtained by differential coverage binning of multiple metagenomes. Nature biotechnology,31(6), 533-538.
2. Blin, K., Shaw, S., Steinke, K., Villebro, R., Ziemert, N., Lee, S. Y., ... & Weber, T. (2019). antiSMASH 5.0: updates to the secondary metabolite genome mining pipeline. Nucleic acids research, 47(W1), W81-W87
3. Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114-2120.
4. Chaumeil, P. A., Mussig, A. J., Hugenholtz, P., & Parks, D. H. (2020). GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database 
5. Kang, D. D., Li, F., Kirton, E., Thomas, A., Egan, R., An, H., & Wang, Z. (2019). MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. PeerJ, 7, e7359.
6. Miller, I. J., Rees, E. R., Ross, J., Miller, I., Baxa, J., Lopera, J., ... & Kwan, J. C. (2019). Autometa: automated extraction of microbial genomes from individual shotgun metagenomes. Nucleic acids research, 47(10), e57-e57.
7. Nurk, S., Meleshko, D., Korobeynikov, A., & Pevzner, P. A. (2017). metaSPAdes: a new versatile metagenomic assembler. Genome research, gr-213959
8. Olm, M. R., Brown, C. T., Brooks, B., & Banfield, J. F. (2017). dRep: a tool for fast and accurate genomic comparisons that enables improved genome recovery from metagenomes through de-replication. The ISME journal, 11(12), 2864.
9. Uritskiy, G. V., DiRuggiero, J., & Taylor, J. (2018). MetaWRAP—a flexible pipeline for genome-resolved metagenomic data analysis. Microbiome, 6(1), 1-13.
10. Wu, Y. W., Simmons, B. A., & Singer, S. W. (2015). MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets. Bioinformatics, 32(4), 605-607.

