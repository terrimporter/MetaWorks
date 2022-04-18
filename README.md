# MetaWorks: A Multi-Marker Metabarcode Pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4741407.svg)](https://doi.org/10.5281/zenodo.4741407)  

MetaWorks consists of a conda environment and Snakemake pipelines that are meant to be run at the command line to bioinformatically processes de-multiplexed Illumina paired-end metabarcodes from raw reads through to taxonomic assignments. MetaWorks currently supports a number of popular marker gene amplicons and metabarcodes: COI (eukaryotes), rbcL (eukaryotes, diatoms), ITS (fungi, plants), 16S (prokaryotes), 18S (eukaryotes, diatoms), 12S (fish, vertebrates), and 28S (fungi).  Taxonomic assignments are made using the RDP classifier that uses a naive Bayesian method to produce taxonomic assignments with a measure of statistical support at each rank. 

## Outline

[Overview](#overview)  

[Available workflows](#available-workflows)  

[Pipeline details](#pipeline-details)  

[Prepare your environment to run the pipeline](#prepare-your-environment-to-run-the-pipeline)    

[Quick start example using COI test data](#quick-start-example-using-COI-test-data)  

[Tutorial](#tutorial)

[Implementation notes](#implementation-notes)  

[How to cite](#how-to-cite)  

[References](#references)  

## Overview

MetaWorks comes with a conda environment file MetaWorks_v1.9.7 that should be activated before running the pipeline.  Conda is an environment and package manager (Anaconda, 2016).  The environment file contains most of the programs and dependencies needed to run MetaWorks.  If pseudogene filtering will be used, then the NCBI ORFfinder program will also need to be installed.  Additional RDP-trained reference sets may need to be downloaded if the reference set needed is not already built in to the RDP classifier (see Table 1 below).

Snakemake is a python-based workflow manager (Koster and Rahmann, 2012) and it requires three sets of files to run the any one of the workflows described described in the next section (Fig 1).


**Fig 1.  Using a conda environment helps to quickly gather programs and dependencies used by MetaWorks.**. 

<img src="/images/conda_env.png" width="500">

The configuration file is edited by the user to specify directory names, indicate the sample and read fields from the sequence filenames, and specify other required pipeline parameters such as primer sequences, marker name, and whether or not pseudogene filtering should be run.

The snakefile describes the pipeline itself and normally does not need to be edited in any way (Fig 2).  

**Fig 2. The snakefile describes the programs used, the settings, and the order in which to run commands.**  
The pipeline takes care of any reformatting needed when moving from one step to another.  In addition to the final results file, a summary of statistics and log files are also available for major steps in the pipeline.

<img src="/images/dataflow.png" width="500">

The pipeline begins with raw paired-end Illumina MiSeq fastq.gz files.  Reads are paired.  Primers are trimmed.  All the samples are pooled for a global analysis.  Reads are dereplicated, denoised, and chimeric sequences are removed producing a reference set of denoised exact sequence variants (ESVs). At this step, the pipeline diverges into several paths:  an ITS specific dataflow, a regular dataflow, and a pseudogene filtering dataflow.  For ITS sequences, flanking rRNA gene regions are removed then they are taxonomically assigned.  For the regular pipeline, the denoised ESVs are taxonomically assigned using the RDP classifier.  If a protein coding marker is being processed you have the option to filter out putative pseudogenes (Porter and Hajibabaei, 2021).  The result is a report containing a list of ESVs for each sample, with read counts, and taxonomic assignments with a measure of bootstrap support (Fig 3).

**Fig 3. The RDP classifier produces a measure of confidence for taxonomic assignments at each rank.**  
Results can be filtered by bootstrap support values to reduce false-positive assignments.  The appropriate cutoffs to use for filtering will depend on the marker/classifier used, query length, and taxonomic rank.  See links in Table 1 for classifier-specific cutoffs to ensure 95-99% accuracy.

<img src="/images/taxonomic_assignments.png" width="500">

## Available workflows

1. The **default workflow** starts with Illumina paired-end demultiplexed fastq files and generates taxonomically assigned exact sequence variants (ESVs).  An adapters.fasta file is required to identify the forward and reverse primers to remove.  An example is available in /testing/adapters.fasta .  Note that the reverse primer should be reverse-complemented in the adapters.fasta file.  Multiple primers sets for the same marker gene can be processed at the same time, E.g. COI_BR5, COI_F230R, COI_mljg.  If any of these amplicons are nested within each other, then the primers should be 'anchored' in the adapters.fasta file, see /testing/adapters_anchored.fasta .  Different marker genes should be processed separately, E.g. COI, ITS, rbcL.

```linux
# quickstart default ESV pipeline
snakemake --jobs 24 --snakefile snakefile_ESV --configfile config_ESV.yaml
```

2. The OTU workflow starts with the taxonomically assigned ESVs from the default dataflow and generates operational taxonomic units (OTUs) based on 97% sequence similarity.

```linux
# quickstart OTU pipeline
snakemake --jobs 24 --snakefile snakefile_OTU --configfile config_OTU.yaml
```

3. The global ESV workflow starts with the taxonomically assigned ESVs from the default dataflow and generates a GLOBAL set of ESV IDs consistent accross all samples *sequenced at different times* to which all ESVs will be mapped. This script may be useful when it is ideal to bioinformatically process samples one season at a time (or one trial at a time, or one year at a time) but still have a consistent set of equivalent ESV IDs project-wide to facilitate multi-season (or multi-trial, or multi-year) comparisons in downstream analyses.

```linux
# quickstart global ESV pipeline
snakemake --jobs 24 --snakefile snakefile_ESV_global --configfile config_ESV_global.yaml
```

4. The global OTU workflow starts with the taxonomically assigned ESVs from the default dataflow and generates a GLOBAL set of OTU IDs consistent across all samples *sequenced at different times* to which all denoised ESVs willb e mapped.  This script may be useful when it is ideal to bioinformatically process samples one season at a time (or one trial at a time, or one year at a time) but still have a consistent set of equivalent OTU IDs project-wide to facilitate multi-season (or multi-trial, or multi-year) comparisons in downstream analyses.

```linux
# quickstart global OTU pipeline
snakemake --jbos 24 --snakefile snakefile_OTU_global --configfile config_OTU_global.yaml
```

5. The single read workflow starts with Illumina paired-end reads that do NOT to overlap.  If you would like to process the R1 and R2 files separately this pipeline can be used with caution.  Note that if your reads overlap, they should be processed using the default pipeline because longer sequences tend to make better taxonomic assignments with higher accuracy.

```linux
# quickstart single read pipeline
snakemake --jobs 24 --snakefile snakefile_ESV_singleRead --configfile config_ESV_singleRead.yaml
```

6. The dual indexed individual sample workflow can handle dual-indexed amplicons sequenced from individuals (as opposed to bulk samples) that were pooled prior to Illumina indexing.  An extra mapping file is needed to sort out individuals and amplicon(s) and should contain the following fields: SampleID, Amplicon, Forward, Reverse.  This should be saved as a comma-separated values (csv) file.  An example of a DualIndexedSamples.csv is available in the 'testing' directory.  The SampleID column should contain the individual number associated with a particular dual-index combination.  The Forward column should contain the 5'-TagForwardPrimer sequence.  The Reverse column should contain the 5'-TagReversePrimer sequence. In our experiments, tags ranged from 5 - 6 bp.  The total TagPrimer length ranged from 25-27bp.  For maximum specificity, the cutadapt -O parameter in the configuration file should be set to the minimum TagPrimer length range, ex. 25bp in our example.

```linux
# quickstart dual indexed sample pipeline
snakemake --jobs 24 --snakefile snakefile_dualIndexedSamples --configfile config_dualIndexedSamples.yaml
```

These workflows will be updated on a regular basis so check for the latest version at https://github.com/terrimporter/MetaWorks/releases .  Note that the OTU and global ESV pipelines above are new and have only been tested with the 16S and ITS markers.

## Pipeline details

Raw paired-end reads are merged using SEQPREP v1.3.2 from bioconda (St. John, 2016).  This step looks for a minimum Phred quality score of 13 in the overlap region, requires at least a 25bp overlap.  These paramters are adjustable.  Using Phred quality cutoff of 13 at this step is actually more stringent than using a Phred quality score cutoff of 20 at this step as more bases will exceed the cutoff when aligning the paired reads and more mismatches (if present) are counted.

Primers are trimmed in two steps using CUTADAPT v3.2 from bioconda (Martin, 2011).  This step now uses the linked adapter approach to remove forward and reverse primers in one step.  Primer sequences need to be specified in an adapters.fasta file and the user may wish to anchor them or not, see CUTADAPT manual for details https://cutadapt.readthedocs.io/en/stable/guide.html?highlight=linked#linked-adapters .  At this step, CUTADAPT looks for a minimum Phred quality score of at least 20 at the ends, no more than 10% errors allowed in the primer, no more than 3 N's allowed in the rest of the sequence, trimmed reads need to be at least 150 bp, untrimmed reads are discarded.  Each of these parameters are adjustable. 

Files are reformatted and samples are combined for a global analysis.

Reads are dereplicated (only unique sequences are retained) using VSEARCH v2.15.2 from bioconda (Rognes et al., 2016).

Zero-radius OTUS (Zotus) are generated using VSEARCH with the unoise3 algorithm (Edgar, 2016).  They are similar to OTUs delimited by 100% sequence similarity, except that they are also error-corrected and rare clusters are removed.  Here, we define rare sequences to be sequence clusters containing only one or two reads (singletons and doubletons) and these are also removed as 'noise'.  We refer to these sequence clusters as exact sequence variants (ESVs) for clarity.  For comparison, the error-corrected sequence clusters generated by DADA2 are referred to as either amplicon sequence variants (ASVs) or ESVs in the literature (Callahan et al., 2016).  Putative chimeric sequences are then removed using the uchime3_denovo algorithm in VSEARCH.

An ESV x sample table that tracks read number for each ESV is generated with VSEARCH using --search_exact .  Note that this in this pipeline this is just an intermediate file and that the final retained set of ESVs and mapped read numbers should be retrieved from the final output file (results.csv).

For ITS, the ITSx extractor is used to remove flanking rRNA gene sequences so that subsequent analysis focuses on just the ITS1 or ITS2 spacer regions (Bengtsson-Palme et al., 2013).

For the standard pipeline (ideal for rRNA genes) performs taxonomic assignments using the Ribosomal Database classifier v2.13 (RDP classifier).  We have provided a list of RDP-trained classifiers that can be used with MetaWorks (Table 1). 

**Table 1.  Where to download trained RDP classifiers for a variety of popular marker gene/metabarcoding targets.**

| Marker | Target taxa | Classifier availability |
| ------ | ----------- | ----------------------- |
| COI | Eukaryotes | https://github.com/terrimporter/CO1Classifier |
| rbcL | Diatoms | https://github.com/terrimporter/rbcLdiatomClassifier |
| rbcL | Eukarytoes | https://github.com/terrimporter/rbcLClassifier |
| 12S | Fish | https://github.com/terrimporter/12SfishClassifier |
| 12S | Vertebrates | https://github.com/terrimporter/12SvertebrateClassifier |
| SSU (18S) | Diatoms | https://github.com/terrimporter/SSUdiatomClassifier |
| SSU (16S) | Vertebrates | https://github.com/terrimporter/16SvertebrateClassifier |  
| SSU (18S) | Eukaryotes | https://github.com/terrimporter/18SClassifier |
| SSU (16S) | Prokaryotes | Built-in to the RDP classifier |  
| ITS | Fungi (Warcup) | Built-in to the RDP classifier |  
| ITS | Fungi (UNITE 2014) | Built-in to the RDP classifier |  
| ITS | Fungi (UNITE 2021) | https://github.com/terrimporter/UNITE_ITSClassifier |
| ITS | Plants (PLANiTS) | https://github.com/terrimporter/PLANiTS_ITSClassifier |
| LSU | Fungi | Built-in to the RDP classifier |

If you are using the pipeline on a protein coding marker that does not yet have a HMM.profile, such as rbcL, then ESVs are translated into every possible open reading frame (ORF) on the plus strand using ORFfinder v0.4.3 (Porter and Hajibabae, 2021).  The longest ORF (nt) is reatined for each ESV.  Putative pseudogenes are removed as outliers with unusually small/large ORF lengths.  Outliers are calcualted as follows:  Sequence lengths shorter than the 25th percentile - 1.5\*IQR (inter quartile range) are removed as putative pseudogenes (ore sequences with errors that cause a frame-shift).  Sequence lengths longer than the 75th percentile + 1.5\*IQR are also removed as putative pseudogenes.

If you use the pipeline on a protein coding marker that has a HMM.profile, such as COI arthropoda, then the ESVs are translated into nucleotide and amino acid ORFs using ORFfinder, the longest orfs are retained and consolidated.  Amino acid sequences for the longest ORFs are used for profile hidden Markov model sequence analysis using HMMER v3.3.2.  Sequence bit score outliers are identified as follows: ORFs with scores lower than the 25th percentile - 1.5\*IQR (inter quartile range) are removed as putative pseudogenes.  This method should remove sequences that don't match well to a profile HMM based on arthropod COI barcode sequences mined from public data at BOLD.

The final output file is results.csv and it has been formatted to specify ESVs from each sample, read counts, ESV/ORF sequences, and column headers.  Additional statistics and log files are also provided for each major bioinformatic step.

*Note 1*: If you choose to further cluster denoised ESVs into OTUs, then this is done using the cluster_smallmem method in VSEARCH with id set to 0.97.  Primer trimmed reads are then mapped to OTUs to create an OTU table using the usearch_global method with id set to 0.97.  The remainder of the pipeline proceeds as described above.  The final output file is results_OTU.csv.

*Note 2*: If you have samples sequenced at diffrent times (multiple seasons, years, or trials), you will likely process these samples right after sequencing resulting in multiple sets of ESVs.  To facilitate downstream processing, it may be advantagous to have a GLOBAL set of ESV ids so that data can be compared at the ESV level accross seasons, years, or trials.  The directories containing the output files processed using the default pipeline need to be in the same directory as the snakefile script. Ex. 16S_trial1, 16S_trial2, 16S_trial3.  In each of these directoreis, there is a cat.denoised.nonchimeras file and a results.csv file.  The denoised (chimera-free) ESVs (or ITSx processed ESVs) are concatenated into a single file, dereplicated, relabelled using the SHA1 method.  This file then becomes the new global ESV sequence database.  A fasta file is generated from the results.csv file and these sequences are clustered with the new global ESV sequence database using the usearch_global methods with the id set to 1.0.  The global ESV that each trial ESV clusters with is parsed and mapped to the final output file called global_results.csv.

*Note 3*: If you have samples sequenced at diffrent times (multiple seasons, years, or trials), you will likely process these samples right after sequencing resulting in multiple sets of ESVs.  To facilitate downstream processing, it may be advantagous to have a GLOBAL set of OTU ids so that data can be compared at the ESV level accross seasons, years, or trials.  The directories containing the output files processed using the default pipeline need to be in the same directory as the snakefile script. Ex. 16S_trial1, 16S_trial2, 16S_trial3.  In each of these directoreis, there is a cat.denoised.nonchimeras file and a results.csv file.  The denoised (chimera-free) ESVs (or ITSx processed ESVs) are concatenated into a single file, dereplicated, relabelled using the SHA1 method, then clustered into OTUs with 97% sequence similarity.  This file then becomes the new global OTU sequence database.  A fasta file is generated from the results.csv file and these sequences are clustered with the new global OTU sequence database using the usearch_global methods with the id set to 0.97 .  The global OTU that each trial ESV clusters with is parsed and mapped to the final output file called global_results.csv.

*Note 4*: If you are using the single read pipeline, the read pairing step with SeqPrep is skipped.  If processing the R1 read, raw read statistics are calculated, then the primer is trimmed using CUTADAPT as described above.  If processing the R2 read, raw read statistics are calculated then the read is reverse-complemented before the primer is trimmed using CUTADAPT as described above.  The final file is results.csv.  When assessing the quality of your taxonomic assignments, be sure to use the appropriate bootstrap support cutoffs for these shorter than usual sequences.

*Note 5*: If you are sorting out tagged individuals, it is advised to set the CUTADAPT -O parameter to maximize the tag+primer overlap with the read.  For example, if the range of tag+primer lengths is 25 - 27bp, then set -O 25 (default is 3 bp).  If the tags are short, then it is especially important to ensure no mismatches between this and the read.  For this reason, it is recommended to set the CUTADAPT -e parameter to allow for zero errors between the tag+primer and read, -e 0 (default is 0.1 or 10% errors).

## Prepare your environment to run the pipeline

1. This pipeline includes a conda environment that provides most of the programs needed to run this pipeline (SNAKEMAKE, SEQPREP, CUTADAPT, VSEARCH, etc.).

```linux
# Create the environment from the provided environment.yml file.  Only need to do this step once.
conda env create -f environment.yml

# Activate the environment
conda activate MetaWorks_v1.9.7

# On the GPSC activate using source
source ~/miniconda/bin/activate MetaWorks_v1.9.7
```

2. The RDP classifier comes with the training sets to classify 16S, fungal LSU or ITS rDNA.  To classify other markers using custom-trained RDP sets, obtain these from GitHub using Table 1 as a guide .  Take note of where the rRNAclassifier.properties file is as this needs to be added to the config.yaml .

```linux
RDP:
    t: "/path/to/CO1Classifier/v4/mydata_trained/rRNAClassifier.properties"
```

3. If doing pseudogene filtering, then download and install the NCBI ORFfinder

The pipeline requires ORFfinder 0.4.3 available from the NCBI at ftp://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/ .  This program should be downloaded, made executable, and put in your conda environment path (ex. ~/miniconda/envs/MetaWorks_v1.9.7/bin).

```linux
# go to your conda environment bin
cd ~/miniconda3/envs/MetaWorks_v1.9.7/bin/.

# download ORFfinder
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/ORFfinder.gz

# decompress
gunzip ORFfinder.gz

# make executable
chmod a+x ORFfinder
```

Run the program to test that it works:

```linux
ORFfinder
```

If you get an error that requries newer GLIBC libraries (libc.so.6) or if you get an error that requires libnghttp2 then follow the instructions at [Add libraries for ORFfinder](#add-libraries-for-orffinder).

4. In most cases, your raw paired-end Illumina reads can go into a directory called 'data' which should be placed in the same directory as the other files that come with this pipeline.

```linux
# Create a new directory to hold your raw data
mkdir data
```

5. Please go through the config.yaml file and edit directory names, filename patterns, etc. as necessary to work with your filenames.

6. Run Snakemake by indicating the number of jobs or cores that are available to run the whole pipeline.  

```linux
snakemake --jobs 24 --snakefile snakefile --configfile config.yaml
```

7. When you are done, deactivate the conda environment:

```linux
conda deactivate
```

## Quick start example using COI test data

This MetaWorks quick start assumes that you have already installed the CO1 Classifier available from https://github.com/terrimporter/CO1Classifier .  If you have not already done so, load the CO1 Classifier by following the quickstart instructions on the CO1 Classifier GitHub page.  It also assumes that you already have conda installed, otherwise see [Prepare your environment to run the pipeline](#prepare-your-environment-to-run-the-pipeline).

```linux
# download latest version of MetaWorks
wget https://github.com/terrimporter/MetaWorks/releases/download/v1.9.7/MetaWorks1.9.5.zip
unzip MetaWorks1.9.5.zip
cd MetaWorks1.9.5

# edit config_testing_COI_data.yaml file to customize path to CO1v4 classifier properties file (line 131)

# create the latest conda environment and activate it
conda env create -f environment.yml
conda activate MetaWorks_v1.9.7

# run the pipeline on the COI test data
snakemake --jobs 1 --configfile config_testing_COI_data.yaml --snakefile snakefile_ESV
```

Once you have the results.csv file, results can be imported into R for further analysis.  If you need an ESV table for downstream analysis this can be generated in R as well [How to create an ESV table](#how-to-create-an-esv-table).

## Tutorial

We have provided a small set of COI paired-end Illumina MiSeq files for this tutorial.  These sequence files contain reads for several pooled COI amplicons, but here we will focus on the COI-BR5 amplicon (Hajibabaei et al., 2012, Gibson et al., 2014).  Same as the quickstart above, but with additional instructions here if needed.

**Step 1.  Prepare your environment for the pipeline.**

Begin by downloading the latest MetaWorks release available at https://github.com/terrimporter/MetaWorks/releases/tag/v1.9.7 by using wget from the command line:

```linux
# download the pipeline
wget https://github.com/terrimporter/MetaWorks/releases/download/v1.9.7/MetaWorks1.9.5.tar.gz

# unzip the pipeline
unzip MetaWorks1.9.5.zip
```

If you don't already have conda on your system, then you will need to install it:

```linux
# Download miniconda3
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Install miniconda3 and initialize
sh Miniconda3-latest-Linux-x86_64.sh

# Add conda to your PATH
# If you do not have a bin folder in your home directory, then create one first.
mkdir ~/bin

# Enter bin folder
cd ~/bin

# Create symbolic link to conda
ln -s ~/miniconda3/bin/conda conda
```

Create then activate the MetaWorks_v1.9.7 environment:

```linux
# Move into the MetaWorks folder
cd MetaWorks1.9.5

# Create the environment from the provided environment.yml file .  Only need to do this step once.
conda env create -f environment.yml

# Activate the environment.  Do this everytime before running the pipeline.
conda activate MetaWorks_v1.9.7

```

To taxonomically assign COI metabarocodes, you will  need to install the RDP-trained COI Classifier from https://github.com/terrimporter/CO1Classifier/releases/tag/v4 .  You can do this at the command line using wget.

```linux
# download the COIv4 classifier
wget https://github.com/terrimporter/CO1Classifier/releases/download/v4/CO1v4_trained.tar.gz

# decompress the file
tar -xvzf CO1v4_trained.tar.gz

# Note the full path to the rRNAClassifier.properties file, ex. mydata_trained/rRNAClassifier.properties
```

If you wish to filter out putative pseudogenes, if you do not already have the NCBI ORFfinder installed on your system, then download it from the NCBI at ftp://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/ .  You can download it using wget then make it executable in your path:

```linux
# download
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/ORFfinder.gz

# decompress
gunzip ORFfinder.gz

# make executable
chmod 755 ORFfinder

# If you do not have a bin folder in your home directory, then create one first.
mkdir ~/bin

# put in your PATH (ex. ~/bin).
mv ORFfinder ~/bin/.
```

**Step 2.  Run MetaWorks using the COI testing data provided.**

The config_testing_COI_data.yaml file has been 'preset' to work with the COI_data files in the testing folder.  You will, however, still need to add the path to the trained COI classifier and save your changes.

```linux
RDP:
# If you are using a custom-trained reference set 
# enter the path to the trained RDP classifier rRNAClassifier.properties file here:
    t: "/path/to/CO1Classifier/v4/mydata_trained/rRNAClassifier.properties"
```    
  
Then you should be ready to run the MetaWorks pipeline on the testing data.
    
```linux
# You may need to edit the number of jobs you would like to run, ex --jobs 1 or --jobs 4, according to how many cores you have available
snakemake --jobs 2 --snakefile snakefile --configfile config_testing_COI_data.yaml
```

**Step 3.  Analyze the output.**
The final output file is called results.csv .  The results are for the COI-BR5 amplicon.  This can be imported into R for bootstrap support filtering, pivot table creation, normalization, vegan analysis, etc.  There are also a number of other output files in the stats directory showing the total number of reads processed at each step as well as the sequence lengths.  Log files are also available for the dereplication, denoising, and chimera removal steps.

If you are done with MetaWorks, deactivate the conda environment:

```linux
conda deactivate
```

## Implementation notes

### Installing Conda

Conda is an open source package and environment management system.  Miniconda is a lightweight version of conda that only contains conda, python, and their dependencies.  Using conda and the environment.yml file provided here can help get most of the necessary programs in one place to run this pipeline.  Snakemake is a Python-based workflow management tool meant to define the rules for running this bioinformatic pipeline.  There is usually no need to edit the snakefile file directly.  Changes to select parameters can be made in the config.yaml file.  If you install conda and activate the environment provided, then you will also get the correct versions of the open source programs used in this pipeline including Snakemake.

Install miniconda as follows:

```linux
# Download miniconda3
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Install miniconda3 and initialize
sh Miniconda3-latest-Linux-x86_64.sh

# Add conda to your PATH, ex. to ~/bin
cd ~/bin
ln -s ~/miniconda3/bin/conda conda

# Activate conda method 1 (working in a container)
source ~/miniconda3/bin/activate MetaWorks_v1.9.7

# Activate conda method 2
conda activate MetaWorks_v1.9.7
```

### Checking program versions

Ensure the program versions in the environment are being used.

```linux
# Create conda environment from file.  Only need to do this once.
conda env create -f environment.yml

# activate the environment
conda activate MetaWorks_v1.9.7

# list all programs available in the environment at once
conda list > programs.list

# you can check that key programs in the conda environment are being used (not locally installed versions)
which SeqPrep
which cutadapt
which vsearch
which perl

# you can also check their version numbers one at a time instead of running 'conda list'
cutadapt --version
vsearch --version
```

Version numbers are also tracked in the snakefile.

### Add libraries for ORFfinder

#### Resolving GLIBC errors

If you have an older version of GLIBC (on Centos6), then you may be missing the libc.so.6 file that ORffinder needs.  This file is already available, but you need to create some links yourself.

Create a symbolic link to the library:

```linux
cd ~/miniconda3/envs/MetaWorks_v1.9.7/lib
ln -s ../glibc-2.14/lib/libc.so.6 libc.so.6
```

Create the shell script file LD_PATH.sh in the following location to set the environment variable: ~/miniconda3/envs/MetaWorks_v1.9.7/etc/conda/activate.d/LD_PATH.sh

Put the following text in the LD_PATH.sh file:

```linux
export LD_LIBRARY_PATH_CONDA_BACKUP=$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
```

Create the file LD_PATH.sh in the following location to unset the environment variable:  
~/miniconda3/envs/MetaWorks_v1.9.7/etc/conda/deactivate.d/LD_PATH.sh

Put the following text in the LD_PATH.sh file:

```linux
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_CONDA_BACKUP
```

Deactivate then reactivate the environment.

Test ORFfinder:

```linux
ORFfinder
```

#### Resolving libnghttp2 errors

The missing library (on Centos7), is already available, only the LD_LIBRARY_PATH needs to be updated as follows:

Create the shell script file LD_PATH.sh in the following location to set the environment variable: ~/miniconda3/envs/MetaWorks_v1.9.7/etc/conda/activate.d/LD_PATH.sh

Put the following text in the LD_PATH.sh file:

```linux
export LD_LIBRARY_PATH_CONDA_BACKUP=$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
```

Create the file LD_PATH.sh in the following location to unset the environment variable:  
~/miniconda3/envs/MetaWorks_v1.9.7/etc/conda/deactivate.d/LD_PATH.sh

Put the following text in the LD_PATH.sh file:

```linux
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_CONDA_BACKUP
```

Deactivate then reactivate the environment.

Test ORFfinder:

```linux
ORFfinder
```

### Use the bold.hmm with the pseudogene removal pipeline

Ensure that bold.hmm as well as the indexed files (bold.hmm.h3p, bold.hmm.h3m, bold.hmm.h3i, and bold.hmm.h3f) are available in the same directory as you snakefile.

### Batch renaming of files

Sometimes it necessary to rename raw data files in batches.  I use Perl-rename (Gergely, 2018) that is available at https://github.com/subogero/rename not linux rename.  I prefer the Perl implementation so that you can easily use regular expressions.  I first run the command with the -n flag so you can review the changes without making any actual changes.  If you're happy with the results, re-run without the -n flag.

```linux
rename -n 's/PATTERN/NEW PATTERN/g' *.gz
```

### Symbolic links

Symbolic links are like shortcuts or aliases that can also be placed in your ~/bin directory that point to files or programs that reside elsewhere on your system.  So long as those scripts are executable (e.x. chmod 755 script.plx) then the shortcut will also be executable without having to type out the complete path or copy and pasting the script into the current directory.

```linux
ln -s /path/to/target/directory shortcutName
ln -s /path/to/target/directory fileName
ln -s /path/to/script/script.sh commandName
```

### Prevent errors caused by lost network connections

Large jobs with many fastq.gz files can take a while to run.  To prevent unexpected network connection interruptions from stopping the MetaWorks pipeline, set the job up to keep running even if you disconnect from your session:

1. You can use nohup (no hangup).

```linux
nohup snakemake --jobs 24 --snakefile snakefile --configfile config.yaml
```

2. You can use screen.

```linux
# to start a screen session
screen
ctrl+a+c
conda activate MetaWorks_v1.9.7
snakemake --jobs 24 --snakefile snakefile --configfile config.yaml
ctrl+a+d

# to list all screen session ids, look at top of the list
screen -ls

# to re-enter a screen session to watch job progress or view error messages
screen -r session_id
```

3. You can submit your job to the queuing system if you use one.

### How to create an ESV table

The results.csv file is the final file in the MetaWorks pipeline.  An ESV table (also known as an OTU table or pivot table) can be created in Excel or R.  In R, changing the results.csv file from a wider to a longer tabular format can be done using different libraries.

```linux
# R using reshape2 library
ESVtable <- reshape2::dcast(df, SampleName ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
```

### How to filter pseudogenes for taxa with different genetic codes

If you are targeting a broad group, such as Metazoa using COI primers, you can still filter out pseudogenes using removal method 1 that uses ORFfinder.  This can be done in two steps, for example by first processing invertebrate phyla, then processing phylum chordata that includes the vertebrata clade (see NCBI taxonomy).  Note that pseudogene removal method 2 that uses HMMer is currently only available for COI arthropoda at this time.

1. Edit the config_ESV.yaml file as follows:  
Line 156, set the taxonomy filter to '-e Annelida -e Arthropoda -e Bryozoa -e Cnidaria -e Echinodermata -e Mollusca -e Nematoda -e Nematomorpha -e Nemertea -e Platyhelminthes -e Porifera' to target invertebrate animal phyla  
Line 163, set to removal method 1 (uses ORFfinder)  
Line 177, set to '5' to use the invertebrate mitochondrial genetic code for translation  
Run snakemake.  Move invertebrate outfiles into their own directory so they do not get over-written: taxon.zotus, chimera.denoised.nonchimeras.taxon, orf.fasta*, rdp.csv.tmp, results.csv .  

2. Edit the config_ESV.yaml file as follows:  
Line 156, set the taxonomy filter to '-e Chordata' to target animals with a notochord (includes the Vertebrata clade) (see NCBI taxonomy)  
Lin 163, keep removal method 1  
Line 177, set to '2' to use the vertebrate mitochondrial genetic code for translation  
Run snakemake.  Move chordata outfiles into their own directory so they do not get over-written: taxon.zotus, chimera.denoised.nonchimeras.taxon, orf.fasta*, rdp.csv.tmp, results.csv .  

The invertebrate and chordata results.csv files can then be combined prior to downstream processing.  

### The final output file isn't created

If you get a warning and see that the last file created was rdp.csv.tmp but not the expected results.csv then you probably are probably trying to process a very large number of sequence files (thousands) and the job ran into memory problems preventing the creation of the final results.csv file.  You have a few options here:
a) Re-run the pipeline with more memory so that the final outfile can be created.  On the GPSC, you will need to create an interactive container, activate the conda environbment, then unlock the snakemake directory first before submitting another job with more memory on a large compute node.
b) Stop here.  You can still use the ESV.table file together with the rdp.csv.tmp file for further processing in R.  Don't forget to filter out sequence clusters with only 1 or 2 reads from the ESV.table (a step that was done automatically by MetaWorks).  For each sequence cluster in the ESV.table, you can grab the taxonomy from the rdp.csv.tmp file using the ZotuID (sequence cluster ID).
c) There is a new MetaWorks workflow being created that handles analyses with very large number of samples more efficiently (in progress).

## How to cite

If you use this dataflow or any of the provided scripts, please cite the MetaWorks preprint:  
Porter, T.M., Hajibabaei, M. 2020.  METAWORKS: A flexible, scalable bioinformatic pipeline for multi-marker biodiversity assessments.  BioRxiv, doi: https://doi.org/10.1101/2020.07.14.202960.

You can also site this repository:
Teresita M. Porter. (2020, June 25). MetaWorks: A Multi-Marker Metabarcode Pipeline (Version v1.9.7). Zenodo. http://doi.org/10.5281/zenodo.4741407 

If you use this dataflow for making COI taxonomic assignments, please cite the COI classifier publication:  
Porter, T. M., & Hajibabaei, M. (2018). Automated high throughput animal CO1 metabarcode classification. Scientific Reports, 8, 4226.  

If you use the pseudogene filtering methods, please cite the pseudogene publication:
Porter, T.M., & Hajibabaei, M. (2021). Profile hidden Markov model sequence analysis can help remove putative pseudogenes from DNA barcoding and metabarcoding datasets. BMC Bioinformatics, 22: 256. 

If you use the RDP classifier, please cite the publication:  
Wang, Q., Garrity, G. M., Tiedje, J. M., & Cole, J. R. (2007). Naive Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Applied and Environmental Microbiology, 73(16), 5261–5267. doi:10.1128/AEM.00062-07  

## References

Anaconda (2016).  Anaconda Software Distribution.  Available: https://anaconda.com 

Banchi, E.; Ametrano, C.G.; Greco, S.; Stanković, D.; Muggia, L.; Pallavicini, A. PLANiTS: a curated sequence reference dataset for plant ITS DNA metabarcoding. Database 2020, 2020.

Bengtsson-Palme J, Ryberg M, Hartmann M, Branco S, Wang Z, Godhe A, et al. (2013) Improved software detection and extraction of ITS1 and ITS2 from ribosomal ITS sequences of fungi and other eukaryotes for analysis of environmental sequencing data.  Methods in Ecology and Evolution, 4: 914–919. doi:10.1111/2041-210X.12073  

Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. Nature Methods, 13(7), 581–583. doi: 10.1038/nmeth.3869

Edgar, R. C. (2016). UNOISE2: improved error-correction for Illumina 16S and ITS amplicon sequencing. BioRxiv. doi:10.1101/081257  

Gergely, S. (2018, January). Perl-rename. Retrieved from https://github.com/subogero/rename  

Gibson, J., Shokralla, S., Porter, T. M., King, I., Konynenburg, S. van, Janzen, D. H., … Hajibabaei, M. (2014). Simultaneous assessment of the macrobiome and microbiome in a bulk sample of tropical arthropods through DNA metasystematics. Proceedings of the National Academy of Sciences, 111(22), 8007–8012. doi: 10.1073/pnas.1406468111

Hajibabaei, M., Spall, J. L., Shokralla, S., & van Konynenburg, S. (2012). Assessing biodiversity of a freshwater benthic macroinvertebrate community through non-destructive environmental barcoding of DNA from preservative ethanol. BMC Ecology, 12, 28. doi: 10.11.8.0472-6785-12-28

Iwasaki W, Fukunaga T, Isagozawa R, Yamada K, Maeda Y, Satoh TP, et al.( 2013) MitoFish and MitoAnnotator: A Mitochondrial Genome Database of Fish with an Accurate and Automatic Annotation Pipeline. Molecular Biology and Evolution, 30:2531–40. 

Koster, J., Rahmann, S. (2012) Snakemake - a scalable bioinformatics workflow engine.  Bioinformatics, 29(19): 2520-2522. 

Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet. Journal, 17(1), pp–10.  

Porter, T. M., & Hajibabaei, M. (2018). Automated high throughput animal CO1 metabarcode classification. Scientific Reports, 8, 4226.  

Porter, T.M., & Hajibabaei, M. (2021). Profile hidden Markov model sequence analysis can help remove putative pseudogenes from DNA barcoding and metabarcoding datasets. BMC Bioinformatics, 22: 256. 

Pruesse E, Quast C, Knittel K, Fuchs BM, Ludwig W, Peplies J, et al.(2007) SILVA: a comprehensive online resource for quality checked and aligned ribosomal RNA sequence data compatible with ARB. Nucleic Acids Research, 35:7188–96.	

Rimet F, Chaumeil P, Keck F, Kermarrec L, Vasselon V, Kahlert M, et al. (2015) R-Syst::diatom: a barcode database for diatoms and freshwater biomonitoring - data sources and curation procedure. INRA Report. 

Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahé, F. (2016). VSEARCH: a versatile open source tool for metagenomics. PeerJ, 4, e2584. doi:10.7717/peerj.2584  

St. John, J. (2016, Downloaded). SeqPrep. Retrieved from https://github.com/jstjohn/SeqPrep/releases  

Wang, Q., Garrity, G. M., Tiedje, J. M., & Cole, J. R. (2007). Naive Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Applied and Environmental Microbiology, 73(16), 5261–5267. doi:10.1128/AEM.00062-07  

Last updated: April 18, 2022
