# MetaWorks: A Multi-Marker Metabarcode Pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4741407.svg)](https://doi.org/10.5281/zenodo.4741407)  

MetaWorks generates exact sequence variants and/or operational taxonomic units and taxonomically assigns them.  Supports a number of popular metabarcoding markers: COI, rbcL, ITS, SSU rRNA, and 12S SSU mtDNA.  See the [MetaWorks website](https://terrimporter.github.io/MetaWorksSite) for quickstart guides, additional pipeline details, FAQs, and a step-by-step tutorial that includes installation.

## Installation

MetaWorks runs at the command-line on linux x86-64 in a conda environment (provided).

Instructions for installing [conda](https://terrimporter.github.io/MetaWorksSite/tutorial/#InstallConda) (if not already installed).

Instructions for installing [ORFfinder](https://terrimporter.github.io/MetaWorksSite/tutorial/#InstallORFfinder) if pseudogene-filtering will be run (optional).

Instructions for installing [MetaWorks](https://terrimporter.github.io/MetaWorksSite/tutorial/#InstallMetaWorks) and activating the [MetaWorks conda environment](https://terrimporter.github.io/MetaWorksSite/tutorial/#ActivateMetaWorksEnv).

Instructions on where to find [custom-trained classifiers](https://terrimporter.github.io/MetaWorksSite/#classifier_table) that can be used with MetaWorks.

## Documentation

A [quickstart](https://terrimporter.github.io/MetaWorksSite/quickstart/) guide to various workflows.

A [detailed](https://terrimporter.github.io/MetaWorksSite/details/) explanation of MetaWorks workflows.

A [tutorial](https://terrimporter.github.io/MetaWorksSite/tutorial/) provides step-by-step instructions on how to prepare your environment and get started quickly using the provided test data.

*NEW* We added answers to some frequently asked questions (FAQs) about MetaWorks and data analysis to the [MetaWorks website](https://terrimporter.github.io/MetaWorksSite/faq).

## How to cite

If you use this dataflow or any of the provided scripts, please cite the MetaWorks preprint:  
Porter, T.M., Hajibabaei, M. 2020.  METAWORKS: A flexible, scalable bioinformatic pipeline for multi-marker biodiversity assessments.  BioRxiv, doi: https://doi.org/10.1101/2020.07.14.202960.

You can also site this repository:
Teresita M. Porter. (2020, June 25). MetaWorks: A Multi-Marker Metabarcode Pipeline (Version v1.10.0). Zenodo. http://doi.org/10.5281/zenodo.4741407 

If you use this dataflow for making COI taxonomic assignments, please cite the COI classifier publication:  
Porter, T. M., & Hajibabaei, M. (2018). Automated high throughput animal CO1 metabarcode classification. Scientific Reports, 8, 4226.  

If you use the pseudogene filtering methods, please cite the pseudogene publication:
Porter, T.M., & Hajibabaei, M. (2021). Profile hidden Markov model sequence analysis can help remove putative pseudogenes from DNA barcoding and metabarcoding datasets. BMC Bioinformatics, 22: 256. 

If you use the RDP classifier, please cite the publication:  
Wang, Q., Garrity, G. M., Tiedje, J. M., & Cole, J. R. (2007). Naive Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Applied and Environmental Microbiology, 73(16), 5261–5267. doi:10.1128/AEM.00062-07  

Last updated: July 30, 2022
