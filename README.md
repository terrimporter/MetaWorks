# MetaWorks: A Multi-Marker Metabarcode Pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4741407.svg)](https://doi.org/10.5281/zenodo.4741407)  

MetaWorks is a multi-marker metabarcode pipeline.  Produces  exact sequence variants and/or operational taxonomic units.  Supports a number of popular metabarcoding markers: COI, rbcL, ITS, SSU rRNA, and 12S SSU mtDNA.  See the [MetaWorks website](https://terrimporter.github.io/MetaWorksSite) for installation instructions, quickstart guides, FAQs, and tutorial.

## Installation

MetaWorks runs at the command-line on linux x86-64 in a conda environment (provided).

Instructions for installing [conda](https://terrimporter.github.io/MetaWorksSite/faq/) (if not already installed).

Instructions for installling [ORFfinder](https://terrimporter.github.io/MetaWorksSite/workflows/) if pseudogene-filtering will be run (optional).

Instructions for loading the [MetaWorks conda environment](https://terrimporter.github.io/MetaWorksSite/workflows/).

## Documentation

A [quickstart](https://terrimporter.github.io/MetaWorksSite/quickstart/) guide to various workflows.

A [detailed](https://terrimporter.github.io/MetaWorksSite/workflows/) explanation of MetaWorks workflows.

A [tutorial](https://terrimporter.github.io/MetaWorksSite/tutorial/) using provided test data to get up and running quickly.

Additional documentation and FAQs are available on the [MetaWorks website](https://terrimporter.github.io/MetaWorksSite/).

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

## Papers & projects that use MetaWorks

Edge TA, Baird DJ, Bilodeau G, Gagné N, Greer C, Konkin D, et al. 2020.  The Ecobiomics project: Advancing metagenomics assessment of soil health and freshwater quality in Canada. Science of The Total Environment, 710: 135906. doi:10.1016/j.scitotenv.2019.135906![image](https://user-images.githubusercontent.com/34171056/181351437-084aa653-91e6-4071-941b-7b211d89429a.png)

Moir, C. 2021. No Stomach, No Problem: an Integrated Morpho-Molecular Approach to Assessing the Diets of the Cunner Wrasse, Tautogolabrus adspersus, among Coastal, Nearshore Regions of Atlantic Canada (Doctoral dissertation, University of Guelph).

Porter, TM & Hajibabaei, M, 2020. METAWORKS: A flexible, scalable bioinformatic pipeline for multi-marker biodiversity assessments. bioRxiv.

Porter, TM, & Hajibabaei, M. (2021). Profile hidden Markov model sequence analysis can help remove putative pseudogenes from DNA barcoding and metabarcoding datasets. BMC Bioinformatics, 22(1): 256. doi: 10.1186/s12859-021-04180-x

Robinson, CV, Baird, DJ, Wright, MTG, Porter, TM, Hartwig, K, Hendriks, E, Maclean, L, Mallinson, R, Monk, WA, Paquette, C and Hajibabaei, M. 2021. Combining DNA and people power for healthy rivers: Implementing the STREAM community-based approach for global freshwater monitoring. Perspectives in Ecology and Conservation, 19(3): 279-285.

Robinson, CV, Porter, TM, Maitland, VC, Wright, MT and Hajibabaei, M. 2021. Multi-marker metabarcoding resolves subtle variations in freshwater condition: Bioindicators, ecological traits, and trophic interactions. bioRxiv.

Rudar, J, Golding, GB, Kremer, SC and Hajibabaei, M. 2022. Decision Tree Ensembles Utilizing Multivariate Splits Are Effective at Investigating Beta-Diversity in Medically Relevant 16S Amplicon Sequencing Data. bioRxiv.

Smenderovac E, Emilson C, Porter T, Morris D, Hazlett P, Diochon A, et al. 2022.  Forest soil biotic communities show few responses to wood ash applications at multiple sites across Canada. Sci Rep., 12: 4171. doi:10.1038/s41598-022-07670-x![image](https://user-images.githubusercontent.com/34171056/181351582-1d51af83-636e-4470-8e9f-d3e7cce33755.png)

Last updated: July 27, 2022
