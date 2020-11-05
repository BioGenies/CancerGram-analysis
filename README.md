[![published in: Pharmaceutics](https://img.shields.io/badge/published%20in-Pharmaceutics-green.svg)](https://doi.org/10.3390/pharmaceutics12111045)

This repository contains data and code necessary to reproduce analysis from the article: Burdukiewicz, M., Sidorczuk, K., Rafacz, D., Pietluch, F., Bąkała, M., Słowik, J., and Gagat, P. (2020). CancerGram: An Effective Classifier for Differentiating Anticancer from Antimicrobial Peptides. Pharmaceutics 12, 1045, https://doi.org/10.3390/pharmaceutics12111045.

The analysis conducted in this article resulted in a predictor of anticancer peptides CancerGram, available as a [*R* package](https://github.com/BioGenies/CancerGram) and a web server (http://biongram.biotech.uni.wroc.pl/CancerGram/).

## How to reproduce the main part of the analysis?

Source **analysis.R**. Computations are time consuming.

## How to generate results and plots for publication?

Run **publication_results.R**. This script requires data generated in the main analysis.

## Repository structure

### data

Data used in the main analysis: 

* **pos_test_main.txt** - ACP sequences for validation dataset, acquired from [Agrawal et al.](https://doi.org/10.1093/bib/bbaa153)
* **neg_test_alternate.txt** - Negative sequences for validation dataset, acquired from [Agrawal et al.](https://doi.org/10.1093/bib/bbaa153)
* **neg_test_main.txt** - AMP sequences for validation dataset, acquired from [Agrawal et al.](https://doi.org/10.1093/bib/bbaa153)
* **pos_train_main.txt** - ACP sequences for training dataset, acquired from [Agrawal et al.](https://doi.org/10.1093/bib/bbaa153)
* **neg_train_alternate.txt** - Negative sequences for training dataset, acquired from [Agrawal et al.](https://doi.org/10.1093/bib/bbaa153)
* **neg_train_main.txt** - AMP sequences for training dataset, acquired from [Agrawal et al.](https://doi.org/10.1093/bib/bbaa153)

Data used for the publication results:

* **mito_ACPs.fasta** - sequences of 12 mitochondrial ACPs not included in the training datasets and used to evaluate CancerGram performance.
* **mitochondrial_ACPs_table.csv** - table containing names and sequences of 12 mitochondrial ACPs mentioned above with their sources.
* **DRAMP_Anticancer_amps.csv** - data on anticancer AMPs downloaded from the [DRAMP database](http://dramp.cpu-bioinfor.org/).
* **DRAMP_Antitumor_amps.csv*** - data on antitumor AMPs downloaded from the [DRAMP database](http://dramp.cpu-bioinfor.org/).
* **all_ACPs.fasta** - sequences of ACPs combined from different sources (DRAMP, APD, and CancerPPD databases).
* **apd_df.csv** - data downloaded from the [APD3 database](http://aps.unmc.edu/AP/main.php).
* **cancerppd_l_natural.txt** - sequences of anticancer peptides downloaded from the [CancerPPD database](http://crdd.osdd.net/raghava/cancerppd/downseq.php).
* **dbamp_df.csv** - data downloaded from the [dbAMP database](http://140.138.77.240/~dbamp/).
* **mACPpred_negative.fasta** - negative dataset used for mACPpred training, acquired from http://thegleelab.org/mACPpred/ACPData.html on 26.10.2020
* **mACPpred_positive.fasta** - positive dataset used for mACPpred training, acquired from http://thegleelab.org/mACPpred/ACPData.html on 26.10.2020
* **pos_test_main_without_mACPpred.fasta** - ACP validation dataset after removal of peptides used for mACPpred training

Data used in the exploratory analysis:

* **pos_test_alternate.txt** - ACP sequences from the alternative validation dataset, acquired from [Agrawal et al.](https://doi.org/10.1093/bib/bbaa153)
* **pos_train_alternate.txt** - ACP sequences from the alternative training dataset, acquired from [Agrawal et al.](https://doi.org/10.1093/bib/bbaa153)


### functions

All functions necessary to repeat the analysis.

### reports

Report summing up the results obtained with the first models. 

### results

* **CancerGram_model.rda** - object containing CancerGram stacked random forest model and important features
* **independent_anticp_local_model1_0.5.csv** - prediction results on the independent dataset obtained with AntiCP 2.0 local version (downloaded from https://github.com/raghavagps/anticp2/ on 26.10.2020) with model 1, SVM cut-off set to 0.5 (default) and window length set to 10 (default)
* **independent_dataset_for_anticp_benchmark.fa** - independent dataset used for benchmarking CancerGram with AntiCP 2.0. Dataset was constructed using experimentally verified ACP sequences downloaded from [DRAMP database](http://dramp.cpu-bioinfor.org/), [CancerPPD database](http://crdd.osdd.net/raghava/cancerppd/downseq.php) and [APD3 database](http://aps.unmc.edu/AP/main.php), as well as experimentally verified AMPs from [dbAMP database](http://140.138.77.240/~dbamp/). To reduce homology of both datasets, we used CD-HIT with identity cut-off 0.95 and 0.6 for ACPs and AMPs respectively. Next, we removed sequences present in the training and validation datasets aquired from [Agrawal et al.](https://doi.org/10.1093/bib/bbaa153) that were used for training CancerGram and/or AntiCP 2.0, obtaining 57 ACPs and 769 AMPs. 
* **mACPpred_predictions_validation.csv** - prediction results on the ACP and AMP validation dataset, from which we removed sequences used for mACPpred training, obtaining 128 ACPs and 170 AMPs. 

### data_exploration

Scripts and functions used during the exploratory analysis, mostly undocumented.


## RSession information

All scripts used in this study are compatible with following version of R:

* **R version 3.6.2 (2019-12-12)**
* **Platform:** x86_64-pc-linux-gnu (64-bit)

Necessary packages and their versions used in the analyses are listed in **renv.lock**

