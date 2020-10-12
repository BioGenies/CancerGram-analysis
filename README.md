# Read me

This repository contains data and code necessary to reproduce analysis from the article: Burdukiewicz, M., Sidorczuk, K., Rafacz, D., Pietluch, F., Bakala, M., Slowik, J. and Gagat, P. CancerGram: an effective classifier for differentiating anticancer from antimicrobial peptides.

The analysis conducted in this article resulted in a predictor of anticancer peptides CancerGram, available as a *R* package and a web server (http://biongram.biotech.uni.wroc.pl/CancerGram/).


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

* **mito_ACPs.fasta** - sequences of 11 mitochondrial ACPs not included in the training datasets and used to evaluate CancerGram performance.
* **mitochondrial_ACPs_table.csv** - table containing names and sequences of 11 mitochondrial ACPs mentioned above with their sources.

Data used in the exploratory analysis:

* **DRAMP_Anticancer_amps.csv** - data on anticancer AMPs downloaded from the [DRAMP database](http://dramp.cpu-bioinfor.org/).
* **DRAMP_Antitumor_amps.csv*** - data on antitumor AMPs downloaded from the [DRAMP database](http://dramp.cpu-bioinfor.org/).
* **all_ACPs.fasta** - sequences of ACPs combined from different sources (DRAMP, APD, and CancerPPD databases).
* **apd_df.csv** - data downloaded from the APD database.
* **cancerppd_l_natural.txt** - sequences of anticancer peptides downloaded from the [CancerPPD database](http://crdd.osdd.net/raghava/cancerppd/downseq.php).
* **dbamp_df.csv** - data downloaded from the [dbAMP database](http://140.138.77.240/~dbamp/).
* **pos_test_alternate.txt** - ACP sequences from the alternative validation dataset, acquired from [Agrawal et al.](https://doi.org/10.1093/bib/bbaa153)
* **pos_train_alternate.txt** - ACP sequences from the alternative training dataset, acquired from [Agrawal et al.](https://doi.org/10.1093/bib/bbaa153)


### functions

All functions necessary to repeat the analysis.

### reports

Report summing up the results obtained with the first models. 

### results

* **CancerGram_model.rda** - object containing CancerGram stacked random forest model and important features
* **benchmark.fasta** - sequences used as a benchmark dataset in the exploratory analysis.
* **benchmark_mc.fasta** - sequences used as a benchmark dataset for multiclass model in the exploratory analysis.

### data_exploration

Scripts and functions used during the exploratory analysis, mostly undocumented.


## RSession information

All scripts used in this study are compatible with following version of R:

* **R version 3.6.2 (2019-12-12)**
* **Platform:** x86_64-pc-linux-gnu (64-bit)

Necessary packages and their versions used in the analyses are listed in **renv.lock**

