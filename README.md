# KBPRNA

KBPRNA is a novel machine learning method based on bulk RNA-seq data to predict kinase activity with a decent accuracy. This model could input bulk RNA-seq dataset of cancer tissues from patients and output corresponding highly representive kinases' activity. This model could replace traditional kinase activity computing methods using phosphorylated proteomics sequencing. The first step of KBPRNA is to leverage L1000 genes' expression levels from original bulk RNA-seq datasets and then apply normalization methods and PCA to deal with LINCS-L1000 genes' expression levels. XGboost (eXtreme Gradient Boosting) which proves an optimized distributed gradient boosting with high accuracy in predicting continuous variable was utilized as the basic framework. Gridsearchcv was used to find the best parameters. The latter part of this model was to predict kinases' activity of our additional test file to show the model's robustness.

## installation
KBPRNA currently supports python >= 3.8, and can be installed from PyPI
pip install KBPRNA

## usage

You can follow our model's usage method through tutorial.ipynb






