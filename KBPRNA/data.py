import os
import scipy
import pandas as pd
from sklearn.decomposition import PCA 
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import ShuffleSplit
import numpy as np

def preprocess_data_all(cancer_type):
    dataMat = pd.read_csv("RNA/%s_RNA.csv" % cancer_type, index_col = 0)
    newdatamat = PCA(n_components=20).fit_transform(dataMat)
    return newdatamat

def preprocess_data_L1000(cancer_type):
    """Helper function to transform a original bulk RNA-seq dataset to a dataMat containing 
    information with LINCS-L1000 genes"""
    L1000 = pd.read_table("L1000.txt")
    dataMat = pd.read_csv("RNA/%s_RNA.csv" % cancer_type, index_col=0)
    newdata = pd.DataFrame()
    for i in range(L1000.shape[0]):
        if L1000.iloc[i].values in dataMat.index.values:
            newdata = newdata.append(pd.DataFrame(dataMat.loc[L1000.iloc[i],:]))
    newdatamat=PCA(n_components=20).fit_transform(newdata)
    return newdatamat   


