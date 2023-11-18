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
    ###obtain high variable genes
    dataMat = pd.read_csv("../RNA/%s_RNA.csv" % cancer_type, index_col = 0)
    #se = dataMat.var(axis=1)
    #feature = se.sort_values(ascending=False).index[0:50].tolist()
    #newdatamat = dataMat.loc[feature,:]
    newdatamat = dataMat
    outdata = (newdatamat - newdatamat.min())/(newdatamat.max() - newdatamat.min())
    return outdata.T

def preprocess_data_L1000(cancer_type):
    """Helper function to transform a original bulk RNA-seq dataset to a dataMat containing 
    information with LINCS-L1000 genes"""
    L1000 = pd.read_table("../L1000.txt")
    dataMat = pd.read_csv("../RNA/%s_RNA.csv" % cancer_type, index_col=0)
    ncol=dataMat.shape[0]
    nrow=dataMat.shape[1]
    newdata = pd.DataFrame()
    for i in range(L1000.shape[0]):
        if L1000.iloc[i].values in dataMat.index.values:
            newdata = pd.concat([newdata,dataMat.loc[L1000.iloc[i], :]],axis=0)
        else:
            genename=''.join(L1000.iloc[i].values)
            newdata.loc[genename,:]=[0.5]*newdata.shape[1]
    outdata = (newdata - newdata.min())/(newdata.max() - newdata.min())
    #se = newdata.var(axis=1)
    #feature = se.sort_values(ascending=False).index[0:50].tolist()
    #newdatamat = newdata.loc[feature,:]
    return outdata.T   
