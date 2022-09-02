import os
import scipy
import pandas as pd
from sklearn.decomposition import PCA 
from sklearn.preprocessing import StandardScaler
import openpyxl
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import ShuffleSplit


def loadfeatset(filename):
    numFeat=len(open(filename).readline().split(','))
    dataMat=[]
    fr=open(filename)
    for line in fr.readlines():
        lineArr=[]
        curline=line.strip().split(',')
        for i in range(numFeat):
            lineArr.append(curline[i])
        dataMat.append(lineArr)
    return dataMat

def preprocess_data_all(cancer_type):
    data = loadfeatset("%d_mRNA.csv" % cancer_type)
    data = pd.Dataframe(data)
    newdatamat = PCA(n_components=20).fit_transform(data)
    return newdatamat

def preprocess_data_L1000(cancer_type):
    """Helper function to transform a original bulk RNA-seq dataset to a dataMat containing 
    information with LINCS-L1000 genes"""
    L1000 = pd.read_table("L1000.txt", sep = "\t")
    data = loadfeatset("%d_mRNA.csv" % cancer_type)
    data = pd.Dataframe(data)
    newmat = []
    for i in len(L1000):
        for j in len(data):
            if data.iloc[j,1] == L1000.iloc[i,1]:
                newmat.append(data.iloc[j,:])
    newdatamat=PCA(n_components=20).fit_transform(newmat)
    return newdatamat   



