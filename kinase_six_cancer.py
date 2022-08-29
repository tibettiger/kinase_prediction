

###用L1000+PCA降维来预测六种癌症的激酶活性

cancer_name=["BC","GBM","HCC","LSCC","LUAD","UCEC"]


from cProfile import label
from cgi import test
from random import random
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA 
from sklearn.preprocessing import StandardScaler
import openpyxl
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import ShuffleSplit
from sklearn.model_selection  import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn import metrics
import xgboost as xgb
from scipy.stats import pearsonr
##读取激酶数据
def loadfeatset(filename,inter):
    numFeat=len(open(filename).readline().split(inter))
    dataMat=[]
    fr=open(filename)
    for line in fr.readlines():
        lineArr=[]
        curline=line.strip().split(inter)
        for i in range(numFeat):
            lineArr.append(eval(curline[i]))
        dataMat.append(lineArr)
    return dataMat

kinase=loadfeatset("C:/学习/激酶预测课题/kinase_dataset/LSCC.csv",",")
kinase=np.asarray(kinase)
kinase=np.delete(kinase,0,axis=0)
kinase=np.delete(kinase,0,axis=1)  ##repeat twice
print(kinase.shape[1])
#将输入文件转换为L1000格式
L1000=[]
f1=open("C:/学习/激酶预测课题/L1000.txt","r")
line=f1.readlines()
for i in line:
    L1000.append(eval(i))
    print(eval(i))

RNAseq=loadfeatset("C:/学习/港中深博一/激酶预测课题/RNA-seq_dataset/LSCC.csv",",")
L1000_RNAseq=[]

import numpy as np
RNAseq=np.asarray(RNAseq)
print(len(RNAseq))
for i in range(len(RNAseq)):
    for j in range(len(L1000)):
        if(RNAseq[i,0]==L1000[j]):
            L1000_RNAseq.append(RNAseq[i])
L1000_RNAseq=np.asarray(L1000_RNAseq)
L1000_RNAseq=np.delete(L1000_RNAseq,0,axis=1)
print(L1000_RNAseq.shape[1])
##构建xgboost模型+PCA降维
sc=StandardScaler()
L1000_RNAseq=PCA(n_components=20).fit_transform(np.transpose(L1000_RNAseq))
print(L1000_RNAseq.shape[1])
kinase=sc.fit_transform(kinase)
print(kinase)
R2=[]
best_params=[]
for i in range(len(kinase)):
    cv=KFold(n_splits=5,shuffle=True,random_state=100)
    xg_reg = xgb.XGBRegressor(nthread=4)
    grid_params=dict(max_depth=[4,5,6],learning_rate=np.linspace(0.03,0.27,9))
    grid_model=GridSearchCV(xg_reg,grid_params,cv=cv,scoring='r2')
    grid_model.fit(L1000_RNAseq,kinase[i])
    best_estimator=grid_model.best_estimator_
    print(grid_model.best_params_)
    best_params.append(grid_model.best_params_)
    R2.append(grid_model.best_score_)
    print(grid_model.best_score_)

R2=pd.DataFrame(R2)
R2.to_excel("C:/学习/激酶预测课题/result/R2_LSCC.xlsx")  
best_params=pd.DataFrame(best_params)
best_params.to_excel("C:/学习/激酶预测课题/result/params_LSCC.xlsx")

##激酶底物作为输入特征





##前二十高相关性基因作为输入特征


