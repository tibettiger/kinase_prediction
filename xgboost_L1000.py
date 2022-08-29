
##xgboost  (L1000 feature) ABL1 prediction

from cProfile import label
from cgi import test
from random import random
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

from sklearn.model_selection  import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn import metrics
import xgboost as xgb
from scipy.stats import pearsonr

def loadfeatset(filename):
    numFeat=len(open(filename).readline().split('\t'))
    dataMat=[]
    fr=open(filename)
    for line in fr.readlines():
        lineArr=[]
        curline=line.strip().split('\t')
        for i in range(numFeat):
            lineArr.append(curline[i])
        dataMat.append(lineArr)
    return dataMat

dataMat=loadfeatset('E:/激酶活性预测课题/L1000_HCC.txt')
kinase=loadfeatset('E:/激酶活性预测课题/HCC_kinase.txt')
dataMat=pd.DataFrame(dataMat)
kinase=pd.DataFrame(kinase)

##PCA降维
from sklearn.decomposition import PCA 
from sklearn.preprocessing import StandardScaler
import openpyxl
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import ShuffleSplit
sc=StandardScaler()
newdataMat=PCA(n_components=20).fit_transform(dataMat)
kinase=sc.fit_transform(kinase)
print(newdataMat.shape[0])
print(newdataMat.shape[1])
print(kinase.shape[0])
print(kinase.shape[1])
R2=[]
best_params=[]
for i in range(185):
    cv=KFold(n_splits=5,shuffle=True,random_state=100)
    xg_reg = xgb.XGBRegressor(nthread=4)
    grid_params=dict(max_depth=[4,5,6,7],learning_rate=np.linspace(0.03,0.3,10),
    n_estimators=[100,200])
    grid_model=GridSearchCV(xg_reg,grid_params,cv=cv,scoring='r2')
    grid_model.fit(newdataMat,kinase[i])
    best_estimator=grid_model.best_estimator_
    print(grid_model.best_params_)
    best_params.append(grid_model.best_params_)
    R2.append(grid_model.best_score_)
    print(grid_model.best_score_)

R2=pd.DataFrame(R2)
R2.to_excel("E:/R2.xlsx")  
best_params=pd.DataFrame(best_params)
best_params.to_excel("E:/params.xlsx")


from xgboost import plot_importance
plot_importance(xg_reg)




##xgboost  (substrate gene)

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

from sklearn.model_selection  import train_test_split
from sklearn.metrics import mean_squared_error

import xgboost as xgb


dataMat=loadfeatset('E:/激酶活性预测课题/L1000_HCC.txt')
kinase=loadfeatset('E:/激酶活性预测课题/HCC_kinase.txt')
dataMat=pd.DataFrame(dataMat)
kinase=pd.DataFrame(kinase)


##用底物基因的表达量预测激酶活性
def get_substrate_gene(filename,kinase_name):
    import csv
    substrate=[]
    with open(filename,'r') as f:
        reader=csv.reader(f)
        for row in reader:
            if(row[0]==kinase_name ):
                substrate.append(row[2])
    return substrate





def build_expressionmatrix(substrate):
    expmat=[]
    import csv
    with open('E:/激酶活性预测课题/HCC_mRNA.csv','r') as f:
        reader=csv.reader(f)
        for row in reader:
            for i in range(len(substrate)):
                if(row[0]==substrate[i]):
                    expmat.append(row)
    return expmat


##xgboost substrate gene
kinase_name=[]

with open('E:/激酶活性预测课题/Cell_HCC_KSEA.csv','r') as f:
    import csv
    reader=csv.reader(f)
    for row in reader:
        kinase_name.append(row[1])
R2=[]
best_params=[]
kinase_name=kinase_name[1:]
print(len(kinase))
for i in range(len(kinase)):
    substrate=get_substrate_gene('E:/激酶活性预测课题/kinase_substrate.csv',kinase_name[i])
    expmat=build_expressionmatrix(substrate)
    expmat=pd.DataFrame(expmat)
    if expmat.empty:
        R2.append('NA')
        best_params.append('NA')
    else:
        expmat.drop(columns=[0],inplace=True)
        expmat_transpose=pd.DataFrame(expmat.values.T,index=expmat.columns,columns=expmat.index)
        expmat_transpose=sc.fit_transform(expmat_transpose)
        cv=KFold(n_splits=5,shuffle=True,random_state=100)
        xg_reg = xgb.XGBRegressor(nthread=4)
        grid_params=dict(max_depth=[4,5,6,7],learning_rate=np.linspace(0.03,0.3,10),
        n_estimators=[100,200])
        grid_model=GridSearchCV(xg_reg,grid_params,cv=cv,scoring='r2')
        grid_model.fit(expmat_transpose,kinase[i,:])
        best_estimator=grid_model.best_estimator_
        print(grid_model.best_params_)
        best_params.append(grid_model.best_params_)
        R2.append(grid_model.best_score_)
        print(grid_model.best_score_)

R2=pd.DataFrame(R2)
R2.to_excel("E:/R2.xlsx")  
best_params=pd.DataFrame(best_params)
best_params.to_excel("E:/params.xlsx")


from xgboost import plot_importance
plot_importance(xg_reg)
print(len(substrate))
print(len(kinase_name))
print(len(kinase.iloc[1,:]))
print(len(R2))
print(len(best_params))






##先做特征筛选，再进行prediction
from cProfile import label
from cgi import test
from random import random
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

from sklearn.model_selection  import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn import metrics
import xgboost as xgb
from scipy.stats import pearsonr
from sklearn.decomposition import PCA 
from sklearn.preprocessing import StandardScaler
import openpyxl
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import ShuffleSplit

#先算每个基因与激酶的相关系数，再选取相关系数最大的前二十个激酶作为输入特征

R2=[]
best_params=[]
for i in range(185):
    cv=KFold(n_splits=5,shuffle=True,random_state=100)
    xg_reg = xgb.XGBRegressor(nthread=4)
    grid_params=dict(max_depth=[4,5,6,7],learning_rate=np.linspace(0.03,0.3,10),
    n_estimators=[100,200])
    grid_model=GridSearchCV(xg_reg,grid_params,cv=cv,scoring='r2')

    grid_model.fit(newdataMat,kinase[i])
    best_estimator=grid_model.best_estimator_
    print(grid_model.best_params_)
    best_params.append(grid_model.best_params_)
    R2.append(grid_model.best_score_)
    print(grid_model.best_score_)

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

##用相关系数前20基因的表达量预测激酶活性+线性回归
expmat=[]
import csv
with open('C:/学习/激酶预测课题/kinase_dataset/BC.csv','r') as f:
    reader=csv.reader(f)
    for row in reader:
        expmat.append(row)
expmat=pd.DataFrame(expmat)
expmat=expmat.drop(0,axis=1)
expmat=expmat.drop(1,axis=1)

def get_top_twenty_gene(kin_num):
    import csv
    pearson=[]
    for i in range(len(expmat)):
        exp=np.array(expmat.iloc[i])
        exp=exp.tolist()
        kin=np.array(kinase.iloc[kin_num])
        kin=kin.tolist()
        for j in range(len(exp)):
            exp[j]=float(exp[j])
            kin[j]=float(kin[j])
        from scipy.stats import pearsonr
        pccs=pearsonr(exp,kin)
        pearson.append(pccs[0])
    import heapq
    re=map(pearson.index,heapq.nlargest(20,pearson))
    return list(re)
R2=[]
best_params=[]
sc=StandardScaler()
kinase=sc.fit_transform(kinase)
import csv
for i in range(190):
    kin_num=get_top_twenty_gene(i)
    newdatamat=[]
    for j in range(len(kin_num)):
        newdatamat.append(expmat.iloc[kin_num[j]])
    newdatamat=pd.DataFrame(newdatamat)
    expmat_transpose=pd.DataFrame(newdatamat.values.T,index=newdatamat.columns,columns=newdatamat.index)
    expmat_transpose=sc.fit_transform(expmat_transpose)
    cv=KFold(n_splits=5,shuffle=True,random_state=100)
    xg_reg = xgb.XGBRegressor(nthread=4)
    grid_params=dict(max_depth=[4,5,6,7],learning_rate=np.linspace(0.03,0.3,10),n_estimators=[100,200])
    grid_model=GridSearchCV(xg_reg,grid_params,cv=cv,scoring='r2')
    grid_model.fit(expmat_transpose,kinase.iloc[i])
    best_estimator=grid_model.best_estimator_
    print(grid_model.best_params_)
    best_params.append(grid_model.best_params_)
    R2.append(grid_model.best_score_)
    print(grid_model.best_score_)
    

kinase=loadfeatset('C:/学习/激酶预测课题/Cell_BC_KA_2.csv')
import pandas as pd
kinase=pd.DataFrame(kinase)
kinase=kinase.drop(0,axis=1)
kinase=kinase.drop(0,axis=0)
print(len(kinase.iloc[2]))
R2=pd.DataFrame(R2)
R2.to_excel("C:/学习/激酶预测课题/R2_BC_top20gene.xlsx")  
best_params=pd.DataFrame(best_params)
best_params.to_excel("C:/学习/激酶预测课题/params_BC_top20gene.xlsx")
import openpyxl
import numpy as np
print(newdatamat.shape[0])
print(kinase.iloc[0])

