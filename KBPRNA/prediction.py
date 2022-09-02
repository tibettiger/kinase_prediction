import os
from xml.dom import INVALID_STATE_ERR
import data
import pandas as pd
import pickle
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.decomposition import PCA 
from sklearn.preprocessing import StandardScaler
import openpyxl
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn import preprocessing
from sklearn.model_selection  import train_test_split
from sklearn.metrics import r2_score
from sklearn import metrics
import xgboost as xgb
from scipy.stats import pearsonr

def xgboost_process(RNAseq, kinaseact, kinasename):
    R2=[]
    best_params=[]
    for i in range(len(kinaseact)):
        cv=KFold(n_splits=5,shuffle=True,random_state=100)
        xg_reg = xgb.XGBRegressor(nthread=4)
        grid_params=dict(max_depth=[4,5,6,7],learning_rate=np.linspace(0.03,0.3,10),n_estimators=[100,200])
        grid_model=GridSearchCV(xg_reg,grid_params,cv=cv,scoring='r2')
        grid_model.fit(RNAseq,kinaseact[i])
        best_estimator=grid_model.best_estimator_
        print(grid_model.best_params_)
        best_params.append(grid_model.best_params_)
        R2.append(grid_model.best_score_)
        s = pickle.dumps(grid_model)
        with open('%d.model' % kinasename[i], 'wb+') as f:
            f.write(s)
    output = pd.concat(R2, best_params)
    output.to_csv("result.csv")


def LR_process(RNAseq, kinaseact, kinasename):
    R2=[]
    for i in range(len(kinaseact)):
        X_train, X_test, y_train, y_test = train_test_split(RNAseq, kinaseact.iloc[i,:], test_size=0.3, random_state=10)
        min_max_scaler = preprocessing.MinMaxScaler()
        X_train = min_max_scaler.fit_transform(X_train)
        y_train = min_max_scaler.fit_transform(y_train.reshape(-1,1))
        X_test = min_max_scaler.fit_transform(X_test)
        y_test = min_max_scaler.fit_transform(y_test.reshape(-1,1))
        lr = LinearRegression()
        lr.fit(X_train, y_train)
        y_test_pred = lr.predict(X_test)
        R2.append(r2_score(y_test, y_test_pred))
        with open('%d.model' % kinasename[i], 'wb+') as f:
            f.write(s)
    R2.to_csv("result.csv")



class KBPRNA:
    """KBPRNA class
    
    
    KBPRNA holds the main modelling functionalites. Its basic usage
    includes:
    1. load data and preprocessing
    2. initialize a model
    3. fit the model
    4. make the predictions of kinase activities
    
    
    Paramters
    --------------
    out_dir : str
        Relative path or name of the desired output directory
    core_model : str
        Core model algorithm. Supported: "Xgboost",  "linear regression"
    in_dir : str
        Relative path or name of the desired input directory
    """

    def __init__(
        self,
        out_dir,
        core_model,
        in_dir
    ):
        self.results = pd.Dataframe()
        self.out_dir = out_dir
        self.core_model = core_model
        self.in_str = in_dir
        self.fitted = False
        self.predicted = False
    
    def fit(self,
            kinase_list,
            cancer_type,
            select_all,
            core_model,
            save = False):
            """Method to fit the model.
            
            Use provided datasets, performing training and fit KBPRNA core_model.
            
            Parameters
            -------------
            kinase_list : list of kinases targeted for prediction
            cancer_type : choose which model of cancer type to be built
            select_all : choose if L1000 genes or total genes are considered
            core_model : we provide two model: xgboost and linear regression for regression prediction
            save : bool
                whether or not to save the fitted core model in KBPRNA out_dir.
                Default: False"""
            
            kinase_data = data.loadfeatset("%d_kinase.csv" % cancer_type)
            kinase_act = []
            for i in range(len(kinase_list)):
                for j in range(len(kinase_data)):
                    if(kinase_list[i,0]==kinase_data[j,0]):
                        kinase_act.append(kinase_data[i,:])
            if not select_all:
                RNAseq = data.preprocess_data_L1000(cancer_type)
            else:
                RNAseq = data.preprocess_data_all(cancer_type)

            if(core_model=="XGboost"):
                xgboost_process(RNAseq,kinase_act)
            if(core_model=="linear_regression"):
                LR_process(RNAseq,kinase_act)
            



    def predict(
        self, file, kinase_name
    ):
        """Method to make predictions.
        
        Make predictions for unknown datasets
        we used pretrained model to predict kinase activities"""
        predict_result = []
        for i in range(len(kinase_name)):
            f = open('%d.model' % kinase_name[i], 'rb')
            s = f.read()
            model = pickle.loads(s)
            data = data.loadfeatset(file)
            predict_result.append(model.predict(data))
        predict_result = pd.concat(kinase_name, predict_result)
        predict_result.to_csv(predict_result)

