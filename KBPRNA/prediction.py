from operator import index
import os
import data
from xml.dom import INVALID_STATE_ERR
import data
import pandas as pd
import pickle
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.svm import SVR
from sklearn.decomposition import PCA 
from sklearn.preprocessing import StandardScaler
import openpyxl
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn import preprocessing
from sklearn.model_selection  import train_test_split
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from sklearn import metrics
import xgboost as xgb
from scipy.stats import pearsonr
import re
from sklearn.feature_selection import RFE

def xgboost_process(RNAseq, kinaseact):
    print("XGboost!")
    R2=[]
    RMSE = []
    MAE = []
    best_params=[]
    kinasename = kinaseact.index.values
    print(kinasename)
    kinaseact = kinaseact.values
    for i in range(len(kinasename)):
        X_train, X_test, Y_train, Y_test = train_test_split(RNAseq, kinaseact[i], test_size=0.25, random_state=10)
        cv=KFold(n_splits=5,shuffle=True,random_state=100)
        xg_reg = xgb.XGBRegressor(nthread=4)
        grid_params=dict(max_depth=[2,3,4,5,6,7],n_estimators=[100,200,500])
        grid_model=GridSearchCV(xg_reg,grid_params,cv=cv,scoring='r2')
        grid_model.fit(X_train,Y_train)
        best_estimator=grid_model.best_estimator_
        feature=best_estimator.feature_importances_
        feature_name=RNAseq.columns
        feature_result = list(zip(feature, feature_name))
        print(sorted(feature_result, key=lambda t: t[0],reverse=True)[0:10])
        ##explainer = shap.TreeExplainer(grid_model.best_estimator_,X_test)
        ##shap_values = explainer.shap_values(X_test)
        ##shap.summary_plot(shap_values, X_test)
        y_pred = best_estimator.predict(X_test)
        print(y_pred)
        ##computing R2
        pred_R2 = r2_score(Y_test,y_pred)
        print(pred_R2)
        R2.append(pred_R2)
        ##computing RMSE
        pred_RMSE = np.sqrt(mean_squared_error(Y_test,y_pred))
        RMSE.append(pred_RMSE)
        ##computing MAE
        pred_MAE = mean_absolute_error(Y_test,y_pred)
        MAE.append(pred_MAE)
        print(grid_model.best_params_)
        best_params.append(grid_model.best_params_)
        s = pickle.dumps(grid_model)
        with open('%s.xgb.model' % kinasename[i], 'wb+') as f:
            f.write(s)
        
    R2 = pd.DataFrame(R2,index = kinasename)
    RMSE = pd.DataFrame(RMSE,index = kinasename)
    MAE = pd.DataFrame(MAE,index = kinasename)
    best_params = pd.DataFrame(best_params, index = kinasename)
    output = pd.concat([R2,RMSE,MAE,best_params], axis = 1)
    output = pd.DataFrame(output)
    output.to_csv("XGBoost_result.csv")


def LR_process(RNAseq, kinaseact):
    print("Linear Regression!")
    R2 = []
    RMSE = []
    MAE = []
    best_params=[]
    kinasename = kinaseact.index.values
    kinaseact = kinaseact.values
    for i in range(len(kinasename)):
        X_train, X_test, Y_train, Y_test = train_test_split(RNAseq, kinaseact[i], test_size=0.25, random_state=10)
        lm = LinearRegression()
        lm.fit(X_train, Y_train)
        rfe = RFE(lm) 
        folds = KFold(n_splits = 5, shuffle = True, random_state = 100)
        param_grid = [{'n_features_to_select': list(range(1, 50))}]
        model_cv = GridSearchCV(estimator = rfe,
                        param_grid = param_grid, 
                        scoring= 'r2', 
                        cv = folds, 
                        verbose = 1,
                        return_train_score=True)   
        model_cv.fit(X_train, Y_train)
        best_estimator=model_cv.best_estimator_
        y_pred = best_estimator.predict(X_test)
        print(y_pred)
        ##compute R2
        pred_R2 = r2_score(Y_test,y_pred)
        R2.append(pred_R2)
        ##compute RMSE
        pred_RMSE = np.sqrt(mean_squared_error(Y_test,y_pred))
        RMSE.append(pred_RMSE)
        ##compute MAE
        pred_MAE = mean_absolute_error(Y_test,y_pred)
        MAE.append(pred_MAE)
        s = pickle.dumps(model_cv)
        with open('%s.model' % kinasename[i], 'wb+') as f:
            f.write(s)
    R2 = pd.DataFrame(R2,index = kinasename)
    RMSE = pd.DataFrame(RMSE,index = kinasename)
    MAE = pd.DataFrame(MAE,index = kinasename)
    best_params = pd.DataFrame(best_params, index = kinasename)
    output = pd.concat([R2,RMSE,MAE,best_params], axis = 1)
    output = pd.DataFrame(output)
    print(output)
    output.to_csv("linear_regression_result.csv")

def SVM_process(RNAseq, kinaseact):
    print("SVM!")
    R2 = []
    RMSE = []
    MAE = []
    best_params=[]
    kinasename = kinaseact.index.values
    kinaseact = kinaseact.values
    for i in range(len(kinasename)):
        X_train, X_test, Y_train, Y_test = train_test_split(RNAseq, kinaseact[i], test_size=0.25, random_state=10)
        regr = SVR() 
        regr.fit(X_train, Y_train)
        folds = KFold(n_splits = 5, shuffle = True, random_state = 100)
        grid_params=dict(C=np.logspace(-2,2,10),gamma=np.logspace(-2,2,10))
        model_cv = GridSearchCV(estimator = regr,
                        param_grid = grid_params, 
                        scoring= 'r2', 
                        cv = folds, 
                        verbose = 1,
                        return_train_score=True)   
        model_cv.fit(X_train, Y_train)
        best_estimator=model_cv.best_estimator_
        best_params.append(model_cv.best_params_)
        y_pred = best_estimator.predict(X_test)
        print(y_pred)
        ##compute R2
        pred_R2 = r2_score(Y_test,y_pred)
        R2.append(pred_R2)
        ##compute RMSE
        pred_RMSE = np.sqrt(mean_squared_error(Y_test,y_pred))
        RMSE.append(pred_RMSE)
        ##compute MAE
        pred_MAE = mean_absolute_error(Y_test,y_pred)
        MAE.append(pred_MAE)
        s = pickle.dumps(model_cv)
        with open('%s.model' % kinasename[i], 'wb+') as f:
            f.write(s)
    R2 = pd.DataFrame(R2,index = kinasename)
    RMSE = pd.DataFrame(RMSE,index = kinasename)
    MAE = pd.DataFrame(MAE,index = kinasename)
    best_params = pd.DataFrame(best_params, index = kinasename)
    output = pd.concat([R2,RMSE,MAE,best_params], axis = 1)
    output = pd.DataFrame(output)
    print(output)
    output.to_csv("SVM_result.csv")

      


def RF_process(RNAseq,kinaseact):
    print("Random Forest!")
    R2 = []
    RMSE = []
    MAE = []
    best_params=[]
    kinasename = kinaseact.index.values
    kinaseact = kinaseact.values
    for i in range(len(kinasename)):
        X_train, X_test, Y_train, Y_test = train_test_split(RNAseq, kinaseact[i], test_size=0.25, random_state=10)
        regr = RandomForestRegressor(random_state=0)
        regr.fit(X_train, Y_train)
        folds = KFold(n_splits = 5, shuffle = True, random_state = 100)
        grid_params=dict(max_depth=[4,5,6],n_estimators=[100,200])
        model_cv = GridSearchCV(estimator = regr,
                        param_grid = grid_params, 
                        scoring= 'r2', 
                        cv = folds, 
                        verbose = 1,
                        return_train_score=True)   
        model_cv.fit(X_train, Y_train)
        best_estimator=model_cv.best_estimator_
        best_params.append(model_cv.best_params_)
        y_pred = best_estimator.predict(X_test)
        print(y_pred)
        ##compute R2
        pred_R2 = r2_score(Y_test,y_pred)
        R2.append(pred_R2)
        print(pred_R2)
        ##compute RMSE
        pred_RMSE = np.sqrt(mean_squared_error(Y_test,y_pred))
        RMSE.append(pred_RMSE)
        ##compute MAE
        pred_MAE = mean_absolute_error(Y_test,y_pred)
        MAE.append(pred_MAE)
        s = pickle.dumps(model_cv)
        with open('%s.model' % kinasename[i], 'wb+') as f:
            f.write(s)
    R2 = pd.DataFrame(R2,index = kinasename)
    RMSE = pd.DataFrame(RMSE,index = kinasename)
    MAE = pd.DataFrame(MAE,index = kinasename)
    best_params = pd.DataFrame(best_params, index = kinasename)
    output = pd.concat([R2,RMSE,MAE,best_params], axis = 1)
    output = pd.DataFrame(output)
    print(output)
    output.to_csv("RandomForest_result.csv")


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
        Core model algorithm. Supported: "Xgboost",  "linear regression", "RF", "SVM", "KNeighbors"
    in_dir : str
        Relative path or name of the desired input directory
    """

    def __init__(
        self,
        core_model,
    ):
        self.results = pd.DataFrame()
        self.core_model = core_model
    
    def prepare(self,
            kinase_list,
            cancer_type,
            select_all,
            core_model,
            ):
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
            
            kinase_act = pd.read_csv("../kinase/%s_kinase.csv" % cancer_type, index_col = 0)
            kin = pd.DataFrame()
            for i in range(len(kinase_list)):
                if kinase_list[i] in kinase_act.index.values:
                    kin = pd.concat([kin, kinase_act.loc[kinase_list[i],:]], axis = 1)
            if not select_all:
                RNAseq = data.preprocess_data_L1000(cancer_type)
            else:
                RNAseq = data.preprocess_data_all(cancer_type)
            print("preprocessed bulk RNA-seq data\n")
            print(RNAseq)
            print("preprocessed kinase data\n")
            kin = kin.T
            print(kin)

            if(core_model=="XGboost"):
                xgboost_process(RNAseq,kin)
            if(core_model=="linear_regression"):
                LR_process(RNAseq,kin)
            if(core_model=="RF"):
                RF_process(RNAseq,kin)
            if(core_model=="SVM"):
                SVM_process(RNAseq,kin)
            



    def predict(
        self, file, 
        kinase_name,
        select_all
    ):
        """Method to make predictions.
        
        Make predictions for unknown datasets
        we used pretrained model to predict kinase activities"""
        predict_result = pd.DataFrame()
        dataMat = pd.read_csv(file, index_col = 0)
        nrow = dataMat.shape[0]
        ncol = dataMat.shape[1]
        
        if not select_all:
            L1000 = pd.read_table("../L1000.txt")
            newdata = pd.DataFrame()
            for i in range(L1000.shape[0]):
                if L1000.iloc[i].values in dataMat.index.values:
                    newdata = pd.concat([newdata,dataMat.loc[L1000.iloc[i], :]],axis=0)
                else:
                    genename=''.join(L1000.iloc[i].values)
                    newdata.loc[genename,:]=[0.5]*newdata.shape[1]
            outdata = (newdata - newdata.min())/(newdata.max() - newdata.min())
            newdata = outdata
            print(newdata)
        else:
            newdata=dataMat
        sample_name = newdata.columns.values
        prediction_result = pd.DataFrame()
        prediction_result = []
        for i in range(len(kinase_name)):
            f = open('%s.model' % kinase_name[i], 'rb')
            s = f.read()
            model = pickle.loads(s)
            result = model.predict(newdata.T)
            prediction_result.append(result)
        prediction_result = pd.DataFrame(prediction_result, columns=sample_name)
        print(prediction_result)
        prediction_result.to_csv("predict_result.csv")

  


