
from xml.dom import INVALID_STATE_ERR
import data
import pandas as pd


def xgboost_process(cancer_type):
    RNA_seq = preprocess









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
        Core model algorithm. Supported: "Xgboost", "randomforest", "linear regression"
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
        self.core_model = init_core_model(core_model)
        self.in_str = in_dir
        self.fitted = False
        self.predicted = False
    
    def fit(self,
            kinase_list,
            save = False):
            """Method to fit the model.
            
            Use provided datasets, performing training and fit KBPRNA core_model.
            
            Parameters
            -------------
            kinase_list : list of kinases targeted for prediction
            save : bool
                whether or not to save the fitted core model in KBPRNA out_dir.
                Default: False"""

    def predict(
        self, name, save=False
    ):
        """Method to make predictions.
        
        Make predictions for unknown datasets
        1. """
