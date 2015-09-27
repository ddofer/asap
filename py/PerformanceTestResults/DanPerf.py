'''
Dan: Get model performance on train, test, FeatureSelection
'''

import sys
sys.path += [r'E:\Dropbox\Dropbox\Protein Cleavage Prediction\asap\py']
sys.path += [r'E:\Dropbox\Dropbox\Protein Cleavage Prediction\asap\py\asap']
sys.path += [r'E:\Dropbox\Dropbox\Protein Cleavage Prediction\asap\py\cleavepred']

import pandas as pd
import asap
from asap.classification import *

#External package:
from mlxtend.classifier import EnsembleClassifier

from sklearn.grid_search import GridSearchCV
from sklearn.svm import SVC
from sklearn.cross_validation import StratifiedKFold,cross_val_score,StratifiedShuffleSplit,cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC,   LinearSVC
from sklearn.feature_selection import VarianceThreshold,SelectFdr
from sklearn.pipeline import Pipeline
from sklearn.externals import joblib
from sklearn.metrics import confusion_matrix, roc_auc_score, f1_score
from sklearn import metrics
from sklearn.linear_model import LogisticRegressionCV
from sklearn.neighbors import KNeighborsClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.cross_validation import cross_val_predict

# Use a constant seed
SEED = 1812
np.random.seed(SEED)
# Silence annoying pandas warnings
pd.options.mode.chained_assignment = None


def get_sps_scores(y_true, y_pred):
    cm = confusion_matrix(y_true, y_pred, labels = [0, 1])
    tn = float(cm[0][0])
    tp = float(cm[1][1])
    fp = float(cm[0][1])
    fn = float(cm[1][0])
    n = tn + fp
    p = tp + fn
    # sensitivity = tp / p
    # precision = tp / (tp + fp)
    specificity = tn / n
    return specificity
    # return sensitivity, precision, specificity

def get_scores(y, y_pred,label=None):
    '''
    From GetResults.py code.
    Returns a dictionary of metrics for a given classification of the data (given by Cross_val_predict).
    y: list
        Class labels
    label: string
        Name of the classifier used
    '''

    precision,recall,fscore,support = metrics.precision_recall_fscore_support(y, y_pred,average='binary')
    print("Precision/PPV: %0.3f " % (precision))
    print("Recall/Sensitivity: %0.3f " % (recall))
    specificity = get_sps_scores(y_true=y, y_pred=y_pred)
    print("Specificity : %0.3f " % (specificity))
    roc_auc = metrics.roc_auc_score(y, y_pred)
    print("roc_auc: %0.3f " % (roc_auc))
    F1 = metrics.f1_score(y, y_pred,average='binary')
    print("F1: %0.4f  " % (F1))
    av_PR = metrics.average_precision_score(y, y_pred) # corresponds to the area under the precision-recall curve
    print("Average Precision (Prec-Recall AUC): %0.3f " % (av_PR))
    accuracy = metrics.accuracy_score(y, y_pred)
    print("Accuracy: %0.3f " % (accuracy))
    mcc = metrics.matthews_corrcoef(y, y_pred)
    print("MCC: %0.3f " % (mcc))

    results_dict = {'ROC_AUC':roc_auc,'F1':F1,'Accuracy':accuracy,
    'Precision':precision,'Recall':recall,
    'Specificity':specificity,
    'Average Precision':av_PR,'mcc':mcc
    }
    results_dict = {k:round(float(v),3) for k, v in results_dict.items()}
    return results_dict


SILLY_NUMBER = 110 #USed as a magic number; when debugging for speed
clf1 = LogisticRegressionCV(Cs=15,class_weight='auto')
clf2 = RandomForestClassifier(n_estimators=int(SILLY_NUMBER*1.5), max_features=45,bootstrap=True,class_weight='auto',n_jobs=2,criterion='gini',random_state=123)
clf3 = SVC(kernel = 'rbf', C = 3.798, gamma = 0.0, cache_size = 2200, class_weight = 'auto')#,probability=True)
# clf4 = KNeighborsClassifier(n_neighbors=5, weights='uniform')

# cclf1,cclf2,cclf3,cclf4 = CalibratedClassifierCV(clf1),CalibratedClassifierCV(clf2,cv=3),CalibratedClassifierCV(clf3,cv=4),CalibratedClassifierCV(clf4)
# clfs_calibrated = [cclf1,cclf2,cclf3,cclf4]
# eclf = EnsembleClassifier(clfs=clfs_calibrated, voting='soft', weights=[1,1,2,1])

eclf = EnsembleClassifier(clfs=[clf1,clf2,clf3], voting='hard')#, weights=[1,1,2,1])

outputFileName="CV_ensemble_bestTrainRes_param"
# CLASSIFIERS_TO_TEST = [SVC(kernel = 'rbf', C = 3.798, gamma = 0.0, cache_size = 2600, class_weight = 'auto')]

SCALER = StandardScaler(copy = False)
model = Pipeline([
    ('VAR_feature_selection',VarianceThreshold(0.04)),
    ('FDR', SelectFdr()),
    # ('svc', SVC(kernel = 'rbf', C = 3.798, gamma = 0.0, cache_size = 2000, class_weight = 'auto')),
    ('eclf',eclf ),
])


if __name__ == '__main__':

    np_simple = pd.read_csv(r'E:\Dropbox\Dropbox\Protein Cleavage Prediction\asap\data\cleavage\neuropred_dataset\window_simple_features.csv')
    np_adv = pd.read_csv(r'E:\Dropbox\Dropbox\Protein Cleavage Prediction\asap\data\cleavage\neuropred_dataset\window_advanced_features.csv')
    uni_simple = pd.read_csv(r'E:\Dropbox\Dropbox\Protein Cleavage Prediction\asap\data\cleavage\uniprot_dataset\window_simple_features.csv')
    uni_adv = pd.read_csv(r'E:\Dropbox\Dropbox\Protein Cleavage Prediction\asap\data\cleavage\uniprot_dataset\window_advanced_features.csv')

    datasets_loc = [np_simple,np_adv,uni_simple,uni_adv]
    datasets_names = ['np_simple','np_adv','uni_simple','uni_adv']
    datasets = zip(datasets_names,datasets_loc)

    def get_best_feature_selection_param(datasets):
        results = {}
        for data_name,data_set in datasets:
            print("\n Dataset: %s" %(data_name))

            feature_names, X, y = get_training_data(windows_data_frame=data_set, drop_only_almost_positives = True, drop_duplicates = True, transformer = DEFAULT_TRANSFORMER)

            # kfold_cv = StratifiedShuffleSplit(y, n_iter=7, test_size=0.4, random_state = SEED)
            kfold_cv = StratifiedKFold(y, n_folds=10, random_state = SEED,shuffle=True)

            'https://github.com/amueller/pydata-nyc-advanced-sklearn/blob/master/Chapter%201%20-%20Combining%20Pipelines%20and%20GridSearchCV.ipynb'

            alphas = np.array([0.3,0.1,0.05])
            # alphas = np.array([0.00000000000001])
            param_grid = {'FDR__alpha':alphas,'eclf__svc__C':[3.798, 11]}
            # param_grid = {'FDR__alpha':alphas}

            grid = GridSearchCV(estimator=model,
            param_grid=param_grid, iid=False,scoring='f1',refit=True,
            n_jobs=-2, cv =kfold_cv )
            grid.fit(X, y)
            print
            # summarize the results of the grid search
            print('Best Score:',grid.best_score_)
            print('grid.best_params_: ',grid.best_params_)
            # print
            # print('grid.best_estimator_ ',grid.best_estimator_)

            joblib.dump(grid.best_estimator_,data_name+'.pkl',compress=2)

            results[data_name] = grid.best_params_
            results[data_name].update({'F1':grid.best_score_})

        print("results dict:")
        print(results)
        res_df = pd.DataFrame(results)
        res_df.to_csv(outputFileName+".tsv", sep='\t')


    def get_test_perf(datasets):
        results={}

        'Load the data and labels for the NeuroPred and UNIprot data:'
        feature_names_np, X_np, y_np = get_training_data(windows_data_frame=np_simple, drop_only_almost_positives = True, drop_duplicates = True, transformer = DEFAULT_TRANSFORMER)
        feature_names_up, X_uni, y_uni = get_training_data(windows_data_frame=uni_simple, drop_only_almost_positives = True, drop_duplicates = True, transformer = DEFAULT_TRANSFORMER)

        feature_names_np_adv, X_np_adv, y_np_adv = get_training_data(windows_data_frame=np_adv, drop_only_almost_positives = True, drop_duplicates = True, transformer = DEFAULT_TRANSFORMER)
        feature_names_up_adv, X_uni_adv, y_uni_adv = get_training_data(windows_data_frame=uni_adv, drop_only_almost_positives = True, drop_duplicates = True, transformer = DEFAULT_TRANSFORMER)


        # 'load a previously trained classifier-pipeline. (Must be of same model type):'
        clf_uni = joblib.load('uni_simple.pkl')
        clf_np = joblib.load('np_simple.pkl')
        clf_uni_adv = joblib.load('uni_adv.pkl')
        clf_np_adv = joblib.load('np_adv.pkl')

        # predict_tests(joblib.load('uni_simple.pkl'),joblib.load('np_simple.pkl'),"Simple")
        # predict_tests(joblib.load('uni_adv.pkl'),joblib.load('np_adv.pkl'),"Advanced")

        modelClass="Simple"
        print(modelClass)
        print("NP Prediction on UniP:")
        y_uni_pred = clf_np.predict(X_uni)
        print("NeuroPred Prediction on UniProt Test Scores:")
        results[modelClass+' NeuroPred on UniProt'] = get_scores(y=y_uni,y_pred= y_uni_pred)

        print
        print("UNIPROT Prediction on NP:")
        y_np_pred = clf_uni.predict(X_np)
        print("\n UniProt  Prediction on NeuroPred:")
        results[modelClass+' UniProt on NeuroPred'] = get_scores(y=y_np, y_pred= y_np_pred)
        print
        print

        modelClass="Adv"
        print(modelClass)
        print("NP Prediction on UniP:")
        y_uni_pred_adv = clf_np_adv.predict(X_uni_adv)
        print("NeuroPred Prediction on UniProt Test Scores:")
        results[modelClass+' NeuroPred on UniProt'] = get_scores(y=y_uni_adv,y_pred= y_uni_pred_adv)

        print
        print("UNIPROT Prediction on NP:")
        y_np_pred_adv = clf_uni_adv.predict(X_np_adv)
        print("\n UniProt Prediction on NeuroPred:")
        results[modelClass+' UniProt on NeuroPred'] = get_scores(y=y_np_adv, y_pred= y_np_pred_adv)

        res_df = pd.DataFrame(results)
        res_df.to_csv('Train_TestVS'+".tsv", sep='\t')


    def get_train_perf(datasets):
        results={}

        model_top = Pipeline([
        ('VAR_feature_selection',VarianceThreshold(0.03)),
        ('FDR', SelectFdr(alpha=0.25)),
        ('eclf',eclf),
    ])
        for data_name,data_set in datasets:
            print("\n Dataset: %s" %(data_name))

            feature_names, X, y = get_training_data(windows_data_frame=data_set, drop_only_almost_positives = True, drop_duplicates = True, transformer = DEFAULT_TRANSFORMER)
            y_pred = cross_val_predict(estimator =model_top,X=X,y=y,cv=10, n_jobs=-2)
            results[data_name] = get_scores(y=y, y_pred= y_pred)

        res_df = pd.DataFrame(results)
        res_df.to_csv('Train_CV'+".tsv", sep='\t')


    # get_best_feature_selection_param(datasets)
    # get_test_perf(datasets)
    get_train_perf(datasets)
