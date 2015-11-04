import datetime
import logging

import numpy as np
import pandas as pd

from sklearn.utils import shuffle
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold, SelectFdr, RFECV
from sklearn.cross_validation import StratifiedKFold
from sklearn.metrics import confusion_matrix, roc_auc_score, f1_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC

from . import util
from . import window_extraction
from . import sklearn_extensions

LOGGER = logging.getLogger('ML')

DEFAULT_CLASSIFIERS = [
    RandomForestClassifier(n_estimators = 180, n_jobs = -2, class_weight = 'auto'),
    SVC(kernel = 'rbf', class_weight = 'auto', probability = True, cache_size = 1500),
]

DEFAULT_TRANSFORMER = StandardScaler(copy = False)

DEFAULT_FEATURE_SELECTOR = sklearn_extensions.FeatureSelectionPipeline([
    VarianceThreshold(0.01),
    SelectFdr(alpha = 0.1),
])

RFECV_FEATURE_SELECTION_DEFAULT_CLASSIFIER = sklearn_extensions.RandomForestClassifierWithCoef(n_estimators = 400, n_jobs = -2, class_weight = 'auto')

# Use a constant seed
SEED = 1812
np.random.seed(SEED)

# Silence annoying pandas warnings
pd.options.mode.chained_assignment = None

class WindowClassifier(object):

    '''
    Classifies windows extracted with their features.
    '''

    def __init__(self, raw_classifier, used_features, transformer = DEFAULT_TRANSFORMER):
        '''
        @param raw_classifier (sklearn classifier): A raw classifier trained against a dataset of windows.
        @param used_features (list of strings): The names of the features used while training the raw classifier (after feature selection
        has been applied).
        @param transformer (sklearn transformer, optional, default sklearn.preprocessing.StandardScaler): The exact same transformer used
        when training the raw classifier (providing any other transformer is expected to give very poor results).
        '''
        self.raw_classifier = raw_classifier
        self.used_features = used_features
        self.transformer = transformer

    def classify_windows(self, windows_data_frame, proba = False):

        '''
        Classifies windows extracted with their features, given in a CSV format. Obviously, these windows don't need to have annotations/labels.
        Even if labels are given, they will be ignored.
        @param windows_data_frame (pandas.DataFrame):
            A data frame of the windows' CSV.
        @param proba (default False):
            Whether to return predicted probabilities (floats from between 0 to 1) or binary labels (0s or 1s).
        @return:
            A numpy array of the predicted labels for the given windows. The length of the returned array will correspond to the number of
            windows in the given data frame.
        '''

        if len(windows_data_frame) == 0:
            return np.empty(shape = 0)
        else:

            X = windows_data_frame[self.used_features].values
            X = self._transform(X)

            if proba:
                return self.raw_classifier.predict_proba(X)[:,1]
            else:
                return self.raw_classifier.predict(X)

    def test_performance(self, windows_data_frame, drop_only_almost_positives = False, drop_duplicates = True, scoring_method = f1_score):
        '''
        Tests the performance of this trained classifier, that was originally trained on a certain dataset, on a new dataset. The given dataset
        should be windows extracted with their features and annotations, given in a CSV format.
        The documentation of this method is partial and lacks some important details, as it's very similar to train_window_classifier, which
        already has a detailed documentation. Therefore, make sure to read the documentation of the other method in order to understand the full
        meaning of all the parameters.
        @param windows_data_frame (pandas.DataFrame):
            A data frame of the windows' CSV.
        @param drop_only_almost_positives (boolean, default False):
            Whether to drop only almost positive windows in the dataset before evaluating the performance of this classifier against it.
        @param drop_duplicates (boolean, default True):
            Whether to drop duplicating windows in the dataset, based on their neighbourhood property, before evaluating the performance
            of this classifier against it.
        @param scoring_method (function, default sklearn.metrics.f1_score):
            A scoring method to evaluate the classifiers by, just like in train_window_classifier.
        @return:
            A tuple of scores measuring the performance of this classifier against the given dataset in the format (score, roc, sensitivity,
            precision, specificity, cm), just like in train_window_classifier.
        '''
        LOGGER.info('Testing ' + str(type(self.raw_classifier)))
        features, X, y = get_windows_data(windows_data_frame, drop_only_almost_positives, drop_duplicates, self.transformer, \
                features = self.used_features)
        LOGGER.info('Predicting %d records...' % len(X))
        y_pred = self.raw_classifier.predict(X)
        return _get_prediction_scores(y, y_pred, scoring_method)

    def _transform(self, X):
        if self.transformer is None:
            return X
        else:
            return self.transformer.fit_transform(X)

class PeptidePredictor(object):

    '''
    Uses a trained window classifier to predicts annotations for new peptide.
    '''

    def __init__(self, window_classifier, window_extraction_params = window_extraction.WindowExtractionParams()):
        '''
        @param window_classifier (WindowClassifier):
            The trained window classifier to use.
        @param window_extraction_params (WindowExtractionParams, default params by default):
            The exact same parameters that have been used to extract the windows on which the window classifier has been trained (providing
            any other set of parameters is expected to result very unpleasant errors).
        '''
        self.window_classifier = window_classifier
        self.window_extraction_params = window_extraction_params

    def predict_annotations(self, seq, extra_tracks_data = {}, proba = False):

        '''
        Predicts the annotations of a peptide.
        @param seq (string):
            The amino-acid sequence of the peptide to predict the annotations for, given in a 20 amino-acid alphabet.
        @param extra_tracks_data (dict, empty by default):
            A dictionary for providing extra tracks of the given peptide. Must receive the data for all the tracks that have been used to
            extract the windows for training this classifier. Specifically, if this predictor relies on a feature that relies on a certain
            track, then this track must be provided here. The given dictionary should map from track names to their sequence.
        @param proba (default False):
            Whether to return predicted probabilities (floats from between 0 to 1) or binary labels (0s or 1s).
        @return:
            If proba = False, will return a binary string (of 0's and 1's) representing the predicted annotations for the given peptide. If
            proba = True, will return a list of floats (between 0 to 1), representing the predicted probabilities. Either way, the length of
            the returned string/list will correspond to the length of the provided peptide sequence.
        '''

        length = len(seq)
        windows_csv = window_extraction.extract_windows_from_seq(seq, extra_tracks_data = extra_tracks_data, \
                window_extraction_params = self.window_extraction_params)
        windows_data_frame = pd.read_csv(windows_csv)
        window_indices = windows_data_frame['window_hot_index'].values

        labels = self.window_classifier.classify_windows(windows_data_frame, proba = proba)
        annotation_mask = [0] * length

        for window_index, label in zip(window_indices, labels):
            if window_index >= 0 and window_index < length:
                annotation_mask[window_index] = label

        if proba:
            return map(float, annotation_mask)
        else:
            return ''.join(map(str, annotation_mask))

def train_window_classifier(windows_data_frame, classifiers = DEFAULT_CLASSIFIERS, drop_only_almost_positives = False, \
        drop_duplicates = True, transformer = DEFAULT_TRANSFORMER, feature_selector = DEFAULT_FEATURE_SELECTOR, n_folds = 5, \
        scoring_method = f1_score, select_best = True):

    '''
    Trains a window classifier using a CSV of windows with extracted features and annotations/labels (obtained by either
    window_extraction.extract_windows_from_file with extract_annotations = True or window_extraction.extract_windows_from_seq with a given
    annotation_mask). The evaluation of the classifiers will be based on the kfold procedure, during which various metrics will be calculated.
    The final training of the classifier will be based on the entire data set.
    @param windows_data_frame (pandas.DataFrame):
        A data frame of the windows' CSV.
    @param classifiers (list of sklearn classifiers, default Gaussian-kernel SVM and random forest):
         A list of classifiers to try training independently, from which the best classifier can be chosen.
    @param drop_only_almost_positives (boolean, default False):
        Whether to drop "only almost positive" windows in the dataset. An only almost positive window is a window with a false label in its
        hot index, but with a true label in either of the flanking indices. In some learning scenarios, the labeling of the residues (i.e.
        annotations) isn't so important in a strict manner, and it only matters whether larger regions contain a positive label. It's especially
        important in cases that the actual used dataset is only accurate up to +/-1 shifts of the labels. In such scenarios, using this parameter
        might enhance performance.
    @param drop_duplicates (boolean, default True):
        Whether to drop duplicating windows in the dataset, based on their neighbourhood property.
    @param transformer (sklearn transformer, optional, default sklearn.preprocessing.StandardScaler):
        A preprocessing transformer to use for the data before starting the kfold evaluation and final training of the classifiers. If None, will
        not perform any preprocessing  transformation.
    @param feature_selector (sklearn feature selector, optional, default a pipeline of VarianceThreshold and SelectFdr):
        A feature selection procedure to apply during both the kfold evaluation and final training of each classifier. If None, will not perform
        feature selection (i.e. will use all features). Note that the given feature selector must implement the get_support method (hence sklearn's
        builtin Pipeline object cannot be used; if you want to pipeline then use FeatureSelectionPipeline of this project).
    @param n_folds (int, default 5):
        The number of folds to use during the kfold evaluation procedure.
    @param scoring_method (function, default sklearn.metrics.f1_score):
        A scoring method to evaluate the classifiers by. Expecting a method that receives two parameters (y_true and y_pred) and returns a float
        score. This score will be calculated for all classifiers, in addition to other metrics. Also, if select_best is set to True, this score
        will be used in order to choose the best classifier.
    @param select_best (boolean, default True):
        Whether to return only the best evaluated classifier or all of them.
    @return:
        For each classifier, will return a tuple of the trained WindowClassifier object and its metrics, as evaluated during the kfold procedure.
        The metrics are also a tuple of floats in the format (score, roc, sensitivity, precision, specificity, cm), where: score is the score
        calculated by scoring_method; roc is Area Under the Curve (AUC); cm stands for the 2X2 confusion matrix of the results. If select_best is
        set to True, will return only the tuple of the best classifier (based on the score). Otherwise, will return a list of tuples for all the
        classifiers, sorted by their score in a descending order.
    '''

    features, X, y = get_windows_data(windows_data_frame, drop_only_almost_positives, drop_duplicates, transformer)
    window_classifiers_and_results = []

    for classifier in classifiers:
        kfold_results = _get_classifier_kfold_results(classifier, X, y, n_folds, feature_selector, scoring_method)
        window_classifier = _get_trained_window_classifier(classifier, features, X, y, feature_selector, transformer)
        window_classifiers_and_results += [(window_classifier, kfold_results)]

    window_classifiers_and_results.sort(key = lambda window_classifier_and_results: window_classifier_and_results[1], reverse = True)

    if select_best:
        best_classifier, best_results = window_classifiers_and_results[0]
        LOGGER.info('The best classifier is %s with score %f.' % (str(type(best_classifier.raw_classifier)), best_results[0]))
        return best_classifier, best_results
    else:
        return window_classifiers_and_results

def get_top_features(windows_data_frame, drop_only_almost_positives = False, drop_duplicates = True, transformer = DEFAULT_TRANSFORMER, \
        classifier = RFECV_FEATURE_SELECTION_DEFAULT_CLASSIFIER, n_folds = 3, step = 0.05, scoring = 'f1'):

    '''
    Using sklearn.feature_selection.RFECV model in order to find the top features of given windows with features, given in a CSV format.
    @param windows_data_frame (pandas.DataFrame):
        A data frame of the windows' CSV.
    @param drop_only_almost_positives (boolean, default False):
        Same as in train_window_classifier.
    @param drop_duplicates (boolean, default True):
        Whether to drop duplicating windows in the dataset, based on their neighbourhood property, prior to RFECV.
    @param transformer (sklearn transformer, optional, default sklearn.preprocessing.StandardScaler):
        A preprocessing transformer to use for the data before applying RFECV. If None, will not perform any preprocessing transformation.
    @param classifier (sklearn classifier, default a special version of random forest suitable for RFECV):
        The classifier to use as the estimator of RFECV.
    @param n_folds (int, default 2):
        The n_folds to use in the kfold cross-validation as part of the RFECV process.
    @param step (default 0.05):
        See sklearn.feature_selection.RFECV
    @param scoring (default 'f1'):
        See sklearn.feature_selection.RFECV
    @return:
        A list of the top features, each represented as a string.
    '''

    features, X, y = get_windows_data(windows_data_frame, drop_only_almost_positives, drop_duplicates, transformer)
    kfold = StratifiedKFold(y, n_folds = n_folds, shuffle = True, random_state = SEED)
    rfecv = RFECV(estimator = classifier, cv = kfold, step = step, scoring = scoring)
    rfecv.fit(X, y)
    return util.apply_mask(features, rfecv.support_)

def get_windows_data(windows_data_frame, drop_only_almost_positives = False, drop_duplicates = True, transformer = DEFAULT_TRANSFORMER, \
        features = None):

    '''
    Extracts numeric vectorial data in numpy format, suitable for applying standard sklearn models on, from a CSV of windows with features.
    @param windows_data_frame (pandas.DataFrame):
        A data frame of the windows' CSV.
    @param drop_only_almost_positives (boolean, default False):
        Same as in train_window_classifier.
    @param drop_duplicates (boolean, default True):
        Whether to drop duplicating windows in the dataset, based on their neighbourhood property.
    @param transformer (sklearn transformer, optional, default sklearn.preprocessing.StandardScaler):
        A transformer to apply on the data (X). If None, will not perform any preprocessing transformation.
    @param features (list of strings, optional):
        The names of the features to extract from each window. If None, will extract all the features that appear in the given CSV.
    @return:
        A tuple comprised of:
        1. features - A list of strings corresponding to the names of the features extracted from the data.
        2. X - A numpy matrix of the extracted data points. Each row in the matrix represents a window, and each column a feature.
        3. y - A numpy array of binary integer values (0s and 1s), corresponding to the label of the extracted data points (windows). The
        length of y is equal to the number of rows in X.
    '''

    LOGGER.info('Given a data frame of %d records X %d columns.' % windows_data_frame.shape)

    if drop_only_almost_positives:
        windows_data_frame = windows_data_frame[windows_data_frame['window_only_almost_positive'] == 0]
        LOGGER.info('Dropped only almost positives. %d records remained.' % len(windows_data_frame))

    if drop_duplicates:
        # When we remove duplicates, we want to give priority to positives
        windows_data_frame.sort(columns = 'window_label', ascending = False, inplace = True)
        windows_data_frame.drop_duplicates(subset = 'window_neighbourhood', inplace = True)
        LOGGER.info('Dropped duplicates. %d records remained.' % len(windows_data_frame))

    if features is None:
        features = [header for header in windows_data_frame.columns if header not in window_extraction.META_WINDOW_HEADERS]

    LOGGER.info('%d features to process.' % len(features))
    X = windows_data_frame[features].values
    y = windows_data_frame['window_label'].values
    X, y = shuffle(X, y, random_state = SEED)

    if transformer is not None:
        X = transformer.fit_transform(X)
        LOGGER.info('Transformed the data.')

    LOGGER.info('Final data: samples = %d, features = %d' % X.shape)
    return features, X, y

def _get_classifier_kfold_results(classifier, X, y, n_folds, feature_selector, scoring_method):

    LOGGER.info('Estimating ' + str(type(classifier)))
    time_before = datetime.datetime.now()

    y_pred = _predict_using_kfold(X, y, classifier, n_folds, feature_selector)
    score, roc, sensitivity, precision, specificity, cm = _get_prediction_scores(y, y_pred, scoring_method)

    time_diff = datetime.datetime.now() - time_before
    LOGGER.info('Finished estimating. Took %d seconds' % int(time_diff.total_seconds()))

    LOGGER.info('score = %f, roc = %f, sensitivity = %f, precision = %f, specificity = %f' % (score, roc, sensitivity, precision, specificity))
    LOGGER.info('Confusion matrix:' + '\n' + str(cm))
    return score, roc, sensitivity, precision, specificity, cm

def _get_trained_window_classifier(classifier, features, X, y, feature_selector, transformer):
    LOGGER.info('Training ' + str(type(classifier)))
    X_reduced = feature_selector.fit_transform(X, y)
    used_features = util.apply_mask(features, feature_selector.get_support())
    classifier.fit(X_reduced, y)
    return WindowClassifier(classifier, used_features, transformer)

def _predict_using_kfold(X, y, classifier, n_folds, feature_selector):

    kfold = StratifiedKFold(y, n_folds = n_folds, shuffle = True, random_state = SEED)
    y_pred = np.zeros(len(y))

    for i, fold in enumerate(kfold):

        LOGGER.info('Running fold %d/%d...' % (i + 1, n_folds))

        train_indices, test_indices = fold
        X_train = X[train_indices]
        X_test = X[test_indices]
        y_train = y[train_indices]

        if feature_selector is not None:
            feature_selector.fit(X_train, y_train)
            X_train = feature_selector.transform(X_train)
            X_test = feature_selector.transform(X_test)
            LOGGER.info('Selected features. Remained with %d features.' % X_train.shape[1])

        classifier.fit(X_train, y_train)
        y_pred[test_indices] = classifier.predict(X_test)

    return y_pred

def _get_prediction_scores(y_true, y_pred, scoring_method):

    cm = confusion_matrix(y_true, y_pred, labels = [0, 1])
    roc = roc_auc_score(y_true, y_pred)
    score = scoring_method(y_true, y_pred)

    tn = float(cm[0][0])
    tp = float(cm[1][1])
    fp = float(cm[0][1])
    fn = float(cm[1][0])
    n = tn + fp
    p = tp + fn

    sensitivity = tp / p
    specificity = tn / n
    precision = tp / (tp + fp)

    return score, roc, sensitivity, precision, specificity, cm
