'''
This module contains extensions required by our project that we wish sklearn supported.
'''

from sklearn.ensemble import RandomForestClassifier

class RandomForestClassifierWithCoef(RandomForestClassifier):

    '''
    A small hack required to make sklearn.ensemble.RandomForestClassifier support sklearn.feature_selection.RFECV.
    '''

    def fit(self, *args, **kwargs):
        '''
        @see sklearn.ensemble.RandomForestClassifier.fit
        '''
        super(RandomForestClassifierWithCoef, self).fit(*args, **kwargs)
        self.coef_ = self.feature_importances_

class FeatureSelectionPipeline(object):
    
    '''
    Like sklearn.pipeline.Pipeline, but suitable for feature selection only.
    Unfortunately we can't use sklearn.pipeline.Pipeline as it is, because it doesn't have the get_support method that we need.
    This class isn't intended for general purpose, and we do not recommend using it outside the context of this project.
    '''
    
    def __init__(self, feature_selectors):
        
        '''
        @param feature_selectors (list of feature selectors):
            The list of feature selectors to pipeline together.
        '''
        
        if len(feature_selectors) == 0:
            raise Exception('Cannot pipeline an empty list of feature selectors')
            
        for feature_selector in feature_selectors:
            for method_name in ['fit', 'transform', 'fit_transform', 'get_support']:
                if not hasattr(feature_selector, method_name):
                    raise Exception('Feature selectors must have a %s method' % method_name)
        
        self.feature_selectors = feature_selectors
        
    def fit(self, X, y):
        
        for feature_selector in self.feature_selectors[:-1]:
            X = feature_selector.fit_transform(X, y)
            
        self.feature_selectors[-1].fit(X, y)
        
    def transform(self, X):
        
        for feature_selector in self.feature_selectors:
            X = feature_selector.transform(X)
            
        return X
        
    def fit_transform(self, X, y):
        
        for feature_selector in self.feature_selectors:
            X = feature_selector.fit_transform(X, y)
        
        return X
        
    def get_support(self):
        
        support = self.feature_selectors[0].get_support()
        
        for feature_selector in self.feature_selectors[1:]:
            support = _embed_vector_in_mask(feature_selector.get_support(), support)
            
        return support
        
def _embed_vector_in_mask(vector, mask):
    
    '''
    Embedding a vector inside the positive indices of a boolean mask.
    For example, if given the vector [x1, x2, x3] and the mask [0, 0, 0, 1, 0, 0, 1, 1], then the returned value will
    be [0, 0, 0, x1, 0, 0, x2, x3].
    Note that the number of 1's in the mask must be equal to the length of the vector.
    '''
    
    _validate_boolean(mask)
    
    if len(vector) != sum(mask):
        raise Exception('Cannot embed a vector of size %d in a mask with %d 1\'s' % (len(vector), sum(mask)))
    
    result = [0] * len(mask)
    vector_index = 0
    
    for i, flag in enumerate(mask):
        if flag:
            result[i] = vector[vector_index]
            vector_index += 1
            
    return result
    
def _validate_boolean(mask):
    for flag in mask:
        if flag not in [0, 1]:
            raise Exception('Expecting a boolean mask, given %s element' % repr(flag))
        
