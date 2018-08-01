
# coding: utf-8

# # Supervised Learning of Drug Response using CORES from Copy Number Log Ratio

# ### Import Python source code

# In[26]:


"""
Created on Thu Jul 26 12:21:38 2018

@author: bbece
"""

from __future__ import division, print_function, unicode_literals
import numpy as np
import os
from IPython.display import display, HTML

from pprint import pprint
np.random.seed(42)

import matplotlib
import matplotlib.pyplot as plt
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
            
import pandas as pd
import scipy

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import Imputer

from sklearn import decomposition

from sklearn.linear_model import Lasso
from sklearn.linear_model import Ridge
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import cross_val_score

from sklearn.tree import export_graphviz

import math



# ### Define method to split training and testing set

# In[2]:

# TODO: Manipulate test_ratio
def split_train_test(training_set, test_ratio = 0.33):
    row_count = training_set.shape[0]
    shuffled_indices = np.random.permutation(row_count)
    test_set_size = int(test_ratio * row_count)
    test_indices = shuffled_indices[:test_set_size]
    train_indices = shuffled_indices[test_set_size:]
    return training_set.iloc[train_indices], training_set.iloc[test_indices]


# ### Load training set matrix

# In[3]:

labeled_matrix_training_set = pd.read_csv("../mlOutput/coreTrainingSet_7_31_2018_1.csv")
labeled_matrix_training_set.columns.values[0] = "sampleId"
labels = list(range(1,6))


# In[4]:

display(labeled_matrix_training_set.copy().head())


# ## Visualize ML Results

# In[6]:

def abline(slope, intercept):
    """Plot a line from slope and intercept"""
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--')


# ### Visualize ML results using Linear Regression

# In[25]:

for label in labels:
    #
    # Select label column
    #
    selected_training_set = labeled_matrix_training_set.iloc[:, list([label]) + list(range(6,labeled_matrix_training_set.shape[1]))].copy()
    selected_training_set = selected_training_set[~np.isnan(selected_training_set.iloc[:,0])]
    #
    # Divide into training set and testing set
    #
    training_set, testing_set = split_train_test(selected_training_set, test_ratio = 0.33) # TODO: Use sklearn's train_test_split

    #
    # Get model training information and preprocess
    #
    model_data = training_set.copy().drop([training_set.columns[0]], axis = 1)
    model_labels = training_set.copy()[[training_set.columns[0]]]

    Ypipeline = Pipeline([
     ('imputer', Imputer(axis=0,strategy="median")),
     ('standardizer', StandardScaler()),
    ])

    model_labels_tr = Ypipeline.fit_transform(model_labels)

    XYpipeline = Pipeline([
            ('imputer', Imputer(axis=0,strategy="median")),
            ('standardizer', StandardScaler()),
            ('lasso_model', Ridge(alpha = 0.8))
    ])

    XYpipeline.fit(model_data, model_labels_tr)

    #
    # TODO: To prevent data leakage, separate the scope after the model has been fit
    #

    #
    # Get model testing information and preprocess
    #
    model_test_data = testing_set.copy().drop([testing_set.columns[0]], axis = 1)
    model_test_labels = testing_set.copy()[[testing_set.columns[0]]]

    model_test_labels_tr = Ypipeline.transform(model_test_labels)
    predictions = XYpipeline.predict(model_test_data)

    def imputer_inverse_transform(pre_data, post_data):
        na_indices = np.where(np.isnan(pre_data))[0]
        pprint(na_indices)
        post_data[na_indices] = float('NaN')
        return post_data


    predictions = Ypipeline.named_steps['standardizer'].inverse_transform(predictions)
    predictions = imputer_inverse_transform(model_test_labels, predictions)

    model_test_labels = model_test_labels.copy().values.flatten()
    model_test_labels = model_test_labels[~np.isnan(model_test_labels)]
    predictions = predictions[~np.isnan(predictions)]

    rmse = np.sqrt(mean_squared_error(model_test_labels, predictions))
    r = scipy.stats.pearsonr(model_test_labels, predictions)
    t = scipy.stats.spearmanr(model_test_labels, predictions)

    print("RMSE: " + str(rmse))
    print("Pearson: " + str(r))
    print("Spearman: " + str(t))

    plt.plot(model_test_labels, predictions, 'bo')
    abline(1,0)
    plt.ylabel("Prediction")
    plt.xlabel("Label")
    plt.show()

    scores = cross_val_score(XYpipeline, model_data, model_labels_tr,
                             scoring = "neg_mean_squared_error", cv=10)
    
    scores = Ypipeline.named_steps['standardizer'].inverse_transform(scores)
    
    print("CV Scores: " + str(scores))
    print("CV Mean: " + str(scores.mean()))
    print("CV STD: " + str(scores.std()))


# ### Visualize ML results using Random Forest Regressor

# In[32]:

for label in labels:
    #
    # Select label column
    #
    selected_training_set = labeled_matrix_training_set.iloc[:, list([label]) + list(range(6,labeled_matrix_training_set.shape[1]))].copy()
    selected_training_set = selected_training_set[~np.isnan(selected_training_set.iloc[:,0])]
    #
    # Divide into training set and testing set
    #
    training_set, testing_set = split_train_test(selected_training_set, test_ratio = 0.33) # TODO: Use sklearn's train_test_split

    #
    # Get model training information and preprocess
    #
    model_data = training_set.copy().drop([training_set.columns[0]], axis = 1)
    model_labels = training_set.copy()[[training_set.columns[0]]]

    Ypipeline = Pipeline([
     ('imputer', Imputer(axis=0,strategy="median")),
     ('standardizer', StandardScaler()),
    ])

    model_labels_tr = Ypipeline.fit_transform(model_labels)

    XYpipeline = Pipeline([
            ('pca', decomposition.PCA(n_components=15)),
            ('imputer', Imputer(axis=0,strategy="median")),
            ('standardizer', StandardScaler()),
            ('rf_model', RandomForestRegressor(n_estimators=1000, max_leaf_nodes=16, n_jobs=4)) # TODO: For now, hardcore the parameters
    ])

    XYpipeline.fit(model_data, model_labels_tr)

    #
    # TODO: To prevent data leakage, separate the scope after the model has been fit
    #

    #
    # Get model testing information and preprocess
    #
    model_test_data = testing_set.copy().drop([testing_set.columns[0]], axis = 1)
    model_test_labels = testing_set.copy()[[testing_set.columns[0]]]

    model_test_labels_tr = Ypipeline.transform(model_test_labels)
    predictions = XYpipeline.predict(model_test_data)

    def imputer_inverse_transform(pre_data, post_data):
        na_indices = np.where(np.isnan(pre_data))[0]
        pprint(na_indices)
        post_data[na_indices] = float('NaN')
        return post_data


    predictions = Ypipeline.named_steps['standardizer'].inverse_transform(predictions)
    predictions = imputer_inverse_transform(model_test_labels, predictions)

    model_test_labels = model_test_labels.copy().values.flatten()
    model_test_labels = model_test_labels[~np.isnan(model_test_labels)]
    predictions = predictions[~np.isnan(predictions)]

    rmse = np.sqrt(mean_squared_error(model_test_labels, predictions))
    r = scipy.stats.pearsonr(model_test_labels, predictions)
    t = scipy.stats.spearmanr(model_test_labels, predictions)

    print("RMSE: " + str(rmse))
    print("Pearson: " + str(r))
    print("Spearman: " + str(t))

    plt.plot(model_test_labels, predictions, 'bo')
    abline(1,0)
    plt.ylabel("Prediction")
    plt.xlabel("Label")
    plt.show()

    #scores = cross_val_score(XYpipeline, model_data, model_labels_tr,
                             #scoring = "neg_mean_squared_error", cv=10)
    #scores = Ypipeline.named_steps['standardizer'].inverse_transform(scores)
    
    #print("CV Scores: " + str(scores))
    #print("CV Mean: " + str(scores.mean()))
    #print("CV STD: " + str(scores.std()))




