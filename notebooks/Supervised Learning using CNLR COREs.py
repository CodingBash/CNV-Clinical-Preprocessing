
# coding: utf-8

# # Supervised Learning of Drug Response using CORES from Copy Number Log Ratio

# ### Import Python source code

# In[1]:


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

from sklearn.linear_model import Lasso
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
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

# In[8]:

labeled_matrix_training_set = pd.read_csv("../mlOutput/coreTrainingSet_7_31_2018_1.csv")
labeled_matrix_training_set.columns.values[0] = "sampleId"
labels = list(range(0,5))


# In[4]:

pprint(labeled_matrix_training_set.copy().head())


# In[18]:

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import Imputer

preprocessing_pipeline = Pipeline([
    ('imputer', Imputer(strategy="median")) # TODO: Temporary fix for scaling
    ,('std_scaler', StandardScaler())
])

processed_training_set = preprocessing_pipeline.fit_transform(labeled_matrix_training_set.copy().drop(["sampleId"], axis=1))


# In[19]:

print(final_training_set)


# In[21]:

final_training_set = pd.DataFrame(data=final_training_set, columns = labeled_matrix_training_set.columns[1:])
print(final_training_set.copy().head())


# In[30]:

def abline(slope, intercept):
    """Plot a line from slope and intercept"""
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--')


# ### Visualize ML results using Linear Regression

# In[33]:

for label in labels:
    # Remove uneeded labels
    selected_training_set = final_training_set.iloc[:, list([label]) + list(range(5,final_training_set.shape[1]))].copy()
    selected_training_set = selected_training_set[np.isfinite(selected_training_set.iloc[:,1])]
    training_set, testing_set = split_train_test(selected_training_set)
    model_data = training_set.copy().drop([training_set.columns[0], training_set.columns[1]], axis = 1)
    model_labels = training_set.iloc[:,1]
    lasso = LinearRegression()
    lasso.fit(model_data, model_labels)
    model_test_data = testing_set.copy().drop([training_set.columns[0], training_set.columns[1]], axis = 1)
    model_test_labels = testing_set.iloc[:,1]
    predictions = lasso.predict(model_test_data)
    mse = mean_squared_error(model_test_labels, predictions)
    rmse = np.sqrt(mse)
    print(rmse)
    r = scipy.stats.pearsonr(model_test_labels, predictions)
    t = scipy.stats.spearmanr(model_test_labels, predictions)
    print("Pearson: " + str(r))
    print("Spearman: " + str(t))
    plt.plot(model_test_labels, predictions, 'bo')
    abline(1,0)
    plt.show()


# ### Visualize ML results using Random Forest Regressor

# In[32]:

for label in labels:
    # Remove uneeded labels
    selected_training_set = final_training_set.iloc[:, list([label]) + list(range(5,final_training_set.shape[1]))].copy()
    selected_training_set = selected_training_set[np.isfinite(selected_training_set.iloc[:,1])]
    training_set, testing_set = split_train_test(selected_training_set)
    model_data = training_set.copy().drop([training_set.columns[0], training_set.columns[1]], axis = 1)
    model_labels = training_set.iloc[:,1]
    lasso = RandomForestRegressor(n_estimators=500, max_leaf_nodes=16, n_jobs=-1)
    lasso.fit(model_data, model_labels)
    model_test_data = testing_set.copy().drop([training_set.columns[0], training_set.columns[1]], axis = 1)
    model_test_labels = testing_set.iloc[:,1]
    predictions = lasso.predict(model_test_data)
    mse = mean_squared_error(model_test_labels, predictions)
    rmse = np.sqrt(mse)
    print(rmse)
    r = scipy.stats.pearsonr(model_test_labels, predictions)
    t = scipy.stats.spearmanr(model_test_labels, predictions)
    print("Pearson: " + str(r))
    print("Spearman: " + str(t))
    plt.plot(model_test_labels, predictions, 'bo')
    abline(1,0)
    plt.show()

