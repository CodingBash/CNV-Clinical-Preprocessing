
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

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import Imputer

from sklearn.linear_model import Lasso
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import cross_val_score

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

pprint(labeled_matrix_training_set.copy().head())


# ## Visualize ML Results

# In[7]:

def abline(slope, intercept):
    """Plot a line from slope and intercept"""
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--')


# ### TODO: Work on cross validation

# In[23]:

label = 1

#
# Select label column
#
selected_training_set = labeled_matrix_training_set.iloc[:, list([label]) + list(range(6,labeled_matrix_training_set.shape[1]))].copy()

#
# Divide into training set and testing set
#
training_set, testing_set = split_train_test(selected_training_set, test_ratio = 0.20)

#
# Get model training information and preprocess
#
model_data = training_set.copy().drop([training_set.columns[0]], axis = 1)
model_labels = training_set.copy()[[training_set.columns[0]]]

training_feature_pipeline = Pipeline([
    ('training_feature_imputer', Imputer(strategy="median")),
    ('training_feature_standardizer', StandardScaler())
])

training_model_data_transformed = training_feature_pipeline.fit_transform(model_data)
model_data = pd.DataFrame(data=training_model_data_transformed, columns = model_data.columns) # TODO: Create a transformer for this

training_label_pipeline = Pipeline([
    ('training_label_imputer', Imputer(strategy="median")),
    ('training_label_standardizer', StandardScaler())
])

training_model_labels_transformed = training_label_pipeline.fit_transform(model_labels)
model_labels = pd.DataFrame(data=training_model_labels_transformed, columns = model_labels.columns)

#
# Fit the RandomForestRegressor model
#
rf_model = RandomForestRegressor(n_estimators=500, max_leaf_nodes=16, n_jobs=-1)
rf_model.fit(model_data, model_labels)

#
# TODO: To prevent data leakage, separate the scope after the model has been fit
#

#
# Get model testing information and preprocess
#
model_test_data = testing_set.copy().drop([testing_set.columns[0]], axis = 1)
model_test_labels = testing_set.copy()[[testing_set.columns[0]]]

testing_feature_pipeline = Pipeline([
    ('testing_feature_imputer', Imputer(strategy="median")),
    ('testing_feature_standardizer', StandardScaler())
])

testing_model_data_transformed = testing_feature_pipeline.fit_transform(model_test_data)
model_test_data = pd.DataFrame(data=testing_model_data_transformed, columns = model_test_data.columns)

testing_label_pipeline = Pipeline([
    ('testing_label_imputer', Imputer(strategy="median")),
    ('testing_label_standardizer', StandardScaler())
])

testing_model_labels_transformed = testing_label_pipeline.fit_transform(model_test_labels)
model_test_labels = pd.DataFrame(data=testing_model_labels_transformed, columns = model_test_labels.columns)



predictions = rf_model.predict(model_test_data)


rmse = np.sqrt(mean_squared_error(model_test_labels.copy().values.flatten(), predictions))
r = scipy.stats.pearsonr(model_test_labels.copy().values.flatten(), predictions)
t = scipy.stats.spearmanr(model_test_labels.copy().values.flatten(), predictions)

print("RMSE: " + str(rmse))
print("Pearson: " + str(r))
print("Spearman: " + str(t))

plt.plot(model_test_labels, predictions, 'bo')
abline(1,0)
plt.ylabel("Prediction")
plt.xlabel("Label")
plt.show()

pprint("Training Set: " + str(list(model_data.index)))
pprint("Test Set: " + str(list(model_test_data.index)))
    

scores = cross_val_score(rf_model, model_data, model_labels,
                         scoring = "neg_mean_squared_error", cv=10)
print(scores)
print(scores.mean())
print(scores.std())

