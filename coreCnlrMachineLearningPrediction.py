# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 12:21:38 2018

@author: bbece
"""

from __future__ import division, print_function, unicode_literals
import numpy as np
import os

from pprint import pprint
np.random.seed(42)

import matplotlib
import matplotlib.pyplot as plt
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
            
import pandas as pd

labeled_matrix_training_set = pd.read_csv("mlOutput/coreTrainingSet_7_26_2018_1.csv")
labeled_matrix_training_set.columns.values[0] = "sampleId"

labels = 1:6
for label in labels:
    # Remove uneeded labels
    selected_training_set = labeled_matrix_training_set.iloc[:, list([0]) + list([label]) + list(range(6,labeled_matrix_training_set.shape[1]))].copy()
    
    


split_train_test()
