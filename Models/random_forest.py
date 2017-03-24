from __future__ import division
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.datasets import make_classification
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import confusion_matrix

def accuracy(pred, actual):
    """Returns percentage of correctly classified labels"""
    return sum(pred==actual) / len(actual)

def main():

	# Import data and split out labels
	train = pd.read_csv('../data/train.csv')
	test = pd.read_csv('../data/test.csv')

	# Split into labels, names and data
	y_train = train['class']
	names_train = train['name']
	X_train = train.drop(['class', 'name', 'sequence'], axis=1)

	y_test = test['class']
	names_train = test['name']
	X_test = test.drop(['class', 'name', 'sequence'], axis=1)

	# Optimize a random forest model using grid search
	rf = RandomForestClassifier()

	param_grid = {
	    'n_estimators': [500], 
	    'max_depth': [15, 20, 25],
	    'max_features': [8, 10, 12],
	    'n_jobs': [30],
	    'min_samples_leaf': [1, 2, 3]
	}

	grid_rf = GridSearchCV(rf, param_grid, cv=2, verbose=3)
	grid_rf.fit(X_train, y_train)
	print("#-------- DONE WITH GRID SEARCH.")
	best_model = grid_rf.best_estimator_
	best_params = grid_rf.best_params_ 
	scores = grid_rf.grid_scores_
	print(best_params)

	# Calculate cross validation accuracy
	rf = RandomForestClassifier(n_estimators=500, criterion='entropy', max_features=8, max_depth=15,
	                            min_samples_leaf = 1, bootstrap=True, oob_score=True, n_jobs=30, random_state=0)
	print(np.mean(cross_val_score(rf, X_train, y_train, cv=5)))

	# # Fit to full training data and table feature importances
	# rf = rf.fit(X_train, y_train)
	# importances = rf.feature_importances_
	# importance = pd.DataFrame(importances, index=X_train.columns, columns=["importance"])
	# # importance.sort('importance', ascending=0)

	# # Print train and test accuracy
	# y_train_pred = rf.predict(X_train)
	# y_test_pred = rf.predict(X_test)
	# print("Training Accuracy = %f" % accuracy(y_train_pred, y_train))
	# print("Test Accuracy = %f" % accuracy(rf.predict(X_test), y_test))

	# confusion_matrix(np.array(y_test), np.array(y_test_pred))

if __name__ == "__main__":
	main()