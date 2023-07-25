
"""
Utility functions for regression
"""
import numpy as np
from sklearn.linear_model import LinearRegression, RidgeCV
from sklearn.preprocessing import StandardScaler
from scipy.stats import pearsonr
import statsmodels.api as sm
from patsy import dmatrices


def coef_determination(targets, predictions):
    """
    Coefficient of determination R^2
    Arguments:
        targets (n_samples, ): ground truth
        predictions (n_samples): of a model
    
    Returns:
        float: coefficient of determination
    """
    residual_variance = np.mean((targets - predictions)**2)
    variance = np.mean((targets - targets.mean())**2)
    return 1 - residual_variance / variance

def regression_predictions(x, y, fitted_regression):
    """
    Fit and plot linear regression to predict y from x
    Just for visualization
    Arguments: 
        x (n_samples, 1): predictor vars
        y (n_samples, ) : target vars
        fitted_regression: scikit learn model fitted to x and y
    Returns:
        x_to_plot (2, ): x points to plot fit
        y_to_plot (2, ): corresponding y points
    """
    x_to_plot = np.array([x.min(), x.max()])
    y_to_plot = fitted_regression.predict(x_to_plot[:,None])
    return x_to_plot, y_to_plot

def leave_one_out_lr(x, y, regularize = False, return_plotvals = True):
    """
    Linear regression, leave one test sample out at a time
    Arguments:
        x (samples, dimensions): inputs 
        y (samples, ): targets
        regularize (bool, ): cross-validated L2 regularization?
    Returns:
        float: coefficient of determination
    """
    if regularize:
        lr = RidgeCV(alphas=(0.1, 1, 10, 100, 1000))
    else:
        lr = LinearRegression()
    n_samples = x.shape[0]
    predictions = np.zeros((n_samples, ))
    for n in range(n_samples):
        xtr = x[np.arange(n_samples)!=n]
        ytr = y[np.arange(n_samples)!=n]
        xte = x[n][None]
        yte = y[n]
        if regularize:
            scaler = StandardScaler()
            xtr = scaler.fit_transform(xtr)
            xte  = scaler.transform(xte)
        lr.fit(xtr, ytr)
        predictions[n] = lr.predict(xte)
    if return_plotvals:
        x_to_plot, y_to_plot = regression_predictions(x, y, lr)
        return coef_determination(y, predictions), x_to_plot, y_to_plot 
    else:
        return coef_determination(y, predictions)


