# Weighted least squares regression
# by Troyer Group, Institute for Theoretical Physics, ETH Zurich (C) 2012 - 2013

import numpy as np
import scipy.linalg as la

def wlstsq(X,y,w):
    '''weighted least squares regression'''
    '''X : matrix of predictor variables, shape (m,p)'''
    '''y : vector of response variables, shape (m,) or (m,k) for multi-valued observations'''
    '''w : vector of weights for the m observations, shape (m,)'''
    assert np.ndim(X) == 2 and np.ndim(w) == 1
    assert len(X) == len(y) and len(X) == len(w)
    # scale predictor and response variables by weights
    X = X * np.reshape(w, (-1,1))
    if np.ndim(y) == 1: y = y * w
    else:               y = y * np.reshape(w, (-1,1))
    fit,resid = la.lstsq(X,y)[:2]
    # coefficient of determination R^2 = 1 - |y - \hat{y}|^2 / |y - \bar{y}|^2
    r2 = 1 - resid / np.var(y, axis=0)
    return fit,r2

def wpolyfit(x,y,w,deg):
    '''Fit a polynomial of degree deg through the points (x[i],y[i]) with weights w[i].'''
    '''returns fit : polynomial coefficients in suitable order for np.polyval'''
    '''         r2 : coefficient of determination R^2'''
    # coefficient matrix for polynomial fit of degree deg
    X = np.transpose([np.power(x,n) for n in range(deg+1)[::-1]])
    fit,r2 = wlstsq(X,np.asanyarray(y),np.asanyarray(w))
    return fit,r2
