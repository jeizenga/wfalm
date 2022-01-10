#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  8 13:51:56 2022

@author: Jordan
"""

import sys
import numpy as np
import scipy.optimize as opt

# method from states, gish, altschul (1991)

def pam_mat(alpha):
    d = 0.99
    o = 0.01 / 3.0
    pam1 = np.array([[d, o, o, o], [o, d, o, o], [o, o, d, o], [o, o, o, d]])
    w, v = np.linalg.eig(pam1)
    return np.dot(v, np.dot(np.diag(np.power(w, alpha)), np.linalg.inv(v)))

    

if __name__ == "__main__":
    
    identity = .995 #float(sys.argv[1])
    
    res = opt.minimize_scalar(lambda a: (pam_mat(a)[0, 0] - identity)**2)
    
    pam = pam_mat(res.x)
    
    score = np.log(4.0 * pam)
    ratio = abs(score[0,0] / score[0,1])
    
    print("match / mismatch ratio for {}% identity should be ~{}".format(identity * 100.0, ratio))
    
    