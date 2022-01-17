"""
# ----------------------------------------------------------------------------
# LAMBDAdemo.py   Demo script to test the LAMBDA software and familiarize with
#                 the associated functions
#
#                 The input ambiguity vectors and variance-covariance matrices
#                 are given below. You need to comment/uncomment the parts
#                 needed. The three data examples have been retrieved by the
#                 Matlab Version 3.0 software (small.mat, sixdim.mat, large.mat)
#
# Created by: Dimitrios Psychas
# Date:       01-NOV-2019
# Modified:   -
# ----------------------------------------------------------------------------
"""

# Import libraries
import numpy as np
import LAMBDA

# Input vc-matrix and ambiguity vector 

#==============================================================================
# small.mat (taken from MATLAB)
#==============================================================================
#Qahat = [[6.2900, 5.9780, 0.5440],[5.9780, 6.2920, 2.3400],[0.5440, 2.3400, 6.2880]]
#Qahat = np.array(Qahat)
#ahat  = [5.4500, 3.1000, 2.9700]
#ahat  = np.array(ahat)
#==============================================================================
#==============================================================================

#==============================================================================
# sixdim.mat (taken from MATLAB)
#==============================================================================
#Qahat = [[ 19068.8559508787, -15783.972282037000,  -17334.200587597500,  14411.9239749603000,  10055.71700893590, -14259.29529038720],\
#         [-15783.9722820370,  59027.703840981500,   38142.692753110200,    562.7173880246450, -13830.08559606760,  27373.42630130190],\
#         [-17334.2005875975,  38142.692753110200,   28177.565389352800,  -7000.5022049704500, -11695.86740593060,  21886.16806305320],\
#         [ 14411.9239749603,    562.717388024645,   -7000.502204970450,  15605.5082283690000,   5039.70281815470,  -9648.96530646004],\
#         [ 10055.7170089359, -13830.085596067600,  -11695.867405930600,   5039.7028181547000,   6820.77250679480,  -6880.24051213224],\
#         [-14259.2952903872,  27373.426301301900,   21886.168063053200,  -9648.9653064600400,  -6880.24051213224,  23246.54896269450]]
#Qahat = np.array(Qahat)         
#ahat  = [-28490.8566886116, 65752.6299198198, 38830.3666554972, 5003.70833517778, -29196.0699104593, -297.658932458787]
#ahat  = np.array(ahat)
#==============================================================================
#==============================================================================

#==============================================================================
# large.mat (taken from MATLAB)
#==============================================================================
Qahat=[[19068.8559508787,	-15783.9722820370,	-17334.2005875975,	14411.9239749603,	10055.7170089359,	-14259.2952903872,	14858.8484050976,	-12299.1993741839,	-13507.1694819930,	11230.0704356810,	7835.62344938376,	-11111.1393808147],\
       [-15783.9722820370,	59027.7038409815,	38142.6927531102,	562.717388024645,	-13830.0855960676,	27373.4263013019,	-12299.1993747356,	45995.6129934030,	29721.5785731468,	438.480887460148,	-10776.6902686912,	21329.9423774758],\
       [-17334.2005875975,	38142.6927531102,	28177.5653893528,	-7000.50220497045,	-11695.8674059306,	21886.1680630532,	-13507.1694826246,	29721.5785738846,	21956.5440705992,	-5454.93697674992,	-9113.66310734779,	17054.1567378091],\
       [14411.9239749603,	562.717388024645,	-7000.50220497045,	15605.5082283690,	5039.70281815470,	-9648.96530646004,	11230.0704356773,	438.480887731461,	-5454.93697653627,	12160.1358938811,	3927.04096307733,	-7518.67445855756],\
       [10055.7170089359,	-13830.0855960676,	-11695.8674059306,	5039.70281815470,	6820.77250679480,	-6880.24051213224,	7835.62344947055,	-10776.6902682086,	-9113.66310687634,	3927.04096320258,	5314.88728015545,	-5361.22656658847],\
       [-14259.2952903872,	27373.4263013019,	21886.1680630532,	-9648.96530646004,	-6880.24051213224,	23246.5489626945,	-11111.1393809211,	21329.9423779274,	17054.1567375591,	-7518.67445829957,	-5361.22656681708,	18114.1936088811],\
       [14858.8484050976,	-12299.1993747356,	-13507.1694826246,	11230.0704356773,	7835.62344947055,	-11111.1393809211,	11578.3237340013,	-9583.79156943782,	-10525.0669778554,	8750.70438611838,	6105.68076067050,	-8658.03053539344],\
       [-12299.1993741839,	45995.6129934030,	29721.5785738846,	438.480887731461,	-10776.6902682086,	21329.9423779274,	-9583.79156943782,	35840.7376978353,	23159.6717654859,	341.673569568934,	-8397.42083743563,	16620.7344703582],\
       [-13507.1694819930,	29721.5785731468,	21956.5440705992,	-5454.93697653627,	-9113.66310687634,	17054.1567375591,	-10525.0669778554,	23159.6717654859,	17108.9956804894,	-4250.60009053988,	-7101.55551676305,	13288.9534523001],\
       [11230.0704356810,	438.480887460148,	-5454.93697674992,	12160.1358938811,	3927.04096320258,	-7518.67445829957,	8750.70438611838,	341.673569568934,	-4250.60009053988,	9475.43086798586,	3060.03207008500,	-5858.70721928591],\
       [7835.62344938376,	-10776.6902686912,	-9113.66310734779,	3927.04096307733,	5314.88728015545,	-5361.22656681708,	6105.68076067050,	-8397.42083743563,	-7101.55551676305,	3060.03207008500,	4141.47090961885,	-4177.57899193454],\
       [-11111.1393808147,	21329.9423774758,	17054.1567378091,	-7518.67445855756,	-5361.22656658847,	18114.1936088811,	-8658.03053539344,	16620.7344703582,	13288.9534523001,	-5858.70721928591,	-4177.57899193454,	14114.9563601479]]
Qahat = np.array(Qahat)       
ahat=[-28490.8566886116,65752.6299198198,38830.3666554972,5003.70833517778,-29196.0699104593,-297.658932458787,-22201.0284440701,51235.8374755528,30257.7809603224,3899.40332138829,-22749.1853575113,-159.278779870217       ]
ahat = np.array(ahat)

# Test the various methods of IAR

# Method 1/a - ILS with 2 cands
afixed,sqnorm,Ps,Qzhat,Z,nfixed,mu = LAMBDA.main(ahat,Qahat,1)
print("--------------")
print("Method 1/a - ILS with 2 integer candidates")
print("--------------")
print("\nVector of fixed integer ambiguities")
print(afixed)
print("\nVector of squared norms")
print(sqnorm)
print("\nAmbiguity success rate")
print(Ps)
print("\nVariance-covariance matrix of decorrelated ambiguities")
print(Qzhat)
print("\nZ-transformation matrix")
print(Z)
print("\nNumber of fixed ambiguities")
print(nfixed)
print("\nThreshold value of Ratio Test")
print(mu)

# Method 1/b - ILS with 3 cands
afixed,sqnorm,Ps,Qzhat,Z,nfixed,mu = LAMBDA.main(ahat,Qahat,1,3)
print("--------------")
print("Method 1/b - ILS with 3 integer candidates")
print("--------------")
print("\nVector of fixed integer ambiguities")
print(afixed)
print("\nVector of squared norms")
print(sqnorm)
print("\nAmbiguity success rate")
print(Ps)
print("\nVariance-covariance matrix of decorrelated ambiguities")
print(Qzhat)
print("\nZ-transformation matrix")
print(Z)
print("\nNumber of fixed ambiguities")
print(nfixed)
print("\nThreshold value of Ratio Test")
print(mu)

# Method 2 - Integer rounding
afixed,sqnorm,Ps,Qzhat,Z,nfixed,mu = LAMBDA.main(ahat,Qahat,2)
print("--------------")
print("Method 2 - IR")
print("--------------")
print("\nVector of fixed integer ambiguities")
print(afixed)
print("\nVector of squared norms")
print(sqnorm)
print("\nAmbiguity success rate")
print(Ps)
print("\nVariance-covariance matrix of decorrelated ambiguities")
print(Qzhat)
print("\nZ-transformation matrix")
print(Z)
print("\nNumber of fixed ambiguities")
print(nfixed)
print("\nThreshold value of Ratio Test")
print(mu)

# Method 3 - Integer bootstrapping
afixed,sqnorm,Ps,Qzhat,Z,nfixed,mu = LAMBDA.main(ahat,Qahat,3)
print("--------------")
print("Method 3 - IB")
print("--------------")
print("\nVector of fixed integer ambiguities")
print(afixed)
print("\nVector of squared norms")
print(sqnorm)
print("\nAmbiguity success rate")
print(Ps)
print("\nVariance-covariance matrix of decorrelated ambiguities")
print(Qzhat)
print("\nZ-transformation matrix")
print(Z)
print("\nNumber of fixed ambiguities")
print(nfixed)
print("\nThreshold value of Ratio Test")
print(mu)

# Method 4/a - Partial ambiguity resolution with default P0=99.5%
afixed,sqnorm,Ps,Qzhat,Z,nfixed,mu = LAMBDA.main(ahat,Qahat,4,2)
print("--------------")
print("Method 4/a - PAR using default P0")
print("--------------")
print("\nVector of fixed integer ambiguities")
print(afixed)
print("\nVector of squared norms")
print(sqnorm)
print("\nAmbiguity success rate")
print(Ps)
print("\nVariance-covariance matrix of decorrelated ambiguities")
print(Qzhat)
print("\nZ-transformation matrix")
print(Z)
print("\nNumber of fixed ambiguities")
print(nfixed)
print("\nThreshold value of Ratio Test")
print(mu)

# Method 4/b - Partial ambiguity resolution with P0=95%
afixed,sqnorm,Ps,Qzhat,Z,nfixed,mu = LAMBDA.main(ahat,Qahat,4,2,0.950)
print("--------------")
print("Method 4/b - PAR with P0=95%")
print("--------------")
print("\nVector of fixed integer ambiguities")
print(afixed)
print("\nVector of squared norms")
print(sqnorm)
print("\nAmbiguity success rate")
print(Ps)
print("\nVariance-covariance matrix of decorrelated ambiguities")
print(Qzhat)
print("\nZ-transformation matrix")
print(Z)
print("\nNumber of fixed ambiguities")
print(nfixed)
print("\nThreshold value of Ratio Test")
print(mu)

# Method 5/a - ILS + FFRT with default P0=0.1%
afixed,sqnorm,Ps,Qzhat,Z,nfixed,mu = LAMBDA.main(ahat,Qahat,5)
print("-------------------")
print("Method 5/a - ILS+FFRT with default P0=0.1%")
print("-------------------")
print("\nVector of fixed integer ambiguities")
print(afixed)
print("\nVector of squared norms")
print(sqnorm)
print("\nAmbiguity success rate")
print(Ps)
print("\nVariance-covariance matrix of decorrelated ambiguities")
print(Qzhat)
print("\nZ-transformation matrix")
print(Z)
print("\nNumber of fixed ambiguities")
print(nfixed)
print("\nThreshold value of Ratio Test")
print(mu)

# Method 5/b - ILS + FFRT with P0=1%
afixed,sqnorm,Ps,Qzhat,Z,nfixed,mu = LAMBDA.main(ahat,Qahat,5,2,0.01)
print("-------------------")
print("Method 5/b - ILS+FFRT with P0=1%")
print("-------------------")
print("\nVector of fixed integer ambiguities")
print(afixed)
print("\nVector of squared norms")
print(sqnorm)
print("\nAmbiguity success rate")
print(Ps)
print("\nVariance-covariance matrix of decorrelated ambiguities")
print(Qzhat)
print("\nZ-transformation matrix")
print(Z)
print("\nNumber of fixed ambiguities")
print(nfixed)
print("\nThreshold value of Ratio Test")
print(mu)