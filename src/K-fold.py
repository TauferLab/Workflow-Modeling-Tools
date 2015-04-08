### Copyright (c)
### 	2015 by The University of Delaware
### 	Contributors: Travis Johnston
### 	Affiliation: Global Computing Laboratory, Michela Taufer PI
### 	Url: http://gcl.cis.udel.edu/, https://github.com/TauferLab
### 
### All rights reserved.
### 
### Redistribution and use in source and binary forms, with or without modification,
### are permitted provided that the following conditions are met:
### 
### 	1. Redistributions of source code must retain the above copyright notice, 
### 	this list of conditions and the following disclaimer.
### 
### 	2. Redistributions in binary form must reproduce the above copyright notice,
### 	this list of conditions and the following disclaimer in the documentation
### 	and/or other materials provided with the distribution.
### 
### 	3. If this code is used to create a published work, one or both of the
### 	following papers must be cited.
### 
### 		Travis Johnston, Mohammad Alsulmi, Pietro Cicotti, Michela Taufer, "Performance
### 		Tuning of MapReduce Jobs Using Surrogate-Based Modeling," International
### 		Conference on Computational Science (ICCS) 2015.
### 
### 		M. Matheny, S. Herbein, N. Podhorszki, S. Klasky, M. Taufer, "Using Surrogate-
### 		based Modeling to Predict Optimal I/O Parameters of Applications at the Extreme
### 		Scale," In Proceedings of the 20th IEEE International Conference on Parallel and
### 		Distributed Systems (ICPADS), December 2014.
### 
### 	4.  Permission of the PI must be obtained before this software is used
### 	for commercial purposes.  (Contact: taufer@acm.org)
### 
### THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
### ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
### WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
### IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
### INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
### BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
### DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
### LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
### OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
### OF THE POSSIBILITY OF SUCH DAMAGE.

import sys
import numpy as np
import numpy.matlib
import random


### This function computes n Choose r

def nCr(n, r):
	if n < r:
		return 0
	retVal = 1
	for i in xrange(r):
		retVal = (retVal*(n-i))/(i+1)
	
	return retVal



### This polynomial returns a float corresponding to the evaluation of a polynomial at a point.
### The function requires the x-value in the polynomial, the y-value in the polynomial,
### the degree of the polynomial and the list of coefficients.
### The coefficients must be in order:
### a_0 + a_1 x + a_2 y + a_3 x^2 + a_4 xy + a_5 y^2 + a_6 x^3 + a_7 x^2 y + a_8 xy^2 + a_9 y^3 + etc.
### In general, this specific order is NOT required, just that the order of the coefficients be internally consistent.

def Evaluate_Poly(x, y, degree, coefficients):
	retVal = 0
	index = 0

	for DEG in xrange(degree+1):
		for yPOW in xrange(DEG+1):
			retVal += coefficients[index]*(x**(DEG-yPOW))*(y**yPOW)
			index += 1
	return float(retVal)




### Seed the random number generator, None causes the seed to be the system time.
random.seed(None)

### This is the filename containing the observed data points.
### The k-fold cross validation will be applied to all the points in this file.
fileName = str(sys.argv[1])


### These lists will contain the x, y, and z (time) data contained in the file, fileName.
data_X = []
data_Y = []
data_Z = []


file = open(fileName, 'r')
for line in file:
	numbers = line.split()
	data_X.append(int(numbers[0]))		### The first entry (integer) is the x value
	data_Y.append(int(numbers[1]))		### The second entry (integer) is the y value
	data_Z.append(float(numbers[2])) 	### The third entry (float) is the runtime (z value)
file.closed


### The user specifies what degree of polynomial surface to build (for the fitting)
### The goal of this code snippet is to determine if this specified degree is better than others.
DEGREE = int(sys.argv[2])

### This user specified parameter dictates how many pieces the data in fileName is partitioned into.
NUM_FOLDS = int(sys.argv[3])

### The number of monomials in a 2-variate polynomial of degree n
nMonomials = nCr(DEGREE+2, 2)


N = len(data_X)			### Number of data points in the file
indices = range(N)
Folds = []


### Generate the random partition of points into the K-folds
for i in xrange(NUM_FOLDS-1):
	Rand_smpl = sorted(random.sample(indices, N/NUM_FOLDS))
	Folds.append(Rand_smpl)
	
	for j in Rand_smpl:
		indices.remove(j)
Folds.append(indices)


### SAVE_Beta will be a list of the coefficients for the different polynomials built.
### Each fold of the data will be used once as a testing set.
### The number of different polynomials we will generate is the same as the number of folds.
SAVE_Beta = []


for testing_fold in xrange(NUM_FOLDS):
	### Use Folds[testing_fold] as the testing fold,
	### Build the best fit surface from the other points.
	nL = N - len(Folds[testing_fold])

	X = np.matlib.ones((nL, nMonomials))
	Beta = np.matlib.zeros((nMonomials, 1))
	Z = np.matlib.zeros((nL, 1))

	data_point = 0
	learning_folds = list(range(NUM_FOLDS))
	learning_folds.remove(testing_fold)
	for i in learning_folds:
		for j in Folds[i]:
			### Add row to X corresponding to values of data point data_X[j], data_Y[j]
			index = 0
			for degree_monomial in xrange(DEGREE+1):
				for yDEG in xrange(degree_monomial+1):
					X[(data_point, index)] = (data_X[j]**(degree_monomial-yDEG))*(data_Y[j]**(yDEG))
					index += 1
			Z[(data_point, 0)] = data_Z[j]
			data_point += 1
	
	### We have now constructed the matrix X and Z
	### Perform the multiplication that gives us the Beta vector (of coefficients)

	Beta = (X.getT()*X).getI()*X.getT()*Z	#(X^T X)^(-1) * X^T * Z
	SAVE_Beta.append(Beta)


### Compute the total SSE for each surface against the testing points, and average the results.

Total_SSE = 0
for testing_fold in xrange(NUM_FOLDS):
	for point in Folds[testing_fold]:
		Total_SSE += (data_Z[point]-Evaluate_Poly(data_X[point], data_Y[point], DEGREE, SAVE_Beta[testing_fold]))**2

print DEGREE, "\t", Total_SSE/NUM_FOLDS 
