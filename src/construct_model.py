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

import scipy.stats as stats
import sys
import numpy as np
import numpy.matlib
import random

### This function return the binomial coefficient n Choose r.
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



### fileName is the name of the file containing the sampled data points.
### In a typical use-case, this is a subset of the entire parameter space.
fileName = str(sys.argv[1])


### These lists will hold the X, Y, and Z data contained in the file.
### It is assumed that each line of the file contains a data point: x-value, y-value, z-value (space separated).
data_X = []
data_Y = []
data_Z = []

file = open(fileName, 'r')
for line in file:
	numbers = line.split()
	data_X.append(int(numbers[0]))		### Read in the x-coordinate (integer)
	data_Y.append(int(numbers[1]))		### Read in the y-coordinate (integer)
	data_Z.append(float(numbers[2]))	### Read in the z-coordinate (observed runtime) float.

file.closed



DEGREE = int(sys.argv[2])		### This is the degree of the polynomial surface to build.
nMonomials = nCr(DEGREE+2, 2)	### This counts the number of monomials that will appear in the polynomial of degree DEGREE in two variables.
N = len(data_X)					### The number of data points
df = N - nMonomials - 1			### Number of degrees of freedom, a parameter to pass to stats.t



### These are the matrices X, Z, and $\beta$ as defined in the paper:
### "Performance Tuning of MapReduce Jobs Using Surrogate-Based Modeling" by Travis Johnston, Mohammad Alsulmi, Pietro Cicotti, and Michela Taufer.
### The description is found in Section 2.2, equations (2).
X = np.matlib.ones((N, nMonomials))
Beta = np.matlib.zeros((nMonomials, 1))
Z = np.matlib.zeros((N, 1))


### Construct the values in matrix X and Z.
for data_point in xrange(N):	
	index = 0
	for degree_monomial in xrange(DEGREE+1):
		for yDEG in xrange(degree_monomial+1):
			X[(data_point, index)] = (data_X[data_point]**(degree_monomial-yDEG))*(data_Y[data_point]**(yDEG))
			index += 1
	Z[(data_point, 0)] = data_Z[data_point]
	

### Perform the multiplication that gives us the Beta vector (of coefficients)
### We want to compute Beta = (X^T X)^-1 X^T Z.
### We compute it in two steps.
X_transX = (X.getT()*X).getI()		### X_transX = (X^T X)^-1 and is computed separately since we will recycle this portion over and over again later.
Beta = X_transX*X.getT()*Z


### The Beta values are the coefficients of the best fit polynomial surface, now we want to apply a confidence interval.
### We have already saved (X^T X)^-1... now we need to compute the SSE
SSE = 0
for i in xrange(len(data_X)):
	### At data_X[i], data_Y[i] we sampled a time data_Z[i]... compare data_Z[i] to the polynomial surface at data_X[i], data_Y[i]
	SSE += ( Z[i] - Evaluate_Poly(data_X[i], data_Y[i], DEGREE, Beta) )**2


### We construct the confidence interval at each point in the feasible region (in this case x in [2, 32] y in [2, 32])
### using the method on pages 14-15 here: http://courses.washington.edu/b515/l6.pdf
x_0 = np.matlib.zeros((nMonomials, 1))
conf_level=.99
for x in xrange(2, 33):
	for y in xrange(2, 33):
		### First, compute x_0 at this point (x,y)
		index = 0
		for degree_monomial in xrange(DEGREE+1):
			for yDEG in xrange(degree_monomial+1):
				x_0[index, 0] = (x**(degree_monomial-yDEG))*(y**yDEG)
				index += 1
		### The [0,0] here is used because SSE and x_0 are 2D matrices with 1 row and 1 column each.
		var = (SSE[0,0]*(x_0.getT()*X_transX*x_0)[0, 0])**.5		### This is the variance of the model at (x,y)
		mean = Evaluate_Poly(x, y, DEGREE, Beta)					### This is the mean of the model at (x,y)... actually the value of the surface itself.

		upper = stats.t(df, loc=mean, scale=var).ppf(conf_level)
		### mean is the value our model predicts, and upper is the larger value in the confidence interval.
		### essentially, we are confident that the actual value would fall somewhere in the confidence interval, 
		### and we say that we are confident that the time required for the run is at most "upper".
		### When we find optimal parameters, we look for the "best worst case," i.e. the smallest "upper."

		### One could use the smallest "mean" value... but there are problems because the surface typically will drop off rapidly on (at least one) boundary of
		### the region.  This leads to negative values.  It is can also lead to strong predictions in areas where little or no sampling has been done.
		print x, y, mean, upper
