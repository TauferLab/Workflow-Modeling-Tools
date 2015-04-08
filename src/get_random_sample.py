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
import random

### Set the random seed (None will use the system time for a seed).
random.seed(None)

### The name of the file containing the data.
### It is expected that this file be in the format of:
###
### x-value_1 y-value_1 z-value_1
### ...
### x-value_i y-value_i z-value_i
### ...
### x-value_n y-value_n z-value_n
###
fileName = str(sys.argv[1])

### These lists will contain the values on each line of the file.
entire_data_X = []
entire_data_Y = []
entire_data_Z = []

file = open(fileName, 'r')
for line in file:
	numbers = line.split()
	entire_data_X.append(int(numbers[0]))	### The first number in the line is the x-value, expected to be an integer.
	entire_data_Y.append(int(numbers[1]))	### The second number on the line is the y-value, expected to be an integer.
	entire_data_Z.append(float(numbers[2]))	### The third number on the line is the z-value (observed run-time), and could be a float.
file.closed


### The user-supplied number of random samples to make.
N = int(sys.argv[2])


### The number of data points in the file.
Ne = len(entire_data_X)
indices = range(Ne)

### Pick N random indices from the list of indices.
Rand_smpl = sorted(random.sample(indices, N))

### Select out the N data points corresponding to the randomly chosen indices and print them.
for i in Rand_smpl:
	print entire_data_X[i], entire_data_Y[i], entire_data_Z[i]
