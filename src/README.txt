Using the code:

DATA: Average_Times.txt

In the data directory, there is a file called Average_Times.txt.
The file is formatted so that each line contains the settings and subsequent
runtime of one (or multiple) experiment(s).  The first column is the 
number of mappers, x, an integer; the second column is the number of reducers,
y, an integer; and the third column is the average runtime, z (floating point).
Every combination of integer values (x, y) with 2 <= x, y <= 32 were 
sampled at least once.  Some of the values were sampled as many as 3 times.
This is due to some overlap in our sample scheduling which was caused
by our computational job periodically timing out.  The data from samples that
had not finished when the job was cancelled have been removed.



IMPORTANT NOTE: All of the python scripts described here work for Python 2.7
and require a few (common) libraries, numpy, stats, matplotlib, for example.




SAMPLING: get_random_sample.py

To simulate taking a random sampling of parameters, one can simply select a 
random set of lines from the data file.  We wrote a small python script to
do this, get_random_sample.py.


	USAGE: (from the command line)

	python get_random_sample <path_to_data>/Average_Times.txt N > Sampled_Data.txt

	N should be the number of points to randomly sample.
The script will print out the selected lines in the same format as they
were read in.  In the example above, the output is captured in a text file
"Sampled_Data.txt" that will be re-used for later illustration.





K-FOLD CROSS VALIDATION: K-fold.py

Once a file containing experimental observations has been created 
(e.g. Sampled_Data.txt) you can apply a K-fold cross validation.
The purpose of this code segment is to determine the ideal degree of a
polynomial surface to fit the data.  For a description of how K-fold
fitting works, see the last half of Section 2.2 in our ICCS 2015 paper
(see LICENSE.txt for full citation).


	USAGE: (from the command line)

	python K-fold.py <fileName> <degree> <number_of_folds>

The code will partition the data in fileName into number_of_folds pieces
of size as equal as possible.  number_of_folds - 1 of these folds will
be used to build a surface of degree degree.  Then, the SSE will
be computed against the reserved fold.  This is repeated so that each
fold is used once for testing.  The code outputs the degree and the average SSE
from each of the computations.

In practice, one will want to run this code several (e.g. 10 or more) times
for each polynomial degree under consideration.  The ideal degree surface,
given the data collected, will be the degree with the lowest average SSE.
In the event that there are two degrees with approximately equal average SSE,
it may be preferable to select the degree with a smaller typical variation in SSE.





BUILDING THE SURROGATE SURFACE: construct_model.py

Once it has been determined what degree surface is desired, one can
generate the full surface using construct_model.py.
It requires the same data file that was used by K-fold.py (e.g. Sampled_Data.txt).

	
	USAGE (from the command line)

	python construct_model.py <fileName> <degree> > Model_Surface.txt

The code uses the data in fileName (assumed to be in format x-val, y-val, z-val)
to generate a surface of degree degree.  It computes the z-values of the surface for
all the (x,y) pairs where 2 <= x, y <= 32.  It applies a 99% confidence interval
to each point.  The result is that the code generates a predicted value z-mean and a bound
from the confidence interval.  The bound from the confidence interval is always larger than z-mean
and reflects an expected worst-case upperbound on that point.  The data is printed out in the form:

x 	y	z-mean	z-upperbound

When we optimize, we look for the smallest value of z-upperbound.  The confidence iterval
penalizes areas in the domain that contain very few sampled points.  The effect of which is
that we are "more confident" of expected runtimes in areas near points that we have sampled.
