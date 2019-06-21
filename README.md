# QuIBL
## Quantifying Introgression via Branch Lengths

``QuIBL`` is a program for detecting introgressed loci and characterizing introgression events in a species network. It takes a set of gene trees as input and estimates the proportion of introgressed loci as well as the timing of an introgression pulse triplet by triplet.

## Installation
Currently, ``QuIBL`` is distributed as a python script. To run it, ensure that you are using Python v2.7, and have installed the following dependencies:
``joblib``, ``ete3``, ``itertools``, ``sys``, ``numpy``, ``math``, ``ConfigParser``, ``csv``, ``joblib``, and ``multiprocessing``.

## Running QuIBL
To test ``QuIBL``, download the repository and type ``python QuIBL.py ./Small_Test_Example/sampleInputFile.txt``. This will run a small provided sample data file. To test your own data, change the ``treefile`` setting to be the path to a file with your Newick trees.

## Input File Settings
``treefile`` : The path to the trees to be analyzed.

``numdistributions`` : The number of branch length distributions in the mixture to test. For now, only two is supported (this corresponds to one ILS and one non-ILS distribution).

``likelihoodthresh`` : The maximum change in likelihood allowed for the gradient ascent search for theta to stop.

``numsteps`` : The number of total EM steps. For thousands of trees, we reccomend trying around 50.

``gradascentscalar`` : The factor to shrink the stepsize when a gradient ascent step fails.

``totaloutgroup`` : The name of the ultimate outrgroup of your sample. All trees are assumed to be rooted using this taxon.

``multiproc`` : Accepts ``True`` or ``False`` and either turns multiprocessing on or off.

``OutputPath`` : Where the output gets written.
