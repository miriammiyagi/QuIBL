# QuIBL
Michael Miyagi (m_miyagi@g.harvard.edu) and Nate Edelman (nedelman@g.harvard.edu)

## Quantifying Introgression via Branch Lengths

``QuIBL`` is a program for detecting introgressed loci and characterizing introgression events in a species network. It takes a set of gene trees as input and estimates the proportion of introgressed loci as well as the timing of an introgression pulse triplet by triplet.

## Installation
Currently, ``QuIBL`` is distributed as a python script. To run it, ensure that you are using Python v2.7, and have installed the following dependencies:
``joblib``, ``ete3``, ``itertools``, ``sys``, ``numpy``, ``math``, ``ConfigParser``, ``csv``, ``joblib``, and ``multiprocessing``.

If you want to use the exact versions of ``ete3``,  ``joblib``, and  ``Cython``we used, type ``pip install -r quibl_requirements.txt`` in terminal while in the main directory. 

## Running QuIBL
To test ``QuIBL``, download and navigate to the repository and type ``python QuIBL.py ./Small_Test_Example/sampleInputFile.txt``. Any dependencies that may not be present can be installed using ``pip``. This will run a small provided sample data file and write an output CSV file to the same directory. To test your own data, change the ``treefile`` setting to be the path to a file with your Newick trees.

## Input File Settings
``treefile`` : The path to the trees to be analyzed.

``numdistributions`` : The number of branch length distributions in the mixture to test. For now, only two is supported (this corresponds to one ILS and one non-ILS distribution).

``likelihoodthresh`` : The maximum change in likelihood allowed for the gradient ascent search for theta to stop.

``numsteps`` : The number of total EM steps. For thousands of trees, we reccomend trying around 50.

``gradascentscalar`` : The factor to shrink the stepsize when a gradient ascent step fails.

``totaloutgroup`` : The name of the ultimate outrgroup of your sample. All trees are assumed to be rooted using this taxon.

``multiproc`` : Accepts ``True`` or ``False`` and either turns multiprocessing on or off.

``OutputPath`` : Where the output gets written.

``maxcores`` : The maximum number of cores QuIBL is allowed to use.

## Reading the Output

Currently ``QuIBL`` prints out a csv file with the following columns:

``triplet`` : The triplet analyzed, separated by underscores.

``outgroup`` : The outgroup for that line.

``C1, C2`` : The time during which two lineages are sequestered for that triplet topology for a 1 distribution and 2 distribution model respectively (see preprint for further details).

``mixprop1, mixprop2`` : The inferred mixing proportions for each distribution. ``mixprop2`` corresponds to the non-ILS component.

``lambda2Dist, lambda1Dist`` : The scaling factor required to go from substitutions per site (the input branch length unit) to coalescent units for a 2 and 1 distribution model. The 1 distribution model corresponds to the inverse of the mean branch length.

``BIC2Dist, BIC1Dist`` : BIC scores for the two models for model selection.

``count`` : The total number of trees in that triplet topology.

Sample scripts to parse the output are provided in the ``analysis`` folder.

## Preprint
For additional details, please see [the preprint](https://www.biorxiv.org/content/10.1101/466292v3).
