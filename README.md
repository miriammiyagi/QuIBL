# QuIBL
Miriam Miyagi (m_miyagi@g.harvard.edu) and Nate Edelman (nedelman@g.harvard.edu)

## Quantifying Introgression via Branch Lengths

``QuIBL`` is a program for detecting introgressed loci and characterizing introgression events in a species network. It takes a set of gene trees as input and estimates the proportion of introgressed loci as well as the timing of an introgression pulse triplet by triplet.

## Installation
Currently, ``QuIBL`` is distributed as a python script. To run it, ensure that you are using Python v2.7, and have installed the following dependencies:
``joblib``, ``ete3``, ``itertools``, ``sys``, ``numpy``, ``math``, ``ConfigParser``, ``csv``, ``joblib``, and ``multiprocessing``.

If you want to use the exact versions of ``ete3`` and ``joblib`` we used, type ``pip install -r quibl_requirements.txt`` in terminal while in the main directory. For specific instructions on the Cython version of QuIBL, see the ``README`` file in that folder.

## Running QuIBL
To test ``QuIBL``, download and navigate to the repository and type ``python QuIBL.py ./Small_Test_Example/sampleInputFile.txt``. Any dependencies that may not be present can be installed using ``pip``. This will run a small provided sample data file and write an output CSV file to the same directory. To test your own data, change the ``treefile`` setting to be the path to a file with your Newick trees. Please note- these trees do not need to be quartets only, and can contain as many terminals as you like, as long as they all have the same set of terminals. For good results, you'll want to have at least a few hundred loci for the triplet topologies of interest.

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

## Interpretation of the Output

Given a species tree (((A,B),C),D); You might be interested in estimating the proportion of introgressed loci from species C into species B. In the output table, you'd see lines that might look something like this:

| triplet | outgroup | C1 | C2 | mixprop1 | mixprop2 | lambda2Dist | lambda1Dist | BIC2Dist | BIC1Dist | count|
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_B_C | A | 0 | 4.5 | 0.10 | 0.90 | 0.002 | 0.01| -1999 | -1599 | 1000 |
| A_B_C | B | 0 | 1.1 | 0.95 | 0.05 | 0.003 | 0.003| -199 | -299 | 100 |
| A_B_C | C | 0 | 6.7 | 0.01 | 0.99 | 0.004 | 0.05| -9999 | -1799 | 100000 |

To look at the putatively introgressed loci, look at the line ``| A_B_C | A |``, which have A as the outgroup. ``C2=4.5`` is the estimate of the time between the introgression event and the speciation event between A and C in coalescent units, and 90% of the loci for this topology are inferred as being likely introgressed, or 900 loci in total. You could compare the BIC values for the one distribution (ILS only) and two distribution models to see how likely this effect is to be real. In this toy example, the two distribution model has a much lower BIC value for this topology and is preferred. Introgression is not supported for loci with B as the outgroup, as the ILS proportion is high (95%) and the BIC value for the introgression model is larger. Finally, the species tree topology has huge support for the non-ILS model due to the branch in the species tree, as well as many more loci. Overall, we would estimate that 900/101100 (0.89%) of all loci or 90/1100 (81%) of loci discordant with the species tree are likely introgressed.

## Testing Species-by-Species
It's important to note that if you are interested in whether a particular species 'A' is introgressed, you only need to test each internal branch between it and the root- not every triplet. This is because many triplets will share internal branches on the species tree, and as a result their results will not be independent.

## Paper Details
For additional details, please see [the paper](https://science.sciencemag.org/content/366/6465/594), and in particular Section 7 of the supplement.
