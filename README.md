# TSAR
Time-Series Analysis tool for Respiratory Viral DREAM Challenge

## Pipeline Overview
(assuming input data is a series of .CEL & .CDF files for each sample)  
1. Quality assessment using [affy R package](https://bioconductor.org/packages/release/bioc/manuals/affy/man/affy.pdf)  
2. Log transform intensities in .CEL files  
3. Batch effect assessment as described in [my previous work](https://github.com/joshuaburkhart/StatisticalMethodsInCompBio/blob/master/HW3/hw3_writeup.pdf)  
4. Correct for batch effects using [snm R package](https://www.bioconductor.org/packages/devel/bioc/vignettes/snm/inst/doc/snm.pdf)
5. Perform Fourier transform on resulting intensities for each sample, generating frequency-domain probe values  
6. Split data into training & test sets and perform cross validation.  
7. Use CCD's [FGS algorithm](http://www.ccd.pitt.edu/wiki/index.php?title=Fast_Greedy_Search_(FGS)_Algorithm_for_Continuous_Variables) to efficiently construct a Bayesian net for the training set's frequency-domain probe values. Try to use OHSU's Exacloud cluster for this. Use a Wilcoxon rank sum method to select values if this still takes too long.  
8. Select the most informative probes using a Markov blanket  
9. Train a small collection of machine learning algorithms on the resulting values: ANN, ELM, SVN, GBM, Random Forrest.  
10. Weight each learner based on respective training set accuracy.  
11. Use probes discovered by Markov blanket for predictions.  
12. Submit predictions to Dream Challenge.  