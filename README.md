# WrittenPrelim
Code for simulation study in written preliminary exam.  Below is a list of the different R scripts and a breif description of what they do.

### Report and Simulation Study
* norwood_prelim.pdf: Literature review for exam.
* norwood_simulation.pdf: Written report of simulation study. 

### R Notebook
* notebook.nb: HTML notebook with all code.
* notebook.Rmd: R Markdown file creating the notebook.

### General Scripts

* funcs.R: Functions used commonly throughout different scripts like calculating the mean reward given certain context and parameters.

* sims.R: Functions to run one simulation replicate.

* run_rims.R: Functions to run M simulation replicates, organize the results, and save the data.

* analyze_sims.R: Script to produce tables and figures of the simulation results.

### Individual Methods

* e_greedy.R: Simulate an experiment with e-greedy sampling.

* IDS_bayes.R: Simulate an experiment with Information-Directed Sampling where we calculate the information ratio via bayesian information gain.

* IDS_freq.R: Simulate an experiment with Information-Directed Sampling where we calculate the information ratio with the determinant of XtX.

* greedy.R: Simualte an experiment with greedy optimization.

* greedy_first.R: Simulate an experiment using the greedy first algorithm.

* pronzato: Simulate an experiment by choosing the action which maximizes the reward plus some information gain penalty.
