Instruction on reproducible materials
----------

### Organization

- the root directory: 
  - R scripts (with suffix `.R`) implements various methods, evaluation metric, simulation studies
  - the bash script (with suffix `.sh`) is a script for conducting all numerical experiments
- the `vote` directory:
  - includes R scripts for real-data cleaning, and real data analysis

### Files description

Particularly, we target of each R script is listed below:

- `evaluation.R`: various evaluation metrics like F-norm, true positive rate, etc
- `method_implementation.R`: implements the baselines (includes RPLE, RISE, logRISE, ELASSO, RLRF)
- `simulation_main.R`: Return the result of one method on one particular dataset

- `simu_degree.R`: conduct experiments for empirical sample complexity analysis on the degree `->` reproducing Figure 1 and Table S1.
- `simu_beta.R`: conduct experiments for empirical sample complexity analysis on the "maximum" signal `->` reproducing Figure 2 and Figure S1.
- `simu_high.R`: conduct experiments for high-dimensional cases `->` reproducing Figure 3 and Figure S2
- `simu_p.R`: empirical sample complexity analysis on the dimension `->` reproducing Figure S3
- `simu_ws.R`: empirical sample complexity analysis on the weakest signal `->` reproduces Figure S4.

In the `vote` directory:
- `clean_senate.R`: to get the binary data for estimating graph
- `fit_graph.R`: to estimate graph with SLIDE
- `senate_analysis.R`: to get Figure 4 in the main text

To reproduce the results, please follow these step: 
1. download data from https://voteview.com/data. Second, run `clean_senate.R`  
2. Run `clean_senate.R`, `fit_graph.R`, `senate_analysis.R` one-by-one.