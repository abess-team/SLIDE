# Instructions for Reproducible Materials

## Organization

- **Bash script** (`batch.sh`): automate the execution of all numerical simulation studies
- **R scripts** (`.R`): implement baseline methods, evaluation metrics, simulation studies, real-data analysis of the congressional voting dataset
- **Data files** (`.csv`):  

## File Descriptions

### Main R scripts
- `simu_degree.R` — empirical sample complexity analysis with respect to the degree.  
- `simu_beta.R` — empirical sample complexity analysis with respect to the maximum signal.  
- `simu_high.R` — experiments for high-dimensional cases.  
- `simu_p.R` — empirical sample complexity analysis with respect to the dimension.  
- `simu_ws.R` — empirical sample complexity analysis with respect to the weakest signal.  
- `DataAnalysis.R` — real-data analysis: data cleaning, estimation of the graphical structure among senators, and visualization.  

#### Utility R scripts (automatically used by the main scripts)
- `simulation_main.R` — runs one method on a given simulated dataset.  
- `method_implementation.R` — implementations of baseline methods (RPLE, RISE, logRISE, ELASSO, RLRF).  
- `evaluation.R` — evaluation metrics (e.g., Frobenius norm, true positive rate).  

## Reproducing Results

The scripts reproduce the results presented in the paper as follows:  

- **Figure 1 and Table S1** → `simu_degree.R`  
- **Figure 2 and Figure S1** → `simu_beta.R`  
- **Figure 3 and Figure S2** → `simu_high.R`  
- **Figure S3** → `simu_p.R`  
- **Figure S4** → `simu_ws.R`  
- **Figure 4** → `DataAnalysis.R`  

The simplest procedure on reproduction:
  1. Use the provided bash scripts (`batch.sh`) to execute the full set of simulation automatically.  
  2. Run `DataAnalysis.R` to reproduce the real-data analysis (Figure 4). 
