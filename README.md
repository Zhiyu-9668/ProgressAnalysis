# ProgressAnalysis
Shared script to perform analyses for the progression project

Create input file as described in the analysis plan. With the input ready, simply download the script `RunProgAnalyses.r` and run 
```
Rscript RunProgAnalyses.r YourInputFile 
```
The script requires dependents `survival`, `dplyr`, and `data.table` in R. They are all quite standard and have most likely been installed in your system, but if not, please install by running
```
install.packages(‘MissingDependent’)
```
