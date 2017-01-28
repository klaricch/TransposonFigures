# Generate Transposon Figures and Tables

#### Set Up Directory Structure
```
mkdir data
mkdir results
mkdir figures
```
#### Plots and tables
Transfer files over from data_for_figures directory on the cluster and obtain complete_mapping_df.Rda and processed_transposons.Rda files. Generate transposon plots and tables. Run Rmd scripts separately.
```
bash generate_figures.sh
bash make_tables.sh
```
