#!/bin/bash
# this script runs the R scripts to generate the transposon figures
# USE: generate_figures.sh  (in R_scripts dir)


############## TE Information ##############
echo "Generating TE Plots..."
# plot contradictory calls and histograms of their differences in read support
Rscript contradictory_calls.R
# plot TE counts vs coverage and plot depth of coverage per strain
Rscript coverage.R
# plot total absences, insertions, and references per transposon family 
Rscript family_freq.R
# plot histogram showing the genomic features in which TEs are located
Rscript gene_interrupt.R
# plot histogram  of the genomic distrubition of transposons (totals across all samples)--change to not totals
Rscript TE_density.R
# plots no. transposon events vs each other for each possible pairing and plot total transposons vs strain per insertions, references, and absences
Rscript Te_totals_distribution.R
# generate site frequency spectrum
Rscript allele_freq.R
#  plots total transposons vs strain per insertions, references, and absences per class and plot histogram of transposons per strain
Rscript TE_vs_DR.R


############## GWAS Mappings ##############
echo "Generating GWAS Plots..."
# plot genomic locations of QTL:
Rscript aggregate_GWAS.R
# for count triats, output new dataframe for QTL that don't meet following requirements: for QTL 100 SNPs away from TE positions and with phenotypes with non-equal medians
Rscript away_counts.R
# for position traits, output new dataframes for QTL 100 SNPs away from TE positions and with phenotypes with non-equal medians
Rscript away.R
# calculate correlations betweeen count, fraction, and activity traits
#Rscript correlations.R
# plot heatmap of correlation between counts and activity traits
#Rscript heatmap.R


############## Simulation Figures ##############
#echo "Generating Simulation Plots..."
# for simulatuon data, plot contradictory calls and histograms of their differences in read support
#Rscript contradictory_calls_sims.R
# plot the TPR and FDR vs the distance cutoff, read support threshold, and popFreq support threshold
#Rscript graph_BEDCOMPARE_MEAN_rounds20-23.R
# plot the TPR and FDR vs the read support threshold for the TEMP absence caller and the TELOCATE reference caller
#Rscript graph_BEDCOMPARE_MEAN_RSV_SIMs.R
# plot the FD, FN, and TP count vs Depth of Coverage, the FD and TP count vs Read Support while displaying the Read Support Threshold and the Read Support vs Coverage for refernce, absence, and insertion calls
#Rscript TPFD_ref.R
#Rscript TPFD_abs.R
#Rscript TPFD_ins.R

############## Rmd ##############
#echo "Running Markdown Scripts..."
#GWAS2.Rmd
#peak_TABLE.Rmd
#cer1_mappings.Rmd

echo "Done"