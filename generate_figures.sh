#!/bin/bash
# this script runs the R scripts to generate the transposon figures
# USE: generate_figures.sh  (in R_scripts dir)

# replace numbers with trait names
Rscript process_new_BF.R #replaces rename.R
# pull out main fig  traits
Rscript subset.R

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
# compare family counts for 15 strains with and without optical duplicates removed
#Rscript OD_comparison.R

############## GWAS Mappings ##############
echo "Generating GWAS Plots..."
# plot histogram of percentage of NAs per TE position
Rscript NA_pos.R
# for count traits, output new dataframe for QTL that have   phenotypes with non-equal medians
Rscript away_counts.R
# plot genomic locations of QTL:
Rscript aggregate_GWAS.R
# plot table of QTL peaks
Rscript peak_TABLE.R
# generates main figure of select GWAS mappings
Rscript main_GWAS_plots.R
# generate fine mapping results
Rscript fine_mappings.R
# snpeff on homologs and genes of lethal RNA phenotype insertions
#REMOVE Rscript homologs.R
# snpeff on lit genes and outlier strains
Rscript lit_genes.R
# snpeff on piRNA genes and outlier strains
XRscript piRNA_genes.R
# check traits of interest for LD between peaks
# REMOVE Rscript LD_check.R


############## Simulation Figures ##############
echo "Generating Simulation Plots..."
# plot TPR comparison of the 3 TE detection programs
Rscript allTPRFDR.R
# plot the TPR and FDR vs the distance cutoff, read support threshold, and popFreq support threshold
Rscript graph_BEDCOMPARE_MEAN_rounds20-23.R
# plot the TPR and FDR vs the read support threshold for the TEMP absence caller and the TELOCATE reference caller
Rscript graph_BEDCOMPARE_MEAN_RSV_SIMs.R
# plot the FD, FN, and TP count vs Depth of Coverage, the FD and TP count vs Read Support while displaying the Read Support Threshold and the Read Support vs Coverage for refernce, absence, and insertion calls
Rscript TPFD_ref.R
Rscript TPFD_abs.R
Rscript TPFD_ins.R

############## Rmd ##############
echo "Run the markdown scripts!"
X# GWAS2.Rmd / TEMP.Rmd
# GWAS2_totals.Rmd
X# IRA_VS.Rmd
X# NAs_per_strain.Rmd

echo "Done"