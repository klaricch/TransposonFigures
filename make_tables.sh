#!/usr/bin/bash



#CHANGE BACK TO JSUT FIGURES
cp ../figures/Peak_Table.txt .
cp ../figures/Pseudogene_Table.txt .
cp ../figures/Severe_Table.txt .
cp ../figures/dEX_TOI_Table.txt .
cp ../figures/Chi_Table.txt .
cp ../figures/Chi_IncDec.txt .
cp ../figures/Raw_Mapping_Table.txt .
cp ../figures/TajimaD_Table.txt .
cp ../figures/Outlier_Table.txt .
cp ../figures/Outlier_family.txt .
cp ../data/pi_variants.txt .

#add in the gene info for regions with high tajima D and TEs
python /Users/kristen/Documents/transposon_figure_data/R_scripts/merge_taj.py
#fix headers for supp tabels 1 and 2
python /Users/kristen/Documents/transposon_figure_data/R_scripts/fix_supp12_headers.py
bash ../R_scripts/transpose_matrix.sh tmp_T_with_monomorphic_kin_matrix_full.txt UT_with_monomorphic_kin_matrix_full.txt
python /Users/kristen/Documents/transposon_figure_data/R_scripts/replace_CE_heavy.py UT_with_monomorphic_kin_matrix_full.txt >tmp && mv tmp UT_with_monomorphic_kin_matrix_full.txt
cat UT_with_monomorphic_kin_matrix_full.txt |sed 's/WBTransposon/WBT/g'  > tmp && mv tmp UT_with_monomorphic_kin_matrix_full.txt;
#fix te mappings names
python /Users/kristen/Documents/transposon_figure_data/R_scripts/fix_T_C_names.py
#fix piRNA table header
cat table_piRNAs.txt|sed '1d' > tmp_piRNA_table.txt
#add info to peak_table
python /Users/kristen/Documents/transposon_figure_data/R_scripts/edit_peak_table.py
mv tmp.txt Peak_Table.txt

#remove BWA and keep BLAST
cat BWA_and_BLAST_table.txt |cut -f1,4-6 > BLAST_table.txt
cat tmp_piRNA_table.txt |grep -E "Transcript|BLAST|Both"|cut -f1-3 > tmp2_piRNA_table.txt 
#mkdir supplemental
#cd supplemental
echo "renaming to supplemental..."
#RENUMERB HERE
bash ../R_scripts/transpose_matrix.sh UT_with_monomorphic_kin_matrix_full.txt T_UT_with_monomorphic_kin_matrix_full.txt
cp T_UT_with_monomorphic_kin_matrix_full.txt Supplemental_Table_S1.txt 
#cp tmp_T_with_monomorphic_kin_matrix_full.txt Supplemental_Table_S1.txt 
cat T_kin_matrix_full.txt| sed 's/Marker/Strain/' >Supplemental_Table_S2.txt 

#cp tmp_T_kin_matrix_full.txt Supplemental_Table_S2.txt 
cp Outlier_Table.txt Supplemental_Table_S3.txt 
cp TajimaD_Table.txt Supplemental_Table_S4.txt 
cp Chi_Table.txt Supplemental_Table_S5.txt 
#cp Pseudogene_Table.txt Supplemental_Table_S6.txt 
cp essentiality_nonredundant_strain_info.txt Supplemental_Table_S6.txt 
cp essentiality_nonredundant_cds_lethal.txt Supplemental_Table_S7.txt 
cp tmp_T_kin_C_matrix_full_reduced.txt Supplemental_Table_S8.txt 
cp Raw_Mapping_Table.txt Supplemental_Table_S9.txt 
cp Peak_Table.txt Supplemental_Table_S10.txt 
cp dEX_TOI_Table.txt Supplemental_Table_S11.txt 

cp BLAST_table.txt Supplemental_Table_S12.txt 
cp tmp2_piRNA_table.txt Supplemental_Table_S13.txt

cp Outlier_family.txt Supplemental_Table_S14.txt
cp pi_variants.txt Supplemental_Table_S15.txt

cp Severe_Table.txt Supplemental_Table_S16.txt 
cp WB_pos_table.txt Supplemental_Table_S17.txt 
cp TE_seqs.txt Supplemental_Table_S18.txt 

cp Chi_IncDec.txt Table2.txt
#??family outlier table???genomic regions of te insertions

# get rid of line with just just _R on line 2
#cat Supplemental_Table_S1.txt |sed '2d' > tmp && mv tmp Supplemental_Table_S1.txt
#cat Supplemental_Table_S1.txt |sed '2d' > tmp && mv tmp Supplemental_Table_S1.txt

#cat Supplemental_Table_S2.txt |sed '2d' > tmp && mv tmp Supplemental_Table_S2.txt




#cp Outlier_Table.txt Supplmental_Table_S6.txt


mkdir final_tables
mv Supplemental* final_tables/
cd final_tables

python ../../R_scripts/total_table.py


for file in /Users/kristen/Documents/transposon_figure_data/tables/final_tables/Supplemental*;do
	#echo "Processing $file";
	cat $file |sed 's/WBTransposon/WBT/g'  > tmp && mv tmp $file;
	cat $file |sed 's/(cu)/(all)/g'  > tmp && mv tmp $file;
	python /Users/kristen/Documents/transposon_figure_data/R_scripts/replace_CE.py $file >tmp && mv tmp $file
done

cp ../Table1.txt .
cp ../Table2.txt .

#echo "converting to excel format..."
#python2 ../../R_scripts/make_excel.py ../Chi_Table.txt 
#python2 ../../R_scripts/make_excel.py ../Peak_Table.txt
#python2 ../../R_scripts/make_excel.py ../Pseudogene_Table.txt
#python2 ../../R_scripts/make_excel.py ../Severe_Table.txt
#python2 ../../R_scripts/make_excel.py ../dEX_TOI_Table.txt
#python2 ../../R_scripts/make_excel.py ../Outlier_Table.txt
#python2 ../../R_scripts/make_excel.py ../Raw_Mapping_Table.txt

#echo "renaming to supplemental..."
#cp Chi_Table.xlsx Supplmental_Table_S1.xlsx 
#cp Peak_Table.xlsx Supplmental_Table_S2.xlsx
#cp Pseudogene_Table.xlsx Supplmental_Table_S3.xlsx
#cp Severe_Table.xlsx Supplmental_Table_S4.xlsx
#cp dEX_TOI_Table.xlsx Supplmental_Table_S5.xlsx
#cp Outlier_Table.xlsx Supplmental_Table_S6.xlsx
#cp Peak_Table.xlsx Supplmental_Table_S7.xlsx
#cp Raw_Mapping_Table.xlsx Supplmental_Table_S8.xlsx