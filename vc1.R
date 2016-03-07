vc1<-function (df, quantile_cutoff_high = 0.9, quantile_cutoff_low = 0.1) 
{
  intervals <- df %>% na.omit() %>% dplyr::distinct(CHROM, 
                                                    startPOS, endPOS) %>% dplyr::distinct(CHROM, startPOS) %>% 
    dplyr::distinct(CHROM, endPOS) %>% dplyr::arrange(CHROM, 
                                                      startPOS)
  strains <- as.character(na.omit(unique(df$strain)))
  intervalGENES <- list()
  for (i in 1:nrow(intervals)) {
    print(paste(100 * signif(i/nrow(intervals), 3), "%", 
                sep = ""))
    nstrains <- data.frame(df) %>% na.omit() %>% dplyr::filter(trait == 
                                                                 as.character(intervals[i, "trait"]))
    nstrains <- length(unique(nstrains$strain))
    chr <- as.character(intervals[i, ]$CHROM)
    left <- intervals[i, ]$startPOS
    right <- intervals[i, ]$endPOS
    region_of_interest <- paste0(chr, ":", left, "-", right)
    snpeff_output <- snpeff(region = region_of_interest,severity = c("LOW", "MODERATE", "HIGH", "MODIFIER"))
    pruned_snpeff_output <- snpeff_output %>% dplyr::filter(strain %in% 
                                                              strains) %>% dplyr::filter(!is.na(impact)) %>% dplyr::distinct(CHROM, 
                                                                                                                             POS, strain, effect, gene_id) %>% dplyr::arrange(effect) %>% 
      dplyr::select(CHROM, POS, REF, ALT, GT, effect, nt_change, 
                    aa_change, gene_name, gene_id, transcript_biotype, 
                    feature_type, strain) %>% dplyr::group_by(CHROM, 
                                                              POS, effect) %>% dplyr::filter(!is.na(GT), GT != 
                                                                                               "HET") %>% dplyr::mutate(num_allele = ifelse(GT == 
                                                                                                                                              "REF", 0, ifelse(GT == "ALT", 1, NA))) %>% dplyr::mutate(num_alt_allele = sum(num_allele, 
                                                                                                                                                                                                                            na.rm = T), num_strains = n()) %>% dplyr::filter(num_alt_allele/num_strains > 
                                                                                                                                                                                                                                                                               0.05) %>% dplyr::filter(num_strains > nstrains * 
                                                                                                                                                                                                                                                                                                         0.8) %>% dplyr::ungroup()
    if (nrow(pruned_snpeff_output) > 0) {
      interval_df <- df %>% dplyr::filter(CHROM == chr, 
                                          startPOS == left, endPOS == right) %>% dplyr::group_by(trait, 
                                                                                                 CHROM, startPOS, endPOS) %>% dplyr::filter(log10p == 
                                                                                                                                              max(log10p)) %>% dplyr::distinct(trait, startPOS, 
                                                                                                                                                                               endPOS, peakPOS, strain) %>% dplyr::ungroup() %>% 
        dplyr::select(trait, startPOS, endPOS, peakPOS, 
                      strain, log10p, pheno_value = value)
      pheno_snpeff_df <- pruned_snpeff_output %>% dplyr::left_join(., 
                                                                   interval_df, by = "strain", copy = TRUE) %>% 
        dplyr::distinct(strain, trait, pheno_value, gene_id, 
                        CHROM, POS, aa_change) %>% dplyr::group_by(trait, 
                                                                   CHROM, POS, effect, feature_type) %>% dplyr::mutate(spearman_cor = cor(pheno_value, 
                                                                                                                                          num_allele, method = "spearman", use = "pairwise.complete.obs")) %>% 
        dplyr::ungroup() %>% dplyr::mutate(abs_spearman_cor = abs(spearman_cor)) %>% 
        dplyr::filter(abs_spearman_cor > quantile(abs_spearman_cor, 
                                                  probs = quantile_cutoff_high, na.rm = T)) %>% 
        dplyr::ungroup() %>% dplyr::arrange(desc(abs_spearman_cor))
      intervalGENES[[i]] <- pheno_snpeff_df
    }
    else {
      intervalGENES[[i]] <- NA
    }
  }
  return(intervalGENES)
}
