compareBC <- function (input, sampleID, variable, level) {
  
  # Summarize taxonomy
  if (level == 2) {
    tax_sum <- summarize_taxonomy(input, level = 2, report_higher_tax = F, relative = F)
    bc <- calc_dm(tax_sum)
  }
  if (level == 3) {
    tax_sum <- summarize_taxonomy(input, level = 3, report_higher_tax = F, relative = F)
    bc <- calc_dm(tax_sum)
  }
  if (level == 4) {
    tax_sum <- summarize_taxonomy(input, level = 4, report_higher_tax = F, relative = F)
    bc <- calc_dm(tax_sum)
  }
  if (level == 5) {
    tax_sum <- summarize_taxonomy(input, level = 5, report_higher_tax = F, relative = F)
    bc <- calc_dm(tax_sum)
  }
  if (level == 6) {
    tax_sum <- summarize_taxonomy(input, level = 6, report_higher_tax = F, relative = F)
    bc <- calc_dm(tax_sum)
  }
  
  # Get comparisons
  bray_mat <- as.matrix(bc)
  bray_mat[upper.tri(bray_mat, diag = TRUE)] <- NA
  bray_df <- as.data.frame(bray_mat)
  bray_df$sampleID <- rownames(bray_df)
  bray_df_long <- melt(bray_df, id.vars = "sampleID")
  bray_df_long <- na.omit(bray_df_long)
  bray_df_long$sampleID <- as.factor(bray_df_long$sampleID)
  vars <- dplyr::select(input_arc_CPM_nz$map_loaded, sampleID, all_of(variable))
  bray_df_long <- inner_join(bray_df_long, vars, by = c("sampleID" = "sampleID"))
  bray_df_long <- inner_join(bray_df_long, vars, by = c("variable" = "sampleID"))
  column1 <- paste(variable, "x", sep = ".")
  column2 <- paste(variable, "y", sep = ".")
  for (i in 1:nrow(bray_df_long)) {
    ifelse(bray_df_long$Environment.x[i] == bray_df_long$Environment.y[i],
           bray_df_long$comparison[i] <- "within",
           bray_df_long$comparison[i] <- "between")
  }
  bray_df_long$comparison <- as.factor(bray_df_long$comparison)
  return(bray_df_long)
}