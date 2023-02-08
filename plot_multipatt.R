# Function to automate plotting multipatt results
# Make plot of indicator correlation coeffiecients and mean abundance across all samples
# Just for taxa indicative of one group
# Can work for Phylum through Genus level
# Can work for factor with up to 10 levels
# Uses CPM for abundances

plot_multipatt <- function(mp_obj, input, tax_sum, group) {
  
  # Make results dataframe
  mp_obj_results <- mp_obj$sign %>%
    mutate(q.value = p.adjust(mp_obj$sign$p.value, method = "fdr"),
           Group = "NA",
           taxon = rownames(mp_obj$sign),
           num_groups = rowSums(mp_obj$sign[,1:(ncol(mp_obj$sign)-3)])) %>%
    filter(q.value < 0.05) %>%
    filter(num_groups == 1)
  
  # Make group variable and assign values
  for (i in 1:nrow(mp_obj_results)) {
    if (mp_obj_results[,1][i] == 1) {
      mp_obj_results$Group[i] <- substring(names(mp_obj$sign)[1], 3)
    }
  }
  for (i in 1:nrow(mp_obj_results)) {
    if (mp_obj_results[,2][i] == 1) {
      mp_obj_results$Group[i] <- substring(names(mp_obj$sign)[2], 3)
    }
  }
  if (ncol(mp_obj$sign)-3 > 2) {
    for (i in 1:nrow(mp_obj_results)) {
      if (mp_obj_results[,3][i] == 1) {
        mp_obj_results$Group[i] <- substring(names(mp_obj$sign)[3], 3)
      }
    }
  }
  if (ncol(mp_obj$sign)-3 > 3) {
    for (i in 1:nrow(mp_obj_results)) {
      if (mp_obj_results[,4][i] == 1) {
        mp_obj_results$Group[i] <- substring(names(mp_obj$sign)[4], 3)
      }
    }
  }
  if (ncol(mp_obj$sign)-3 > 4) {
    for (i in 1:nrow(mp_obj_results)) {
      if (mp_obj_results[,5][i] == 1) {
        mp_obj_results$Group[i] <- substring(names(mp_obj$sign)[5], 3)
      }
    }
  }
  if (ncol(mp_obj$sign)-3 > 5) {
    for (i in 1:nrow(mp_obj_results)) {
      if (mp_obj_results[,6][i] == 1) {
        mp_obj_results$Group[i] <- substring(names(mp_obj$sign)[6], 3)
      }
    }
  }
  if (ncol(mp_obj$sign)-3 > 6) {
    for (i in 1:nrow(mp_obj_results)) {
      if (mp_obj_results[,7][i] == 1) {
        mp_obj_results$Group[i] <- substring(names(mp_obj$sign)[7], 3)
      }
    }
  }
  if (ncol(mp_obj$sign)-3 > 7) {
    for (i in 1:nrow(mp_obj_results)) {
      if (mp_obj_results[,8][i] == 1) {
        mp_obj_results$Group[i] <- substring(names(mp_obj$sign)[8], 3)
      }
    }
  }
  if (ncol(mp_obj$sign)-3 > 8) {
    for (i in 1:nrow(mp_obj_results)) {
      if (mp_obj_results[,9][i] == 1) {
        mp_obj_results$Group[i] <- substring(names(mp_obj$sign)[9], 3)
      }
    }
  }
  if (ncol(mp_obj$sign)-3 > 9) {
    for (i in 1:nrow(mp_obj_results)) {
      if (mp_obj_results[,10][i] == 1) {
        mp_obj_results$Group[i] <- substring(names(mp_obj$sign)[10], 3)
      }
    }
  }
  
  # Check number in each group
  table(mp_obj_results$Group)
  
  # Arrange by Group then taxon, later use this to order the taxon factor for the heatmap
  mp_obj_results <- mp_obj_results %>%
    arrange(Group, taxon)
  
  # Get dataframe with indicator correlation coefficients, to plot
  mp_obj_corrs <- as.data.frame(mp_obj$str) %>%
    dplyr::select(1:length(levels(input$map_loaded$Environment))) %>%
    mutate(taxon = rownames(.)) %>%
    filter(taxon %in% mp_obj_results$taxon) %>%
    set_names(c(levels(input$map_loaded$Environment), "taxon"))
  
  # Add abundances to dataframe
  mean_CPM <- data.frame("CPM" = rowSums(tax_sum)/nrow(input$map_loaded),
                         "taxon" = rownames(tax_sum))
  mp_obj_results <- mp_obj_results %>%
    left_join(., mean_CPM, by = "taxon")
  
  # Melt dataframe for heatmap
  hm.melted <- mp_obj_corrs %>%
    melt(., id.vars = c("taxon")) %>%
    mutate(taxon = factor(taxon,
                          levels = mp_obj_results$taxon))
  
  # Plot heatmap with abundance barplot
  hm <- ggplot(data = hm.melted, 
               aes(x = factor(taxon), y = variable, fill = value)) + 
    geom_tile() + 
    scale_fill_distiller(name = "Indicator correlation index", 
                         palette = "RdBu", direction = -1, 
                         na.value = "transparent", 
                         type = "div", 
                         limits = c(-1, 1)) +
    scale_x_discrete(breaks = unique(hm.melted$taxon), 
                     labels = unique(hm.melted$taxon),
                     limits = rev(levels(hm.melted$taxon))) +
    scale_y_discrete(position = "right") +
    theme_minimal() +
    theme(legend.position = "bottom", 
          legend.direction = "horizontal",
          legend.title = element_text(size = 8), 
          legend.key.height = unit(0.4,"cm"),
          legend.key.width = unit(0.8, "cm"), 
          legend.text = element_text(size = 6)) +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
    coord_flip()
  
  # Get legend
  input.l <- get_legend(hm)
  
  # Clean
  hm.clean <- hm +
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_text(margin = margin(c(0,-5,0,0)), size = 5),
          axis.ticks.y = element_blank(), 
          axis.ticks.x = element_blank(), 
          panel.grid.major = element_blank(), 
          legend.position = "none",
          axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)), 
                                         angle = 45, hjust = 0), 
          plot.margin = margin(c(0,-2,0,5))) +
    labs(y = group)
  
  # Abundance plot
  bp.y <- ggplot(data = mp_obj_results, aes(x = taxon, y = CPM)) + 
    geom_bar(stat = "identity", fill = "grey") + 
    scale_y_continuous(position = "right") +
    scale_x_discrete(limits = rev(levels(mp_obj_results$taxon))) +
    coord_flip() + 
    theme_minimal() +
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.position = "none", 
          plot.margin = margin(c(0,5,0,0)),
          axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
    labs(y = "Mean CPM")
  
  top <- plot_grid(hm.clean, bp.y, ncol = 2, 
                   rel_widths = c(0.7, 0.3), align = "h")
  plot_grid(top, input.l, nrow = 2, rel_heights = c(0.90, 0.10))
}