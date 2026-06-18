gometh_dotplot <- function(gometh_res,
                           showCategory = 20, 
                           fdr.cutoff = 0.05, 
                           pval.cutoff = NA, 
                           label_format = 30) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tibble)
    library(stringr)  
  })

  # 1. Pre-process and calculate GeneRatio FIRST (Correct Way)
  processed_data <- gometh_res |>
    rownames_to_column(var = "ID") 
  
  # 2. Apply filters to the pathways afterward
  if (!is.null(fdr.cutoff)) {
    processed_data <- processed_data |> filter(FDR < fdr.cutoff)
  }
  
  if (!is.null(pval.cutoff)) {
    processed_data <- processed_data |> filter(P.DE < pval.cutoff)
  }
  
  # 3. Select top categories based on GeneRatio
  plot_data <- processed_data |>
    arrange(desc(GeneRatio)) |>
    slice_head(n = showCategory)
  
  if (nrow(plot_data) == 0) {
    stop("No pathways passed the specified significance cutoffs.")
  }
  
  # 4. Format Labels (Wrap long text strings)
  plot_data <- plot_data |>
    mutate(WrappedTERM = stringr::str_wrap(TERM, width = label_format))
  
  plot_data$WrappedTERM <- factor(plot_data$WrappedTERM, levels = rev(unique(plot_data$WrappedTERM)))
  
  # 5. Build the ggplot
  p <- ggplot(plot_data, aes(x = GeneRatio, y = WrappedTERM)) +
    theme_bw(base_size = 12) +
    geom_point(aes(size = N, fill = FDR), shape = 21, color = "black", stroke = 0.5) +
    scale_fill_gradient(low = "#e41a1c", high = "#377eb8", name = "p.adjust") +
    scale_size_continuous(name = "Count") +
    labs(
      title = "dotplot for ORA",
      x = "GeneRatio",
      y = NULL
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      panel.grid.major = element_line(color = "gray92"),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(color = "black", lineheight = 0.8), 
      axis.text.x = element_text(color = "black")
    )
  
  return(p)
}