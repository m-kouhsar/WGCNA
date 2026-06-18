gometh_dotplot <- function(gometh_res,
                           showCategory = 20, 
                           fdr.cutoff = NA, 
                           pval.cutoff = NA, 
                           label_format = 35,
                           plot.title = NA) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    # library(tibble)
    # library(stringr)  
  })
  
  plot_data <- gometh_res 
  
  if (!is.na(fdr.cutoff)) {
    plot_data <- plot_data |> filter(FDR < fdr.cutoff)
  }
  
  if (!is.na(pval.cutoff)) {
    plot_data <- plot_data |> filter(P.DE < pval.cutoff)
  }
  
  # 3. Select top categories based on GeneRatio
  plot_data <- plot_data |>
    arrange(P.DE) |>
    slice_head(n = showCategory)
  
  plot_data <- plot_data |>
    arrange(desc(GeneRatio))
  
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
    geom_point(aes_string(size = "N", fill = ifelse(is.na(fdr.cutoff),"P.DE","FDR")), shape = 21, color = "black", stroke = 0.5) +
    scale_fill_gradient(low = "#e41a1c", high = "#377eb8", name = ifelse(is.na(fdr.cutoff),"pvalue","p.adjust")) +
    scale_size_continuous(name = "Count") +
    labs(
      title = ifelse(is.na(plot.title),paste0("Enrichment Analysis: Top ",showCategory," significant results."),plot.title),
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