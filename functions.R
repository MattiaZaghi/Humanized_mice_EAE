# Function with control cluster centering (OL_1 as baseline)
calculate_signed_score <- function(seurat_obj, geneset_up, geneset_down, 
                                   control_cluster = "OL_1",
                                   cluster_col = "MOL_fine_clustering",
                                   assay = "RNA") {
  
  # Step 1: Get normalized counts
  counts <- GetAssayData(seurat_obj, assay = assay, slot = "data")
  
  # Step 2: Log2 transform: Log2(normCount + 1)
  log_expr <- log2(counts + 1)
  
  # Step 3: Find upregulated genes present in data
  up_genes_present <- intersect(geneset_up, rownames(log_expr))
  
  if (length(up_genes_present) == 0) {
    warning("No upregulated genes found in data!")
    up_score <- rep(0, ncol(log_expr))
  } else {
    cat(sprintf("  Up-regulated: %d/%d genes found\n", 
                length(up_genes_present), length(geneset_up)))
    # Calculate average expression of upregulated genes
    up_score <- Matrix::colMeans(log_expr[up_genes_present, , drop = FALSE])
  }
  
  # Step 4: Find downregulated genes present in data
  down_genes_present <- intersect(geneset_down, rownames(log_expr))
  
  if (length(down_genes_present) == 0) {
    warning("No downregulated genes found in data!")
    down_score <- rep(0, ncol(log_expr))
  } else {
    cat(sprintf("  Down-regulated: %d/%d genes found\n", 
                length(down_genes_present), length(geneset_down)))
    # Calculate average expression of downregulated genes
    down_score <- Matrix::colMeans(log_expr[down_genes_present, , drop = FALSE])
  }
  
  # Step 5: Compute signed average: up genes (+1) minus down genes (-1)
  signed_score <- up_score - down_score
  
  # Step 6: Control-center by subtracting the CONTROL CLUSTER average (OL_1)
  # Get cells in control cluster
  control_cells <- colnames(seurat_obj)[seurat_obj@meta.data[[cluster_col]] == control_cluster]
  
  if (length(control_cells) == 0) {
    warning(paste("No cells found in control cluster:", control_cluster))
    control_mean <- mean(signed_score)
  } else {
    # Calculate mean score in control cluster
    control_mean <- mean(signed_score[control_cells])
    cat(sprintf("  Control cluster (%s) mean: %.4f\n", control_cluster, control_mean))
  }
  
  # Center by subtracting control cluster mean
  centered_score <- signed_score - control_mean
  
  return(centered_score)
}


# Function that handles ENSEMBL suffixes (e.g., "MBP" matches "MBP-ENSG00001234")
calculate_signed_score_with_ensembl <- function(seurat_obj, geneset_up, geneset_down, 
                                                control_cluster = "OL_1",
                                                cluster_col = "type_fine",
                                                assay = "RNA") {
  
  # Step 1: Get normalized counts
  counts <- GetAssayData(seurat_obj, assay = assay, slot = "data")
  
  # Step 2: Log2 transform: Log2(normCount + 1)
  log_expr <- log2(counts + 1)
  
  # Step 3: Extract gene symbols from rownames (remove ENSEMBL suffix)
  all_genes_in_data <- rownames(log_expr)
  gene_symbols_in_data <- sapply(strsplit(all_genes_in_data, "-"), `[`, 1)
  
  # Step 4: COLLAPSE genes with same symbol by summing their expression
  # Create mapping: symbol -> list of ENSEMBL versions
  symbol_to_ensembl <- split(all_genes_in_data, gene_symbols_in_data)
  
  # Helper function to get collapsed expression for a gene set
  get_collapsed_expression <- function(gene_set) {
    # Find all ENSEMBL versions for each gene in the set
    ensembl_genes <- unlist(symbol_to_ensembl[gene_set])
    ensembl_genes <- ensembl_genes[!is.na(ensembl_genes)]
    
    if (length(ensembl_genes) == 0) {
      return(rep(0, ncol(log_expr)))
    }
    
    # Get expression matrix for these genes
    expr_subset <- log_expr[ensembl_genes, , drop = FALSE]
    
    # For genes with multiple ENSEMBL versions, sum them
    # Group by original symbol
    gene_groups <- sapply(strsplit(rownames(expr_subset), "-"), `[`, 1)
    
    # Sum expression per symbol, then average across symbols
    collapsed_per_symbol <- rowsum(as.matrix(expr_subset), group = gene_groups)
    mean_score <- Matrix::colMeans(collapsed_per_symbol)
    
    return(mean_score)
  }
  
  # Step 5: Calculate scores for upregulated and downregulated genes
  up_score <- get_collapsed_expression(geneset_up)
  down_score <- get_collapsed_expression(geneset_down)
  
  # Step 6: Compute signed average
  signed_score <- up_score - down_score
  
  # Step 7: Control-center by subtracting the CONTROL CLUSTER average
  control_cells <- colnames(seurat_obj)[seurat_obj@meta.data[[cluster_col]] == control_cluster]
  
  if (length(control_cells) == 0) {
    warning(paste("No cells found in control cluster:", control_cluster))
    control_mean <- mean(signed_score)
  } else {
    control_mean <- mean(signed_score[control_cells])
    cat(sprintf("  Control cluster (%s) mean: %.4f\n", control_cluster, control_mean))
  }
  
  # Center by subtracting control cluster mean
  centered_score <- signed_score - control_mean
  
  return(centered_score)
}
