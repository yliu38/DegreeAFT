#' Build an inverse-degree & Laplacian matrix from a STRING-PPI sub-network
#'
#' @param genes character vector of HGNC gene symbols present in STRING.
#'
#' @return A list with elements  
#' * **A**       – weighted adjacency (`dgCMatrix`),  
#' * **DegInv**  – inverse-degree diagonal (`dgCMatrix`),  
#' * **Lap**     – unnormalised Laplacian (`dgCMatrix`),  
#' * **genes**   – vector of genes kept after mapping.
#' @export
build_string_net <- function(genes) {
  string_db <- STRINGdb::STRINGdb$new(
                 version = "12.0", species = 9606,
                 score_threshold = 200, network_type = "full",
                 link_data = "combined_only")
  mapped   <- string_db$map(data.frame(gene = genes),
                            "gene", removeUnmappedRows = TRUE)
  id2gene  <- setNames(mapped$gene, mapped$STRING_id)

  int_raw  <- string_db$get_interactions(mapped$STRING_id)
  int      <- dplyr::filter(int_raw,
                            from %in% mapped$STRING_id,
                            to   %in% mapped$STRING_id) |>
              dplyr::transmute(from_gene = id2gene[from],
                               to_gene   = id2gene[to],
                               weight    = combined_score / 1e3)

  g        <- igraph::graph_from_data_frame(int, directed = FALSE)
  A        <- igraph::as_adjacency_matrix(g, attr = "weight",
                                          sparse = TRUE)

  k        <- pmax(Matrix::rowSums(A), .Machine$double.eps)

  list(A      = A,
       DegInv = Matrix::Diagonal(x = 1 / k),
       Lap    = Matrix::Diagonal(x = k) - A,
       genes  = rownames(A))
}

#' Construct an unsigned WGCNA adjacency + penalties
#'
#' @param expr_mat `p × n` expression matrix (genes × samples); **no NAs**.
#'
#' @return Same list structure as [build_string_net()].
#' @export
build_wgcna_net <- function(expr_mat) {
  powers <- 1:20
  sft    <- WGCNA::pickSoftThreshold(t(expr_mat),
                                     powerVector = powers,
                                     verbose = 0)
  softPwr <- sft$powerEstimate
  A <- WGCNA::adjacency(t(expr_mat), power = softPwr,
                        type = "unsigned",
                        corFnc = "cor",
                        corOptions = "use = 'p'")
  A <- as(Matrix::Matrix(A, sparse = TRUE), "dgCMatrix")
  k <- pmax(Matrix::rowSums(A), .Machine$double.eps)

  list(A      = A,
       DegInv = Matrix::Diagonal(x = 1 / k),
       Lap    = Matrix::Diagonal(x = k) - A,
       genes  = rownames(expr_mat))
}

#' Is a selected gene isolated in the underlying graph?
#'
#' @param sel_genes character vector of selected genes.
#' @param A_mat     sparse adjacency matrix (`dgCMatrix`).
#'
#' @return character vector `"isolated"` / `"connected"` matching `sel_genes`.
#' @export
isol_status <- function(sel_genes, A_mat) {
  deg <- Matrix::rowSums(A_mat != 0)
  names(deg) <- rownames(A_mat)
  ifelse(deg[sel_genes] == 0, "isolated", "connected")
}

