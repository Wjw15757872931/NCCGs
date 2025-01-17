runStemness <- function(X, stem.sig = NULL, species = "human"){
  message("[", Sys.time(), "] -----: stemness score calculation")
  if(is.null(stem.sig)){
    stem.sig <- read.delim("./pcbc-stemsig.tsv", header = FALSE, row.names = 1)
    if(species == "mouse"){
      sig.genes <- getMouseGene(rownames(stem.sig), bool.name = T)
      stem.sig <- stem.sig[names(sig.genes), , drop=F]
      rownames(stem.sig) <- sig.genes
    }
  }
  
  common.genes <- intersect(rownames(stem.sig), rownames(X))
  X <- X[common.genes, ]
  stem.sig <- stem.sig[common.genes, ]
  
  s <- apply(X, 2, function(z) {cor(z, stem.sig, method = "sp", use = "complete.obs")})
  names(s) <- colnames(X)
  
  s <- s - min(s)
  s <- s / max(s)
  
  return(s)
}