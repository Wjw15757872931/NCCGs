### get go/kegg genes

source("./script/getGoTerm.R")
GO_DATA <- get_GO_data("org.Hs.eg.db", "ALL", "SYMBOL")
save(GO_DATA, file = "./Rdata/GO_DATA_human.RData")

getGO <- function(ID){
  
  if(!exists("GO_DATA"))
    load("./Rdata/GO_DATA.RData")
  allNAME = names(GO_DATA$PATHID2EXTID)
  if(ID %in% allNAME){
    geneSet = GO_DATA$PATHID2EXTID[ID]
    names(geneSet) = GO_DATA$PATHID2NAME[ID]
    return(geneSet)     
  } else{
    cat("No results!\n")
  }
}

load("./Rdata/GO_DATA_human.RData") 
getGO("GO:0006959") 

getKEGG <- function(ID){
  
  library("KEGGREST")
  
  gsList = list()
  for(xID in ID){
    
    gsInfo = keggGet(xID)[[1]]
    if(!is.null(gsInfo$GENE)){
      geneSetRaw = sapply(strsplit(gsInfo$GENE, ";"), function(x) x[1])   
      xgeneSet = list(geneSetRaw[seq(2, length(geneSetRaw), 2)])          
      NAME = sapply(strsplit(gsInfo$NAME, " - "), function(x) x[1])
      names(xgeneSet) = NAME
      gsList[NAME] = xgeneSet 
    } else{
      cat(" ", xID, "No corresponding gene set in specific database.\n")
    }
  }
  return(gsList)
}
  genelist <- as.data.frame(getKEGG('bta00020'))
  write.csv(genelist,"bta00020.csv")