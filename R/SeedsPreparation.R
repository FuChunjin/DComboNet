
# essentialgene_maker = function(drugseeds, cell){
#
#   genepair = read.csv('genepairs.csv',header = T,stringsAsFactors = F)
#   drugtarget = read.csv(paste0(outputdir,'oneil/drugtarget_oneil.csv'), sep = ',',header = T, stringsAsFactors = F)[-2]
#   dt_lib = read.csv('data/drugtarget.csv')
#   dt = unique(rbind(drugtarget,dt_lib))
#   dt$Drug = toupper(dt$Drug)
#   #dt[dt$Drug == drugseeds,]$Target
#   #dt[which(dt$Drug %in% candidate),]$Target
#   drugseed_target = dt[dt$Drug %in% drugseeds,]$Target
#   if(length(genepair[genepair$drugseed_target %in% drugseed_target, ]$candidate_target) != 0 ){
#     essentialgene = unique(c(genepair[genepair$drugseed_target %in% drugseed_target, ]$candidate_target))#,
#                       #'BRAF','TRRAP','C4orf22','OR10C1','CDKN2A','MTOR','GABRG1','PRPS1L1','DDI1','CD1D','TTN','ZNF804A','STAC','RETNLB','GPR98','SLC27A6','PLEC','SETX','PIK3C2G','OR6C70','SLC26A10'))
#   }else{
#     essentialgene =c('BRAF','TRRAP','C4orf22','OR10C1','CDKN2A','MTOR','GABRG1','PRPS1L1','DDI1','CD1D','TTN','ZNF804A','STAC','RETNLB','GPR98','SLC27A6','PLEC','SETX','PIK3C2G','OR6C70','SLC26A10')
#   }
#   #essentialgene = unique(c(essentialgene, cancer))
#
#   return(essentialgene)
# }



#' Generate drugseed-score file
#' @description The function \code{seedscore} is to generate the seed-score file
#' @param seeds strings, randow walk with restart process will start from the
#'   provided drug seed which should be included in drug network
#' @param eta numeric parameter to controls the probability of restarting in the
#'   corresponding network
#' @return \code{seeds_score} drugseed-score \code{data.frame}
#' @export
#' @examples
#' seeds = "Sorafenib"
#' seedscore(seeds = seeds)

seedscore <- function(seeds,
                      eta = 1){

  drugs = seeds[1]
  genes = seeds[-1]
  if(length(drugs)!= 0 && length(genes) ==0){
    seeds_score <- data.frame(Seeds_ID = seeds,
                              Score = eta/length(drugs))
  }else if(length(drugs)!= 0 && length(genes) !=0){
    drugseed_score = eta/length(drugs)
    geneseed_score = rep((1-eta)/length(genes),length(genes))
    seeds_score <- data.frame(Seeds_ID = seeds,
                              Score = c(drugseed_score,geneseed_score),
                              stringsAsFactors = FALSE)
  }
  return(seeds_score)
}


