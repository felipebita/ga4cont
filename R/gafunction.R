#' Calculate GA-I
#' @description This function calculates the GA-I score of a sample.
#' @param array1 array containing the progenitor's genotypes
#' @param array2 array containing the sample's genotypes
#' @export

calculateGAI <- function(array1,array2){ return(sum(array1==array2,na.rm = T)/length(array1))}

#' Calculate GA-II
#' @description This function calculates the GA-II score of a sample.
#' @param array1 array containing the progenitor's genotypes.
#' @param array2 array containing the sample's genotypes.
#' @export

calculateGAII <- function(array1,array2){ return(sum(array1==array2 & array1==0,na.rm = T)/sum(array1==0,na.rm = T))}

#' Calculate GA-III
#' @description This function calculates the GA-III score of a sample.
#' @param ploidy population ploidy value, 4 for tetraploid and 6 for hexaploid.
#' @param array1 array containing the progenitor1's genotypes.
#' @param array2 array containing the progenitor2's genotypes.
#' @param array3 array containing the sample's genotypes
#' @export
#'
calculateGAIII <- function(ploidy,array1,array2,array3){
  if(ploidy == 4){
    count <- 0
    for(i in 1:length(array1)){
      if(is.na(array1[i]) | is.na(array2[i]) | is.na(array3[i])){ next}
      if((array1[i]==array2[i]) & (array1[i]==0)){count = count + (array3[i] == 1 | array3[i] == 2 | array3[i] == 3 | array3[i] == 4)}
      else if((array1[i]==array2[i]) & (array1[i]==1)){count = count + (array3[i] == 3 | array3[i] == 4)}
      else if((array1[i]==array2[i]) & (array1[i]==3)){count = count + (array3[i] == 0 | array3[i] == 1)}
      else if((array1[i]==array2[i]) & (array1[i]==4)){count = count + (array3[i] == 0 | array3[i] == 1 | array3[i] == 2 | array3[i] == 3)}
      else if((array1[i]==0 & array2[i]==1) | (array1[i]==1 & array2[i]==0)){count = count + (array3[i] == 2 | array3[i] == 3 | array3[i] == 4)}
      else if((array1[i]==0 & array2[i]==2) | (array1[i]==2 & array2[i]==0)){count = count + (array3[i] == 3 | array3[i] == 4)}
      else if((array1[i]==0 & array2[i]==3) | ((array1[i]==3 & array2[i]==0))){count = count + (array3[i] == 0 | array3[i] == 3 | array3[i] == 4)}
      else if((array1[i]==0 & array2[i]==4) | (array1[i]==4 & array2[i]==0)){count = count + (array3[i] == 0 | array3[i] == 1 | array3[i] == 3 | array3[i] == 4)}
      else if((array1[i]==1 & array2[i]==2) | (array1[i]==2 & array2[i]==1)){count = count + (array3[i] == 4)}
      else if((array1[i]==1 & array2[i]==3) | (array1[i]==3 & array2[i]==1)){count = count + (array3[i] == 0 | array3[i] == 4)}
      else if((array1[i]==1 & array2[i]==4) | (array1[i]==4 & array2[i]==1)){count = count + (array3[i] == 0 | array3[i] == 1 | array3[i] == 4)}
      else if((array1[i]==2 & array2[i]==3) | (array1[i]==3 & array2[i]==2)){count = count + (array3[i] == 0)}
      else if((array1[i]==2 & array2[i]==4) | (array1[i]==4 & array2[i]==2)){count = count + (array3[i] == 0 | array3[i] == 1)}
      else if((array1[i]==3 & array2[i]==4) | (array1[i]==4 & array2[i]==3)){count = count + (array3[i] == 0 | array3[i] == 1 | array3[i] == 2)}
    }
    return(count/length(array1))
  }
  else if(ploidy == 6){
    count <- 0
    for(i in 1:length(array1)){
      if(is.na(array1[i]) | is.na(array2[i]) | is.na(array3[i])){ next}
      if((array1[i]==array2[i]) & (array1[i]==0)){count = count + (array3[i] == 1 | array3[i] == 2 | array3[i] == 3 | array3[i] == 4 | array3[i] == 5 | array3[i] == 6)}
      else if((array1[i]==array2[i]) & (array1[i]==1)){count = count + (array3[i] == 3 | array3[i] == 4 | array3[i] == 5 | array3[i] == 6)}
      else if((array1[i]==array2[i]) & (array1[i]==2)){count = count + (array3[i] == 5 | array3[i] == 6)}
      else if((array1[i]==array2[i]) & (array1[i]==4)){count = count + (array3[i] == 0 | array3[i] == 1)}
      else if((array1[i]==array2[i]) & (array1[i]==5)){count = count + (array3[i] == 0 | array3[i] == 1 | array3[i] == 2 | array3[i] == 3)}
      else if((array1[i]==array2[i]) & (array1[i]==6)){count = count + (array3[i] == 0 | array3[i] == 1 | array3[i] == 2 | array3[i] == 3 | array3[i] == 4 | array3[i] == 5)}
      else if((array1[i]==0 & array2[i]==1) | (array1[i]==1 & array2[i]==0)){count = count + (array3[i] == 2 | array3[i] == 3 | array3[i] == 4 | array3[i] == 5 | array3[i] == 6)}
      else if((array1[i]==0 & array2[i]==2) | (array1[i]==2 & array2[i]==0)){count = count + (array3[i] == 3 | array3[i] == 4| array3[i] == 5| array3[i] == 6)}
      else if((array1[i]==0 & array2[i]==3) | (array1[i]==3 & array2[i]==0)){count = count + (array3[i] == 4 | array3[i] == 5 | array3[i] == 6)}
      else if((array1[i]==0 & array2[i]==4) | (array1[i]==4 & array2[i]==0)){count = count + (array3[i] == 0 | array3[i] == 4 | array3[i] == 5 | array3[i] == 6)}
      else if((array1[i]==0 & array2[i]==5) | (array1[i]==5 & array2[i]==0)){count = count + (array3[i] == 0 | array3[i] == 1 | array3[i] == 4 | array3[i] == 5 | array3[i] == 6)}
      else if((array1[i]==0 & array2[i]==6) | (array1[i]==6 & array2[i]==0)){count = count + (array3[i] == 0 | array3[i] == 1 | array3[i] == 2 | array3[i] == 4 | array3[i] == 5 | array3[i] == 6)}
      else if((array1[i]==1 & array2[i]==2) | (array1[i]==2 & array2[i]==1)){count = count + (array3[i] == 4 | array3[i] == 5 | array3[i] == 6)}
      else if((array1[i]==1 & array2[i]==3) | (array1[i]==3 & array2[i]==1)){count = count + (array3[i] == 5 | array3[i] == 6)}
      else if((array1[i]==1 & array2[i]==4) | (array1[i]==4 & array2[i]==1)){count = count + (array3[i] == 0 | array3[i] == 5 | array3[i] == 6)}
      else if((array1[i]==1 & array2[i]==5) | (array1[i]==5 & array2[i]==1)){count = count + (array3[i] == 0 | array3[i] == 1 | array3[i] == 5 | array3[i] == 6)}
      else if((array1[i]==1 & array2[i]==6) | (array1[i]==6 & array2[i]==1)){count = count + (array3[i] == 0 | array3[i] == 1 | array3[i] == 2 | array3[i] == 5 | array3[i] == 6)}
      else if((array1[i]==2 & array2[i]==3) | (array1[i]==3 & array2[i]==2)){count = count + (array3[i] == 6)}
      else if((array1[i]==2 & array2[i]==4) | (array1[i]==4 & array2[i]==2)){count = count + (array3[i] == 0 | array3[i] == 6)}
      else if((array1[i]==2 & array2[i]==5) | (array1[i]==5 & array2[i]==2)){count = count + (array3[i] == 0 | array3[i] == 1 | array3[i] == 6)}
      else if((array1[i]==2 & array2[i]==6) | (array1[i]==6 & array2[i]==2)){count = count + (array3[i] == 0 | array3[i] == 1 | array3[i] == 2 | array3[i] == 6)}
      else if((array1[i]==3 & array2[i]==4) | (array1[i]==4 & array2[i]==3)){count = count + (array3[i] == 0)}
      else if((array1[i]==3 & array2[i]==5) | (array1[i]==5 & array2[i]==3)){count = count + (array3[i] == 0 | array3[i] == 1)}
      else if((array1[i]==3 & array2[i]==6) | (array1[i]==6 & array2[i]==3)){count = count + (array3[i] == 0 | array3[i] == 1 | array3[i] == 2)}
      else if((array1[i]==4 & array2[i]==5) | (array1[i]==5 & array2[i]==4)){count = count + (array3[i] == 0 | array3[i] == 1 | array3[i] == 2)}
      else if((array1[i]==4 & array2[i]==6) | (array1[i]==6 & array2[i]==4)){count = count + (array3[i] == 0 | array3[i] == 1 | array3[i] == 2 | array3[i] == 3)}
      else if((array1[i]==5 & array2[i]==6) | (array1[i]==6 & array2[i]==5)){count = count + (array3[i] == 0 | array3[i] == 1 | array3[i] == 2 | array3[i] == 3 | array3[i] == 4)}
    }
    return(count/length(array1))
  }
}
















