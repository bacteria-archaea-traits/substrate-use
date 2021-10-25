# FUNCTIONS TO CREATE PCA WITH SOME LEVEL OF GRAPHICAL CONTROL

pca_format <- function(data) {
  PC <- prcomp(data, center = TRUE, scale. = TRUE)
  data <- data.frame(obsnames = row.names(PC$x), PC$x)
  return(data)
}


pca_vectors <- function(data, var) {
  #This is the same as the pca_convert - but is required 
  #To calculate and output vectors individually from pca
  PC <- prcomp(data, center = TRUE, scale. = TRUE)
  
  x <- var[1]
  y <- var[2]

  data <- data.frame(obsnames=row.names(PC$x), PC$x)
  
  #Calculate vectors
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  return(datapc)
}

pca_prop_explained <- function(dat,var) {
  pca <- prcomp(dat, center = TRUE, scale. = TRUE)
  eigs <- pca$sdev^2
  sum_eigs <- sum(eigs)
  prop <- c()
  
  #Remove "PC" from var values
  var <- as.numeric(gsub("PC","",var))
  
  for(i in 1:length(var)) {
    if(is.numeric(var[i])) {
      this_prop <- eigs[var[i]] / sum_eigs * 100
      prop[i] <- round(this_prop, digits = 1)
    } else {
      print("ERROR: Var must be integer.")
    }
  }
  return(prop)
}