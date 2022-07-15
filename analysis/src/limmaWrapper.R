library(limma)

checkInputs <- function(measurments, targets)
{
  if(class(measurments) != "data.frame")
  {
    error_message <- paste("The measurments argument should be a data.frame. It's currently a", paste(class(measurments), ".",sep = ""))
    return(list(FALSE, error_message))
  }
  else
  {
    if(dim(measurments)[1] == 0)
    {
      error_message <- "The measurments dataframe doesn't seem to contain any measurments..."
      return(list(FALSE, error_message))
    }
    else
    {
      if(dim(measurments)[2] == 0)
      {
        error_message <- "The measurments dataframe doesn't seem to contain any samples..."
        return(list(FALSE, error_message))
      }
      else
      {
        if(class(as.matrix(measurments)[,1]) != "numeric")
        {
          return(list(FALSE, "The measurments dataframe should contain only numerical values (or NAs)."))
        }
        else
        {
          if(class(targets) != "data.frame")
          {
            error_message <- paste("The targets argument should be a data.frame. It's currently a", paste(class(targets), ".",sep = ""))
            return(list(FALSE, error_message))
          }
          else
          {
            if(dim(targets)[2] < 2)
            {
              return(list(FALSE,"The targets dataframe should have at least two columns, sample names and conditions."))
            }
            else
            {
              if(dim(targets)[1] != dim(measurments)[2])
              {
                error_message <- paste("The targets dataframe should have as many samples (targets rows) as the measurements (measurments columns). Currently, the targets dataframe has", paste(dim(targets)[1], "samples and the measurements have", paste(dim(measurments)[2],"samples.")))
                return(list(FALSE, error_message))
              }
              else
              {
                
              }
            }
          }
        }
      }
    }
  }
  return(list(TRUE, "All seems to be in order..."))
}

normalise <- function(measurments, targets, batch = NULL, method = 3)
{
  input_check <- checkInputs(measurments, targets)
  if (input_check[[1]])
  {
    
  }
  else
  {
    print(input_check[[2]])
    return(input_check[[1]])
  }
}

makeContrastsAlt <- function(targets, comparisons)
{
  cont.matrix <- matrix(0,nrow = length(unique(targets$condition)), ncol = length(comparisons))
  i <- 1
  for (comparison in comparisons)
  {
    for (j in 1:length(comparison))
    {
      cont.matrix[abs(comparison[j]),i] <- cont.matrix[abs(comparison[j]),i]+(comparison[j]/abs(comparison[j]))
    }
    i <- i + 1
  }
  return(cont.matrix)
}

poolContrasts <- function(cont.matrix, pool)
{
  max_size <- 1
  for (pooler in pool)
  {
    max_size <- ifelse(max_size > length(pooler), max_size, length(pooler))
  }
  
  max_iterations <- (max_size-1)*length(pooler)
  
  for (k in 1:max_iterations)
  {
    for (pooler in pool)
    {
      cont.matrix <- apply(cont.matrix, 2, function(x)
      {
        if (x[pooler[1]] == 1)
        {
          x[pooler[2]] <- -1
          return(x)
        }
        if (x[pooler[1]] == -1)
        {
          x[pooler[2]] <- 1
          return(x)
        }
        return(x)
      })
    }
  }
  return(cont.matrix)
}

runLimma <- function(measurements, targets, comparisons = NULL, pool = NULL, regress_out = NULL)
{
  input_check <- checkInputs(measurements, targets)
  if (input_check[[1]]) #input has correct format
  {
    if (!is.null(comparisons))
    {
      if (!is.null(regress_out))
      {
        for (regressor in regress_out)
        {
          measurements <- removeBatchEffect(measurements, targets[,regressor])
        }
      }
      
      cont.matrix <- makeContrastsAlt(targets, comparisons)
      
      if (!is.null(pool))
      {
        cont.matrix <- poolContrasts(cont.matrix, pool)
      }
      
      cont.matrix <- as.data.frame(cont.matrix)
      row.names(cont.matrix) <- unique(targets$condition)
      cont.matrix <- as.matrix(cont.matrix)
      
      fcond <- factor(targets$condition, levels = unique(targets$condition))
      
      design <- model.matrix(~0+fcond)
      design <- as.data.frame(design)
      names(design) <- unique(targets$condition)
      design <- as.matrix(design)
      
      print(cont.matrix)
      
      fit <- lmFit(measurements, design)
      fit2 <- contrasts.fit(fit, cont.matrix)
      fit2 <- eBayes(fit2)
      
      return(list(fit2, cont.matrix, fit))
    }
  }
  else
  {
    print(input_check[[2]])
    return(input_check[[1]])
  }
}


ttopFormatter <- function(ttop)
{
  ttop$ID <- row.names(ttop)
  ttop <- ttop[,c(7,1,2,3,4,5,6)]
  ttop <- ttop[complete.cases(ttop),]
  return(ttop)
}