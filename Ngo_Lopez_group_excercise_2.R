# Gavin Lopez, and Hieu Ngo

# Your function will have the following inputs.
# 
# * x - a numeric input vector
# * y - a numeric response
#
# Note span and degree are shown with their default values. (Read about this in the description)
# * degree should be 1 or 2 only
# * span can be any value in interval (0, 1) non-inclusive.
#
# If show.plot = TRUE then you must show a plot of either the final fit

myloess <- function(x, y, span = 0.5, degree = 1, show.plot = TRUE) {
  
  N_total <- length(x)
  yhat <- numeric(N_total)
  
  # Capture variable names for axis labeling
  x_label <- deparse(substitute(x))
  y_label <- deparse(substitute(y))
  
  # Neighborhood size
  k <- ceiling(span * N_total)
  
  for (i in 1:N_total) {
    dists <- abs(x - x[i])
    
    neighbors <- order(dists)[1:k]
    x_local <- x[neighbors]
    y_local <- y[neighbors]
    
    D <- max(abs(x_local - x[i]))
    u <- abs(x_local - x[i]) / D
    w <- ifelse(u < 1, (1 - abs(u)^3)^3, 0)
    
    if (degree == 1) {
      X <- cbind(1, x_local)
      x_i_row <- c(1, x[i])
    } else if (degree == 2) {
      X <- cbind(1, x_local, x_local^2)
      x_i_row <- c(1, x[i], x[i]^2)
    }
    
    W <- diag(w)
    beta <- solve(t(X) %*% W %*% X) %*% (t(X) %*% W %*% y_local)
    
    yhat[i] <- x_i_row %*% beta
  }
  
  SSE <- sum((y - yhat)^2)
  MSE <- SSE / N_total
  
  df_plot <- data.frame(x = x, y = y, yhat = yhat)
  loessplot <- ggplot(df_plot, aes(x = x)) +
    geom_point(aes(y = y)) +
    geom_line(aes(y = yhat), color = "blue", size = 1) +
    labs(
      title = paste("LOESS - Degree:", degree, ", Span:", round(span, 2)),
      x = x_label,
      y = y_label
    )
  
  if (show.plot) {
    print(loessplot)
  }
  
  return(list(
    span = span,
    degree = degree,
    N_total = N_total,
    SSE = SSE,
    MSE = MSE,
    loessplot = loessplot
  ))
}

# Your function should return a named list containing the following:
# span: proportion of data used in each window (controls the bandwidth)
# degree: degree of polynomial
# N_total: total number of points in the data set
# SSE: Error Sum of Squares (Tells us how good of a fit we had).
# MSE: The Mean Squared Error for the Predictions
# loessplot: An object containing the ggplot so that we can see the plot later. 
#  We want this even if show.plot = FALSE
#  Note: you are NOT allowed to simply use stat_smooth() or geom_smooth() 
#         to have it automatically do LOESS.
#  You should use geom_line() or similar to plot your final the LOESS curve.

# Make sure you can access the objects properly using the $ notation.

mykNN <- function(train, test, y_train, y_test, k = 3, weighted = TRUE){
  
  # If the inputs are data frames, convert categorical predictors to dummy variables.
  if(is.data.frame(train)){
    train <- model.matrix(~ . - 1, data = train)
  }
  
  
  
  # Ensure test and train are matrices (if test is a single observation, convert it to a row matrix)
  if (is.null(dim(test))) {
    test <- matrix(test, nrow = 1)
  }
  if (is.null(dim(train))) {
    train <- matrix(train, nrow = 1)
  }
  # Standardize training data and apply the same scaling to test data.
  t_means <- colMeans(train)
  train_sds <- apply(train, 2, sd)
  train <- scale(train, center = t_means, scale = train_sds)
  test <- scale(test, center = t_means, scale = train_sds)
  
  
  n_test <- nrow(test)
  
  
  # For classification, use a character vector,
  # for regression, use a numeric vector.
  classification <- is.factor(y_train)
  if(classification){
    yhat <- vector("character", n_test)
  } else {
    yhat <- numeric(n_test)
  }
  
  # Loop over each test observation.
  for(i in 1:n_test){
    
    test_point <- test[i, ]
    # Compute Euclidean distances 
    distances <- sqrt(rowSums((train - matrix(test_point, 
                                              nrow = nrow(train), 
                                              ncol = ncol(train), 
                                              byrow = TRUE))^2))
    
    # Identify the indices of the k nearest neighbors.
    n_idx <- order(distances)[1:k]
    n_distances <- distances[n_idx]
    n_responses <- y_train[n_idx]
    
    # Classification: y_train is a factor.
    if(classification){
      if(weighted){
        # Handle any zero distances by directly using those neighbors.
        if(any(n_distances == 0)){
          exact_matches <- n_responses[n_distances == 0]
          prediction <- names(which.max(table(exact_matches)))
        } else {
          # Compute weights as inverse distances.
          weights <- 1 / n_distances
          # Sum weights for each class.
          class_weights <- tapply(weights, n_responses, sum)
          prediction <- names(which.max(class_weights))
        }
      } else {
        # Unweighted: simple majority vote.
        prediction <- names(which.max(table(n_responses)))
      }
      yhat[i] <- prediction
      
    } else {
      # Regression option
      if(weighted){
        if(any(n_distances == 0)){
          
          prediction <- mean(n_responses[n_distances == 0])
        } else {
          weights <- 1 / n_distances
          
          weights <- weights / sum(weights)
          prediction <- sum(weights * n_responses)
        }
      } else {
        
        prediction <- mean(n_responses)
      }
      yhat[i] <- prediction
    }
  }
  
  # Return the output for each type
  if(classification){
    # Convert predictions to factor with same levels as y_train.
    yhat <- factor(yhat, levels = levels(y_train))
    # Calculate accuracy and error rate.
    accuracy <- sum(yhat == y_test) / length(y_test)
    error_rate <- 1 - accuracy
    # Create confusion matrix.
    confusion_mat <- table(Predicted = yhat, Actual = y_test)
    
    return(list(yhat = yhat,
                accuracy = accuracy,
                error_rate = error_rate,
                confusion_matrix = confusion_mat,
                k = k))
  } else {
    # For regression, compute residuals and errors.
    residuals <- y_test - yhat
    SSE <- sum(residuals^2)
    MSE <- SSE / length(y_test)
    
    return(list(yhat = yhat,
                residuals = residuals,
                SSE = SSE,
                n_test = length(y_test),
                MSE = MSE,
                k = k))
  }
}