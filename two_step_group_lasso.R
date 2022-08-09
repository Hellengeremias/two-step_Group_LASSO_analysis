# -----------------------------------------------------------------------------------------------------
# Packages
# -----------------------------------------------------------------------------------------------------

library(glmnet)
library(gglasso)
library(tidyverse)
library(parallel)
library(doParallel)

# -----------------------------------------------------------------------------------------------------
# data organization - variable selection (gglasso algorithm)
# -----------------------------------------------------------------------------------------------------

# SNP data coded as dummy variables
data <- readRDS("case_control_coord_h38_recodeA_semNA_dummies.txt")
data <- data.frame(data)

# outcome
df_y0 <- data.table::fread("als_info.csv")
y0 <- df_y0$PHENOTYPE
y0 <- ifelse(y0 == 2, 1, -1)

# aux objects for bootstrap analysis with random candidate variables 
vec_f.i <- readRDS("vec_f.txt")
vec_f.i <- vec_f.i %>% 
  mutate_if(is.factor, as.character)

vec2 <- vec_f.i %>%
  dplyr::filter(numLevels == 2) %>% 
  dplyr::select(vec_c, index)
names(vec2)[1] <- "vec"

vec <- vec_f.i %>% 
  dplyr::select(vec_c, index)
names(vec)[1] <- "vec"

# group vector for group-lasso
g2 <- (nrow(vec2))/2
g3 <- (nrow(dplyr::filter(vec_f.i, numLevels == 3)))/3
group <- c(rep(1:g2, each = 2), rep((g2 + 1):(g2 + g3), each = 3))  

# number of predictors
p <- g2 + g3

# -----------------------------------------------------------------------------------------------------
# data organization - pairwise interaction
# -----------------------------------------------------------------------------------------------------

# SNP data
X1 <- readRDS("case_control_coord_h38_recodeA_semNA.txt")
X1 <- data.frame(X1)

# outcome
y1 <- df_y0$PHENOTYPE
y1 <- as.numeric(ifelse(y1 == 2, 1, 0))

# aux objects for variable and categories names
vec.i <- readRDS("vec.txt")
vec.i <- vec.i %>% 
  mutate_if(is.factor, as.character)

vec2.i <- readRDS("vec2.txt")
vec2.i <- vec2.i %>% 
  mutate_if(is.factor, as.character)

# -----------------------------------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------
# [F1] pairwise interaction selection
# -----------------------------------------------------------------------------------------------------

procedure_int <- function(X1 = X1, y1 = y1, 
                          vec.i = vec.i, vec_f.i = vec_f.i, vec2.i = vec2.i, 
                          candidates = names_candidates, 
                          bst.sample.idx = bst.sample.idx){
  
  
  library(tidyverse)
  
  # fit
  cv.fit <- glinternet::glinternet.cv(X = as.matrix(X1[bst.sample.idx,
                                                       vec.i$index.i[vec.i$vec %in% candidates]]),
                                      Y = y1[bst.sample.idx],
                                      numLevels = vec.i$numLevels[vec.i$index.i[vec.i$vec %in% candidates]],
                                      nFolds = 10, family = "binomial")
  
  
  # lambda
  int_lambda <- data.frame(cv.fit[8])
  colnames(int_lambda)[1] <- "s0"
  row.names(int_lambda)[1] <- "lambda"
  
  # interpect
  int_b0 <- data.frame(cv.fit[[5]][[1]][1])  
  colnames(int_b0)[1] <- "s0"
  row.names(int_b0)[1] <- "b0"
  
  # main effect
  selected <- colnames(X1)[vec.i$index.i[vec.i$vec %in% candidates]]
  selected <- as.data.frame(selected)
  selected$mainEffects <- rep(1:nrow(selected))
  
  selected$coef_0 <- vector(length = nrow(selected))
  selected$coef_0 <- NA
  selected$coef_1 <- vector(length = nrow(selected))
  selected$coef_1 <- NA
  selected$coef_2 <- vector(length = nrow(selected))
  selected$coef_2 <- NA
  
  cod_selected <- data.frame(mainEffects = coef(cv.fit)$mainEffects$cat)
  cod_selected$index <- seq(1:nrow(cod_selected))
  
  selected <- dplyr::left_join(selected, cod_selected, by = "mainEffects")
  
  i = 1
  for(i in 1:nrow(selected)){
    if(i %in% cod_selected$mainEffects){
      if(selected$selected[i] %in% vec2.i$vec){
        selected$coef_0[i] <-  coef(cv.fit)$mainEffectsCoef$cat[[selected$index[i]]][1]
        selected$coef_1[i] <-  coef(cv.fit)$mainEffectsCoef$cat[[selected$index[i]]][2]
        selected$coef_2[i] <-  NA
        
      }
      if(!selected$selected[i] %in% vec2.i$vec){
        selected$coef_0[i] <-  coef(cv.fit)$mainEffectsCoef$cat[[selected$index[i]]][1]
        selected$coef_1[i] <-  coef(cv.fit)$mainEffectsCoef$cat[[selected$index[i]]][2]
        selected$coef_2[i] <-  coef(cv.fit)$mainEffectsCoef$cat[[selected$index[i]]][3]
      }
    }
  }  
  
  selected <- selected %>%
    tidyr::gather(coef_0, coef_1, coef_2, key = "coef", value = "value") %>% 
    dplyr::arrange(selected)
  
  df_coef <- data.frame(selected = vec_f.i$vec, selected_c = vec_f.i$vec_c)
  df_coef$coef <- c(rep(c("coef_0", "coef_1"), length(unique(vec2.i$index.i))), 
                        rep(c("coef_0", "coef_1", "coef_2"), 
                            length(unique(vec.i$index.i)) - length(unique(vec2.i$index.i))))
  
  selected <- dplyr::left_join(selected, df_coef,  by = c("selected", "coef"))
  
  df_coef <- dplyr::left_join(df_coef,
                              dplyr::select(selected, selected_c, value),
                              by = "selected_c")
  row.names(df_coef) <- df_coef$selected_c
  df_coef <- dplyr::select(df_coef, value)
  colnames(df_coef)[1] <- "s0"
  
  # Interaction effect
  ID_int <- coef(cv.fit)$interactions$catcat
  
  if(!is.null(ID_int)){
    l <- nrow(coef(cv.fit)$interactions$catcat)
    
    df_int0 <- as.data.frame(coef(cv.fit)$interactionsCoef$catcat[[1]])
    
    df_int <- df_int0 %>% 
      tidyr::pivot_longer(colnames(df_int0), names_to = "key_col", values_to = "coef")
    
    df_int$col_names_temp <- as.numeric(stringr::str_sub(
      stringr::str_remove(string = df_int$key_col, pattern = "cat"), end = -3))
    
    df_int$key_row <- rep(row.names(df_int0), each = ncol(df_int0)) 
    
    df_int$row_names_temp <- as.numeric(stringr::str_sub(
      stringr::str_remove(string = df_int$key_row, pattern = "cat"), end = -3))
    
    if(l > 1){
      i = 2
      for (i in 2:l){
        
        df_int_temp0 <- as.data.frame(coef(cv.fit)$interactionsCoef$catcat[[i]])
        
        df_int_temp <- df_int_temp0 %>% 
          tidyr::pivot_longer(colnames(df_int_temp0), names_to = "key_col", values_to = "coef")
        
        df_int_temp$col_names_temp <- as.numeric(stringr::str_sub(
          stringr::str_remove(string = df_int_temp$key_col, pattern = "cat"), end = -3))
        
        df_int_temp$key_row <- rep(row.names(df_int_temp0), each = ncol(df_int_temp0))
        
        df_int_temp$row_names_temp <- as.numeric(stringr::str_sub(
          stringr::str_remove(string = df_int_temp$key_row, pattern = "cat"), end = -3))
        
        df_int <- rbind(df_int, df_int_temp)
      }
    }
    
    variables_list <- selected %>% 
      dplyr::select(selected, mainEffects) %>% 
      dplyr::distinct(, .keep_all = TRUE)
    
    names(df_int)[3] <- "mainEffects"
    df_int <- dplyr::left_join(df_int, variables_list, by = "mainEffects") 
    names(df_int)[3] <- "cod_col"
    names(df_int)[6] <- "v1"
    
    names(df_int)[5] <- "mainEffects"
    df_int <- dplyr::left_join(df_int, variables_list, by = "mainEffects")
    names(df_int)[5] <- "cod_row"
    names(df_int)[7] <- "v2"
    
  }
  
  if(is.null(ID_int)){
    df_int <- "no interaction"
  }
  
  # final results
  res_int_bootstrap <- vector(mode = "list", length = 2)
  res_int_bootstrap[[1]] <- rbind(int_lambda, int_b0, df_coef)
  res_int_bootstrap[[2]] <- df_int
  
  rm(selected, df_coef)
  gc(verbose = FALSE)
  
  return(res_int_bootstrap)
}


# -----------------------------------------------------------------------------------------------------
# [F2] bootstrap analysis with random candidate variables [variable selection]
# -----------------------------------------------------------------------------------------------------

procedure <- function(X0 = data, y0 = y0, q = q, 
                      vec = vec, vec2 = vec2,
                      bst.predictor.idx = bst.predictor.idx, 
                      bst.sample.idx = bst.sample.idx,
                      X1 = X1, y1 = y1, 
                      vec.i = vec.i, vec2.i = vec2.i, vec_f.i = vec_f.i){
  
  library(tidyverse)
  
  # Estimate coefficients for each bootstrap analysis
  
  # Initialize beta
  gl_beta <- vector(length = nrow(vec))
  gl_beta[] <- NA
  names(gl_beta) <- vec$vec
  
  # index of selected predictors to create group arg on gglasso
  index <- vec$index[which(vec$index %in% bst.predictor.idx)]
  
  i = 1
  j = 1
  group <- NULL
  for(k in 1:q){
    if((index[i] %in% vec2$index) == TRUE){
      group[i] <- j
      group[i + 1] <- j
      j <- j + 1
      i <- i + 2
    }else{
      group[i] <- j
      group[i + 1] <- j
      group[i + 2] <- j
      j <- j + 1
      i <- i + 3
    }
  }
  
  # lambda estimate
  gl_lambda_cv <- gglasso::cv.gglasso(x = as.matrix(X0[bst.sample.idx, 
                                                        names(X0)[which(colnames(X0) %in% 
                                                                          vec$vec[which(vec$index %in% 
                                                                                          bst.predictor.idx)])]]),
                                      y = y0[bst.sample.idx], intercept = TRUE,
                                      group = group, 
                                      nfolds = 10, 
                                      loss = "logit",
                                      pred.loss = "loss")$lambda.1se
  
  # fit
  list_gglasso <- gglasso::gglasso(x = as.matrix(X0[bst.sample.idx, 
                                                    names(X0)[which(colnames(X0) %in% 
                                                                      vec$vec[which(vec$index %in% 
                                                                                      bst.predictor.idx)])]]), 
                                   y = y0[bst.sample.idx], group = group, intercept = TRUE, 
                                   lambda = gl_lambda_cv,
                                   loss = "logit") [c(1,2,5)] 
  
  gc(verbose = FALSE)
  
  # lambda fit
  gl_lambda <- data.frame(list_gglasso[3])
  colnames(gl_lambda)[1] <- "s0_gl"
  row.names(gl_lambda)[1] <- "lambda"
  
  # estimated interpect
  gl_b0 <- data.frame(list_gglasso[1])  
  colnames(gl_b0)[1] <- "s0_gl"
  row.names(gl_b0)[1] <- "b0"
  
  # estimated coefficients
  gl_beta[which(names(gl_beta) %in% 
                  row.names(data.frame(list_gglasso[2])))] <- list_gglasso[2]$beta
  
  gl_beta <- as.data.frame(gl_beta)
  colnames(gl_beta)[1] <- "s0_gl" 
  
  candidates <- dplyr::filter(dplyr::filter(gl_beta, s0_gl != 0), !is.na(s0_gl))

  # candidate variables for interaction model
  names_candidates <- row.names(candidates)
  names_candidates <- stringr::str_remove(names_candidates, "_0")
  names_candidates <- stringr::str_remove(names_candidates, "_1")
  names_candidates <- stringr::str_remove(names_candidates, "_2")
  names_candidates <- unique(names_candidates)
  
  # function for pairwise interaction seletcion [F1]
  matrix_int <- procedure_int(X1 = X1, y1 = y1, 
                              vec.i = vec.i, vec2.i = vec2.i, vec_f.i = vec_f.i,
                              candidates = names_candidates, 
                              bst.sample.idx = bst.sample.idx)
  
  # list for results objects
  res_complete <- vector(mode = "list", length = 3)
  
  # results from first step: variable selection
  res_complete[[1]] <- rbind(gl_lambda, gl_b0, gl_beta)
  
  # results from the second step: main and interaction effects
  res_complete[[2]] <- matrix_int[[1]]
  res_complete[[3]] <- matrix_int[[2]]
  
  rm(list_gglasso, gl_beta, matrix_int)
  gc(verbose = FALSE)
  
  return(res_complete)
}


# -----------------------------------------------------------------------------------------------------
# Random Procedure
# -----------------------------------------------------------------------------------------------------

cl <- parallel::makePSOCKcluster(4) 
doParallel::registerDoParallel(cl)

a = 1
for(a in 1:6){
  
  trials <- 12
  #system.time({
  matrix_bootstrap.1 <- foreach::foreach(icount(trials)) %dopar% {
    procedure(X0 = data, y0 = y0, q = round(.25*p, 0),
              vec = vec, vec2 = vec2,
              bst.sample.idx = sample(x = nrow(data), size = nrow(data), replace = TRUE),
              bst.predictor.idx = sample(x = p, size = round(.25*p, 0), replace = FALSE),
              X1 = X1, y1 = y1, 
              vec.i = vec.i, vec2.i = vec2.i, vec_f.i = vec_f.i)
  }
  #   })
  
  k = 1
  for(k in 1:trials){
    
    # random ID
    r <- round(runif(1, 1, 2000), 2)
    
    # variable selection from gglasso
    selected <- matrix_bootstrap.1[[k]][[1]]
    saveRDS(selected, paste("selected_", k, r, ".txt", sep = "_"))
    
    rm(selected)
    gc()
    
    # main effect from glinternet
    main_effect <- matrix_bootstrap.1[[k]][[2]]
    saveRDS(main_effect, paste("main_effect_", k, r,".txt", sep = "_"))
    
    rm(main_effect)
    gc()
    
    # interaction effect from glinternet
    int_effect <- matrix_bootstrap.1[[k]][[3]]
    saveRDS(int_effect, paste("int_effect_", k, r, ".txt", sep = "_"))
    
    rm(int_effect)
    gc()
    
  }
  
  rm(matrix_bootstrap.1)
  gc()
  
}

stopImplicitCluster()