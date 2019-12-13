###############################################################
# Define the test
###############################################################

testbaseline = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT",  "TP")
  contrast.names = c("TP+", "TP-")
  
  c.1        = c(0, 1)
  c.2        = c(0, -1)
  
  cont.mat = rbind(c.1, c.2)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  
  # if you want to select certain volumes only
  X = subset(X, GROUP == 'exp' & TP %in% c(1, 2))
  
  tryCatch({ 
    model = lmer(y ~ 1 + TP + (1 |SUBJECT), data = X)
    
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    coefs = fixef(model)
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  }
  )
  
}


testlinearhalf = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT",  "TRAINING")
  contrast.names = c("TRAINING+", "TRAINING-")
  
  c.1        = c(0, 1)
  c.2        = c(0, -1)
  
  cont.mat = rbind(c.1, c.2)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  
  # if you want to select certain volumes only
  X = subset(X, GROUP == 'exp' & TP %in% c(1, 2, 3, 4))
  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  
  tryCatch({ 
    model = lmer(y ~ 1 + TRAINING + (1 |SUBJECT), data = X)
    
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    coefs = fixef(model)
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  }
  )
  
}

testgrouplinearhalf = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "GROUP", "TRAINING", "GROUP_x_TRAINING")
  contrast.names = c("GROUP_x_TRAINING+", "GROUP_x_TRAINING-")
  
  c.1        = c(0, 0, 0, 1)
  c.2        = c(0, 0, 0, -1)
  
  cont.mat = rbind(c.1, c.2)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  
  # if you want to select certain volumes only
  X = subset(X, TP %in% c(1, 2, 3, 4))
  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  
  tryCatch({ 
    model = lmer(y ~ 1 + GROUP*TRAINING + (1 |SUBJECT), data = X)
    
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    coefs = fixef(model)
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  }
  )
  
}


testlinearleft = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "GROUP", "TRAINING", "GROUP_x_TRAINING")
  contrast.names = c("GROUP_x_TRAINING+", "GROUP_x_TRAINING-")
  
  c.1        = c(0, 0, 0, 1)
  c.2        = c(0, 0, 0, -1)
  
  cont.mat = rbind(c.1, c.2)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  
  # if you want to select certain volumes only
  X = subset(X, LATERALIZATION == 'left')
  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  
  tryCatch({ 
    model = lmer(y ~ 1 + GROUP*TRAINING + (1 |SUBJECT), data = X)
    
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    coefs = fixef(model)
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  }
  )
  
}


testlinear = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "EXPERIMENT", "TRAINING")
  contrast.names = c("TRAINING+", "TRAINING-")
  
  c.1        = c(0, 0, 1)
  c.2        = c(0, 0, -1)
  
  cont.mat = rbind(c.1, c.2)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  
  # if you want to select certain volumes only
  X = subset(X, GROUP == 'exp')
  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  
  tryCatch({ 
    model = lmer(y ~ 1 + EXPERIMENT + TRAINING + (1 + TRAINING|SUBJECT), data = X)
    
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    coefs = fixef(model)
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  }
  ) 
}


testlinearFA = function(y, X, ALTERNATIVE = 'greater')
  {
    var.names = c("INTERCEPT", "GROUP", "TRAINING", "EXPERIMENT")
    contrast.names = c("GROUP+", "GROUP-", "TRAINING+" , "TRAINING-")
    
    c.1        = c(0, 1, 0, 0)
    c.2        = c(0, -1, 0, 0)
    c.3        = c(0, 0, 1, 0)
    c.4        = c(0, 0, -1, 0)
    
    cont.mat = rbind(c.1, c.2, c.3, c.4)
    colnames(cont.mat) = var.names
    rownames(cont.mat) = contrast.names
    
    tags = c(paste0(var.names, '_coef'),
             paste0(contrast.names, '_coef'),
             paste0(var.names, '_p'),
             paste0(contrast.names, '_p'))
    
    X$y = y
    levels(X$GROUP) <- c(levels(X$GROUP), 'control')
    X$GROUP[X$SUBJECT < 8] = 'control'

    tryCatch({ 
      model = lm(y ~ 1 + GROUP + TRAINING + EXPERIMENT, data = X)
      
      pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
      coefs = coefficients(model)
      glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
      contrast.pvalues = summary(glh)$test$pvalues
      contrast.coefs = coef(glh)
      val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
      names(val) = tags
      return(val)
      
    },
    
    error = function(cond){
      # something went wrong, save special values
      pvalues = rep(2, length(var.names))
      coefs = rep(0, length(var.names))
      contrast.pvalues = rep(2, length(contrast.names))
      contrast.coefs = rep(0, length(contrast.names))
      val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
      names(val) = tags
      return(val)
      
    }
    )
    
  }


testlinearcontrols = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "TRAINING")
  contrast.names = c("TRAINING+", "TRAINING-")
  
  c.1        = c(0, 1)
  c.2        = c(0, -1)
  
  cont.mat = rbind(c.1, c.2)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  
  # if you want to select certain volumes only
  X = subset(X, GROUP == 'control')
  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  
  tryCatch({ 
    model = lmer(y ~ 1 + TRAINING + (1 |SUBJECT), data = X)
    
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    coefs = fixef(model)
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  }
  ) 
}


testquadratic = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "EXPERIMENT", "TRAINING", "TRAINING.Q", "EXPERIMENT_x_TRAINING")
  contrast.names = c("TRAINING.Q+", "TRAINING.Q-")
  
  c.1        = c(0, 0, 0, 1, 0)
  c.2        = c(0, 0, 0, -1, 0)
  
  cont.mat = rbind(c.1, c.2)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  
  # if you want to select certain volumes only
  X = subset(X, GROUP == 'exp')
  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  X$TRAINING.Q = scale(X$TRAINING.Q, center = T, scale = T) 
  
  tryCatch({ 
    model = lmer(y ~ 1 + EXPERIMENT*TRAINING + TRAINING.Q + (1 + TRAINING|SUBJECT), data = X)
    
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    coefs = fixef(model)
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  }
  )
  
}

testasymptotic = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "EXPERIMENT", "TRAINING", "TRAINING.A", "EXPERIMENT_x_TRAINING")
  contrast.names = c("TRAINING.A+", "TRAINING.A-")
  
  c.1        = c(0, 0, 0, 1, 0)
  c.2        = c(0, 0, 0, -1, 0)
  
  cont.mat = rbind(c.1, c.2)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  
  # if you want to select certain volumes only
  X = subset(X, GROUP == 'exp')
  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  X$TRAINING.A = scale(X$TRAINING.A, center = T, scale = T) 
  
  tryCatch({ 
    model = lmer(y ~ 1 + EXPERIMENT*TRAINING + TRAINING.A + (1 + TRAINING|SUBJECT), data = X)
    
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    coefs = fixef(model)
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  }
  )
  
}

testlateralization.asymptotic = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "EXPERIMENT", "TRAINING", "LATERALIZATION", "TRAINING.A", 
                "EXPERIMENT_x_TRAINING", "LATERALIZATION_x_TRAINING.A")
  contrast.names = c("TRAINING.A+", "TRAINING.A-",
                     "LATERALIZATION_x_TRAINING.A+", "LATERALIZATION_x_TRAINING.A-")
  
  c.1        = c(0, 0, 0, 0, 1, 0, 0)
  c.2        = c(0, 0, 0, 0, -1, 0, 0)
  c.3        = c(0, 0, 0, 0, 0, 0, 1)
  c.4        = c(0, 0, 0, 0, 0, 0, -1)
  
  cont.mat = rbind(c.1, c.2, c.3, c.4)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  
  # if you want to select certain volumes only
  X = subset(X, GROUP == 'exp')
  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  X$TRAINING.A = scale(X$TRAINING.A, center = T, scale = T) 
  
  tryCatch({ 
    model = lmer(y ~ 1 + EXPERIMENT*TRAINING + LATERALIZATION*TRAINING.A + (1 + TRAINING|SUBJECT), data = X)
    
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    coefs = fixef(model)
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  }
  )
  
}


testlateralization.quadratic = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "EXPERIMENT", "TRAINING", "LATERALIZATION", "TRAINING.Q", 
                "EXPERIMENT_x_TRAINING", "LATERALIZATION_x_TRAINING.Q")
  contrast.names = c("LATERALIZATION_x_TRAINING.Q+", "LATERALIZATION_x_TRAINING.Q-")
  
  c.1        = c(0, 0, 0, 0, 0, 0, 1)
  c.2        = c(0, 0, 0, 0, 0, 0, -1)
  
  cont.mat = rbind(c.1, c.2)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  
  # if you want to select certain volumes only
  X = subset(X, GROUP == 'exp')
  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  X$TRAINING.Q = scale(X$TRAINING.Q, center = T, scale = T) 
  
  tryCatch({ 
    model = lmer(y ~ 1 + EXPERIMENT*TRAINING + LATERALIZATION*TRAINING.Q + (1 + TRAINING|SUBJECT), data = X)
    
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    coefs = fixef(model)
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  }
  )
  
}


testlateralization.linear = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "EXPERIMENT", "TRAINING", "LATERALIZATION",  
                "EXPERIMENT_x_TRAINING", "EXPERIMENT_x_LATERALIZATION", "TRAINING_x_LATERALIZATION",
                "EXPERIMENT_x_TRAINING_x_LATERALIZATION")
  
  contrast.names = c("LATERALIZATION_x_TRAINING+", "LATERALIZATION_x_TRAINING-")
  
  c.1        = c(0, 0, 0, 0, 0, 0, 1, 0)
  c.2        = c(0, 0, 0, 0, 0, 0, -1, 0)
  
  cont.mat = rbind(c.1, c.2)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  
  # if you want to select certain volumes only
  X = subset(X, GROUP == 'exp')
  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  
  tryCatch({ 
    model = lmer(y ~ 1 + EXPERIMENT*TRAINING*LATERALIZATION + (1 + TRAINING|SUBJECT), data = X)
    
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    coefs = fixef(model)
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  }
  )
  
}


modelcomparison = function(y, X)
{
  
  X$y = y
  
  # if you want to select certain volumes only
  X = subset(X, GROUP == 'exp')
  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  X$TRAINING.A = scale(X$TRAINING.A, center = T, scale = T) 
  X$TRAINING.Q = scale(X$TRAINING.Q, center = T, scale = T) 
  
  tryCatch({ 
    model.q = lmer(y ~ 1 + EXPERIMENT*TRAINING + TRAINING.Q + (1 + TRAINING|SUBJECT), REML = F, data = X)
    model.a = lmer(y ~ 1 + EXPERIMENT*TRAINING + TRAINING.A + (1 + TRAINING|SUBJECT), REML = F, data = X)
    
    AIC.val = AIC(model.q, model.a)
    BIC.val = BIC(model.q, model.a)    
    val = c(AIC.val$AIC[1], AIC.val$AIC[2], AIC.val$AIC[1]-AIC.val$AIC[2],
            BIC.val$BIC[1], BIC.val$BIC[2], BIC.val$BIC[1]-BIC.val$BIC[2])
    names(val) = c("AIC_quadratic", "AIC_asymptotic", "AIC_quad-asympt", "BIC_quadratic", "BIC_asymptotic", "BIC_quad-asympt")    
    return(val)
    
  },
  
  error = function(cond){
    val = rep(0, 6)
    names(val) = c("AIC_quadratic", "AIC_asymptotic", "AIC_quad-asympt", "BIC_quadratic", "BIC_asymptotic", "BIC_quad-asympt")    
    return(val)
  }
  )
  
}

testgroupasymptotic = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("INTERCEPT", "TRAINING", "GROUP", "TRAINING.A", "GROUP_x_TRAINING.A")
  contrast.names = c("GROUP_x_TRAINING.A+", "GROUP_x_TRAINING.A-")
  
  c.1        = c(0, 0, 0, 0, 1)
  c.2        = c(0, 0, 0, 0, -1)
  
  cont.mat = rbind(c.1, c.2)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  
  # if you want to select certain volumes only
  X = subset(X, EXPERIMENT == '2018')
  X$TRAINING = scale(X$TRAINING, center = T, scale = T) 
  X$TRAINING.A = scale(X$TRAINING.A, center = T, scale = T) 
  
  tryCatch({ 
    model = lmer(y ~ 1 + TRAINING + GROUP*TRAINING.A + (1 + TRAINING| SUBJECT), data = X)
    
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    coefs = fixef(model)
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)
    
  },
  
  error = function(cond){
    # something went wrong, save special values
    pvalues = rep(2, length(var.names))
    coefs = rep(0, length(var.names))
    contrast.pvalues = rep(2, length(contrast.names))
    contrast.coefs = rep(0, length(contrast.names))
    val = c(coefs, contrast.coefs, pvalues, contrast.pvalues)
    names(val) = tags
    return(val)    
  }
  )
  
}

