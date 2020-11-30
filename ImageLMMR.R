library(lmerTest)
library(multcomp)
library(reshape)
library(neurobase)
library(doParallel)
library(foreach)
library(dplyr)
mypackages =  c('lmerTest', 'multcomp', 'dplyr', 'reshape', 'irr')

# do shuffling properly
# add acceleration
#save num_outliers
# mass stat, use abs value?
##############################################################################
# Example

# basic inputs (see function for more options)
# IMAGING_FILE: 4-D file with VBM images
# OUTPUT_DIR: output directory
# DATA: matrix with the independent variables, same order as IMAGING_FILE
# MASK_FILE: which voxels to test
# MYTESTS: a function with the tests that should be done (see below)

# vbanalysis(
#  IMAGING_FILE, 
#   OUTPUT_DIR, 
#   DATA, 
#   MASK_FILE, 
#   MYTESTS)
###############################################################################


# define the tests as function to pass to vbanalysis
# this is just an example, variable names should correspond to those in the data matrix
# this is where lmer is called (any other test can be done: lm, etc...)

MYTESTS = function(y, X, ALTERNATIVE = 'greater')
{
  var.names = c("Intercept",  "training", "GROUPearly", "training_x_GROUPearly") #REALWEEK
  
  #  contrast.names = c("contrast_GROUP", "contrast_training_x_GROUPearly")
  contrast.names = c("contrast_training_x_GROUPearly")
  
  #  c.1        = c(0, 0, 1, 0)
  c.2        = c(0, 0, 0, 1)
  
  #  cont.mat = rbind(c.1, c.2)
  cont.mat = rbind(c.2)
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  
  # if you want to select certain volumes only
  X = subset(X, PHASE==1 & WEEK.MEAN != 2)
  
  model = lmer(y ~ 1 + t.TRAINING.1*GROUP + (1 |SUBJECT), data= X)
  
  pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
  coefs = fixef(model)
  
  glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE) 
  contrast.pvalues = summary(glh)$test$pvalues
  contrast.coefs = coef(glh)
  
  val = c(coefs, contrast.coefs, pvalues,  contrast.pvalues)
  names(val) = tags
  return(val)
  
}

###############################################################################

# functions
time.old = 0

get_time = function(set=F){
  
  now = proc.time()[3]
  time.diff = now - time.old
  time.old <<- now
  if (!set) print(paste("Duration:", format(time.diff, digits = 3, nsmall = 3), "seconds."))
}

test_test = function(do_test){
  print(summary(do_test( seq(nrow(data)), data)
  )
  )
}

markoutliersIQR = function(x){
  y = ifelse(x %in% boxplot(x, plot = F)$out, NA, x)
  return(y)
}

find_outliers = function(y, X){

  X$y = y
  if(is.null(X$DEPTH)){
  X = X %>% group_by(SUBJECT, SYSTEM) %>% dplyr::mutate(y = markoutliersIQR(y))
  } else {
  X = X %>% group_by(SUBJECT, SYSTEM, DEPTH) %>% dplyr::mutate(y = markoutliersIQR(y))
  }
  return(X$y)  
}


find_outliers_old = function(y, X){
  
  # is it a global outlier  
  #  y = imaging_data[, i]
  #  m = mean(imaging_data)
  #  s = sd(imaging_data)
  
  m = mean(y)
  s = sd(y)
  
  global_outlier = abs(y - m) > 3*s
  
  # is it a local outlier (outlier for the subject)  
  #  local_outlier = global_outlier
  #  for (subject in unique(X$SUBJECT)){
  #    z = y[X$SUBJECT == subject]
  #    m = mean(z)
  #    s = sd(z)
  #    local_outlier[X$SUBJECT == subject] = (abs(z - m) > 3*s)
  #  }
  
  #y[global_outlier] = NA
  
  z = y
  dz = y
  for (subject in unique(X$SUBJECT)){
    z[X$SUBJECT == subject] = y[X$SUBJECT == subject] - mean(y[X$SUBJECT == subject], na.rm=T)
    #    dz = c(0, diff(z))
  }
  s = sd(z, na.rm = T)
  #print(s)
  local_outlier = global_outlier
  for (subject in unique(X$SUBJECT)){
    local_outlier[X$SUBJECT == subject] = (abs(z[X$SUBJECT == subject]) > 3*s)
    #  if (sum(X$SUBJECT == subject) == 2) local_outlier[X$SUBJECT == subject] = T 
    
  }
  
  y[local_outlier] = NA
  
  return(y)
}

gifti_convert = '/home/share/Software/HCP/workbench/bin_rh_linux64/wb_command -metric-convert -from-nifti '

save_result = function(x, mask, tag, REF_FILE, to_gifti = '', flip = F){

  if (is.character(mask) & length(mask) == 1) mask = fast_readnii(mask)
  print(paste('Saving file', tag))
  n = nrow(x)
  if (is.null(n)) n = 1
  d = dim(mask)
  if (is.na(d[2])) d[2] = 1
  if (is.na(d[3])) d[3] = 1
  
  d = c(d[1:3], n)    

  data.out = nifti(array(0, d))
  data.out = copyNIfTIHeader(img = mask, data.out, drop_slots = "Dimension")
  #pixdim(data.out) <- pixdim(mask)
  #data.out@pixdim[5] = 1  
  datatype(data.out) <- 16
  bitpix(data.out) <- 32
  #data.out@xyzt_units <- mask@xyzt_units
  
  aux = mask
  if (n > 1){
    for (i in seq(n)){
      print(i/n*100)
      aux[mask > 0] <- x[i, ]
      data.out@.Data[, , , i] <- aux    
    }
  }
  else {
    data.out[mask > 0] <- x
    #data.out@qform_code <- mask@qform_code
    #data.out@sform_code <- mask@sform_code
    
  }

  #suppressWarnings(writeNIfTI(data.out, filename = tag))
  #system(paste("fslcpgeom", REF_FILE, tag))
  writenii(data.out, filename = tag)
  if (flip) system(paste("fslswapdim", tag, "-x y z", tag))
  
  if (to_gifti != ''){ 
    command = paste0(gifti_convert, tag, '.nii.gz ', to_gifti, ' ', tag, '.func.gii') 
    system(command)
    unlink(paste0(tag, '.nii.gz'))
  }
  
}  



vbanalysis  = function(IMAGING_FILE, OUTPUT_DIR, data, MASK_FILE, do_tests, excluded=NULL, to_gifti='', remove_outliers=T, 
                       contrast.names = NULL, contrast.mat = NULL, upsample = NULL, OUT_SUBDIR = NULL, alpha = 0.05, flip = F){
  
  # ensure we have an intercept
  data$Intercept = 1
  
  # open imaging data
  get_time(T)
  
  print(paste('Reading data from', IMAGING_FILE))  
  
  imaging <- fast_readnii(IMAGING_FILE)
  mask <- fast_readnii(MASK_FILE)
  
  shape <- dim(imaging)
  n_scans <- shape[4]
  n_voxels <- sum(mask > 0)
  print(paste('Voxels in mask:', n_voxels))  
  
  imaging.mat <- matrix(0, n_scans, n_voxels) 
  
  get_time()
  print('Reshaping data')  
  for (i in seq(n_scans)) {
    #print(i)
    imaging.mat[i, ] <- imaging[, , , i][mask > 0]
  }
  
  rm(imaging)
  # remove bad cases 
  if (!is.null(excluded)){
    imaging.mat = imaging.mat[-excluded, ]
    data = data[-excluded, ]
  }

  imagdim = dim(imaging.mat)
  print(paste('Final imaging data size', imagdim[1], imagdim[2]))
  
  get_time()
  
  dev = apply(imaging.mat, 2, sd)
  if(remove_outliers){
    print('Finding outliers')  
    
    cl <- makeCluster(NPROCS)
    registerDoParallel(cl)
    imaging.mat = foreach(i = seq(n_voxels), 
                          .export = c('find_outliers', 'markoutliersIQR'),
                          .packages = mypackages, 
                          .combine = 'cbind') %dopar% find_outliers(imaging.mat[, i], data)
    stopCluster(cl)
    get_time()
  }
  
  num_outliers = colSums(is.na(imaging.mat))
  print(paste('Max number of outliers =', max(num_outliers)))
  print('Testing model')  
  
  # test model
  cl <- makeCluster(NPROCS)
  registerDoParallel(cl)
  
  if (is.null(contrast.mat)){
    results = foreach(i = seq(n_voxels), 
                      .packages = mypackages,
                      .export = c('do_tests'), 
                      .combine = 'rbind') %dopar% do_tests(imaging.mat[, i], data)
  } else {
    
    results = foreach(i = seq(n_voxels), 
                      .packages = mypackages,
                      .export = c('do_tests'), 
                      .combine = 'rbind') %dopar% do_tests(imaging.mat[, i], data, 
                                                           contrast.names = contrast.names, contrast.mat = contrast.mat)
  }
  
  stopCluster(cl)
  tags = colnames(results) 
  
  #FDR adjustment
  p.where = grep( "\\_p$", tags)
  print(p.where)
  if (length(p.where) > 0){
    pvalues = results[, p.where, drop = F]
    p.fdr = apply(pvalues, 2, function(x) p.adjust(x, method="fdr"))
    colnames(p.fdr) = paste(colnames(pvalues), "fdr", sep="_")
    results = cbind(results, p.fdr)
    
    tags = colnames(results) 
    results[, grep( "\\_p$", tags)] = 1 - results[, grep( "\\_p$", tags)]
    results[, grep( "\\_fdr$", tags)] = 1 - results[, grep( "\\_fdr$", tags)]
    
    print("Num of significant voxels - uncorrected:")
    print(colSums(pvalues < alpha))
    
    print("Num of significant voxels - FDR corrected:")
    print(colSums(p.fdr < alpha))
  } else {
    tags = colnames(results)
    pvalues = NULL
  }
  
  get_time()
  
  #write out results in images
  #browser()
  if (!is.null(OUT_SUBDIR)) OUTPUT_DIR = file.path(OUTPUT_DIR, OUT_SUBDIR)
  print(paste('Writing to:', OUTPUT_DIR))
  
  if (file.exists(OUTPUT_DIR)){
    setwd(OUTPUT_DIR)
    unlink(file.path(OUTPUT_DIR, '*'))
  } else {
    dir.create(OUTPUT_DIR)
    setwd(OUTPUT_DIR)
  }  
  for (i in seq(length(tags))){
    save_result(results[, i], mask, tags[i], MASK_FILE, to_gifti, flip = flip)
    if (!is.null(upsample)){
      system(paste0('mri_vol2vol --mov ', tags[i], '.nii.gz --targ ', upsample, ' --o ', tags[i], '.nii.gz --regheader'),
             ignore.stdout = T, ignore.stderr = T)
    }
  }
  
  
  if(remove_outliers){
     save_result(num_outliers, mask, "num_outliers", MASK_FILE, to_gifti)
     save_result(dev, mask, "std", MASK_FILE, to_gifti)
     
     if (!is.null(upsample)){
         system(paste0('mri_vol2vol --mov std.nii.gz --targ ', upsample, ' --o std.nii.gz --regheader'))
          system(paste0('mri_vol2vol --mov num_outliers.nii.gz --targ ', upsample, ' --o num_outliers.nii.gz --regheader'))
     }
  }   
  print('Done!')  
  print('------------------------------------------------------------')  
  
  return(list (pvalues = pvalues, imaging.mat = imaging.mat, stats = results, mask = mask, excluded = excluded))
}


shuffle_data = function(data, by = NULL){
  # get the group for each subject
  data.shuff = data
  
  data.pattern = data %>% group_by(SUBJECT) %>% summarise(PATTERN = paste(c(TP), collapse = '-'))
  
  var.between = c("GROUP")
  var.within = c("TRAINING", "TRAINING.Q", "TRAINING.A")

  # shuffle by SUBJECT
  if ('BETWEEN' %in% by){
    for (p in unique(data.pattern$PATTERN)){
      # shuffle subjects with the same pattern
      mysubjects = data.pattern$SUBJECT[data.pattern$PATTERN == p]
      new.index = sample(length(mysubjects), replace = F)
      
      for (index in seq(length(mysubjects))){
          data.shuff[data$SUBJECT == mysubjects[index], c("SUBJECT", var.between, var.within)] = 
            data[data$SUBJECT == mysubjects[new.index[index]], c("SUBJECT", var.between, var.within)]
          
          if ('WITHIN' %in% by){
            # shuffle data within subject
            mysample = sample(seq(sum(data$SUBJECT == mysubjects[index])), replace = F)
            data.shuff[data$SUBJECT == mysubjects[index], var.within] = 
              data.shuff[data$SUBJECT == mysubjects[index], var.within][mysample, ] 
          }  
      }
    }
  }
  
  # change the sign randomly
  if ('INTERCEPT' %in% by){
    data.shuff$Intercept = 2*rbinom(nrow(data), 1, 0.5) - 1
  } 

  return(data.shuff)
}


vbanalysis_perm  = function(IMAGING_FILE, OUTPUT_DIR, data, MASK_FILE, do_tests, excluded=NULL, to_gifti='', remove_outliers=T, 
                            NPERMS=100, MINCLUSSIZE = 10, cluster.thr = 0.05, alpha.FWE = 0.05, 
                            contrast.names = NULL, contrast.mat = NULL, upsample = NULL, OUT_SUBDIR = NULL, shuffle_by = NULL,
                            statistic = 'Mass'){
  
  # ensure we have an intercept
  data$Intercept = 1
  
  # open imaging data
  
  get_time(T)
  
  print(paste('Reading data from', IMAGING_FILE))  
  #browser()
  imaging <- fast_readnii(IMAGING_FILE)
  mask <- fast_readnii(MASK_FILE)

#  mask[mask > 0] = c(rep(1, 2000), rep(0, sum(mask > 0) - 2000))
  shape <- dim(imaging)
  n_scans <- shape[4]
  n_voxels <- sum(mask > 0)
  print(paste('Voxels in mask:', n_voxels))  
  
  imaging.mat <- matrix(0, n_scans, n_voxels) 
  
  get_time()
  print('Reshaping data')  
  for (i in seq(n_scans)) {
    #print(i)
    imaging.mat[i, ] <- imaging[, , , i][mask > 0]
  }
  
  rm(imaging)
  
  # remove bad cases 
  if (!is.null(excluded)){
    imaging.mat = imaging.mat[-excluded, ]
    data = data[-excluded, ]
  }
  
  get_time()
  
  dev = apply(imaging.mat, 2, sd)
  if(remove_outliers){
    print('Finding outliers')  
    
    cl <- makeCluster(NPROCS)
    registerDoParallel(cl)
    imaging.mat = foreach(i = seq(n_voxels),
                          .packages = mypackages,
                          .export = c('find_outliers', 'markoutliersIQR'),
                          .combine = 'cbind') %dopar% find_outliers(imaging.mat[, i], data)
    stopCluster(cl)
  }
  num_outliers = colSums(is.na(imaging.mat))
  print(paste('Max number of outliers =', max(num_outliers)))
  
  
  if (!is.null(OUT_SUBDIR)) OUTPUT_DIR = file.path(OUTPUT_DIR, OUT_SUBDIR)  
  
  if (file.exists(OUTPUT_DIR)){
    setwd(OUTPUT_DIR)
    unlink(file.path(OUTPUT_DIR, '*'))
  } else {
    dir.create(OUTPUT_DIR)
    setwd(OUTPUT_DIR)
  }  
  
  print('Testing model')  
  print(paste('Analyzing', nrow(data), 'observations'))  
  
  # test model
  cl <- makeCluster(NPROCS)
  registerDoParallel(cl)

  cluster.stats = NULL
  for (perm in seq(NPERMS)){
    if (perm != 1){
      print(paste('Computing permutation ', perm, 'of', NPERMS))
      data.shuff = shuffle_data(data, by = shuffle_by)
    } else {
      data.shuff = data
    }
    
    if (is.null(contrast.mat)){
      results.perm = foreach(i = seq(n_voxels), 
                        .packages = mypackages,
                        .export = c('do_tests'), 
                        .combine = 'rbind') %dopar% do_tests(imaging.mat[, i], data.shuff)
    } else {
      
      results.perm = foreach(i = seq(n_voxels), 
                        .packages = mypackages,
                        .export = c('do_tests'), 
                        .combine = 'rbind') %dopar% do_tests(imaging.mat[, i], data.shuff, 
                                                             contrast.names = contrast.names, contrast.mat = contrast.mat)
    }
    
    # print out and find clusters
    tags = colnames(results.perm) 
    unlink('clusters.sum')
    cluster.stat = clusters = NULL
    clusters.perm = read.csv(text = "Cluster,Size,Volume,X,Y,Z,Max,index,test_name")
    
    tstats = results.perm[, grep( "\\_tstat$", tags), drop = F]
    coefs = results.perm[, grep( "\\_coef$", tags), drop = F]
    pvalues = results.perm[, grep( "\\_p$", tags), drop = F]
    ptags = colnames(pvalues)
    
    for (i in seq(length(ptags))){
      save_result(1 - pvalues[, i], mask, ptags[i], MASK_FILE, to_gifti)
      unlink('clusters.sum')
      
      if (to_gifti != ""){
        hemi = ifelse(grep("rh.", to_gifti), "rh", "lh")
        system(paste('mri_surfcluster --in', paste0(ptags[i],'.func.gii'), '--thmin', 1 - cluster.thr,  '--minarea', 
                     MINCLUSSIZE, '--hemi', hemi, '--ocn clusters.nii.gz --sum clusters.sum --subject fsaverage'))
        system('if [ `cat clusters.sum | wc -l` -eq 34 ]; then rm clusters.sum; fi') # no clusters
        
      } else { 
        system(paste('mri_volcluster --in', paste0(ptags[i],'.nii.gz'), '--thmin', 1 - cluster.thr,  '--minsizevox', 
                     MINCLUSSIZE, '--ocn clusters.nii.gz --sum clusters.sum'))
        system('if [ `cat clusters.sum | wc -l` -eq 27 ]; then rm clusters.sum; fi') # no clusters
      }
      system(paste0('cp clusters.nii.gz ', ptags[i], '_cluster.nii.gz'))                   
      
      if (file.exists('clusters.sum')) {
        clusters.sum = read.table('clusters.sum')        
        if (to_gifti != "") clusters.sum = clusters.sum[c(1, 8, 4, 5, 6, 7, 2)] 
        colnames(clusters.sum) = c("Cluster", "Size", "Volume", "X", "Y", "Z", "Max")
        clusters.sum$Mass = 0
        clusters.sum$index = i
        clusters.sum$test_name = ptags[i]
        clusters.vol = fast_readnii('clusters.nii.gz')
        clusters.mat = clusters.vol[mask > 0]
        for (j in seq(max(clusters.mat))){
          clusters.sum$Mass[j] = sum((clusters.mat == j) * (abs(tstats[, i]))) 
        }
        
        # choose statistic to use
        cluster.stat[i] = max(clusters.sum[, statistic])        
        clusters = cbind(clusters, clusters.mat)       
        clusters.perm = rbind(clusters.perm, clusters.sum)
        
      } else { 
        cluster.stat[i] = 0        
        clusters = cbind(clusters, rep(0, sum(mask)))
      }
    }
    cluster.stats = rbind(cluster.stats, cluster.stat) # maximum for each perm; rows = perms, cols = tests
    colnames(cluster.stats) = ptags
    print(cluster.stats)
    print(clusters.perm)
    get_time()
    write.table(cluster.stats, row.names = F, col.names = T, file = 'clusters.stat')
    
    if (perm == 1) {
      results = results.perm
      clusters.results = clusters.perm   
    }
    
  } # perms
  stopCluster(cl)
  
  # last iteration corresponds to unshuffled data
  #results = results.perm
  #clusters.results = clusters.perm
  
  # find pvalues and threshold
  if (nrow(clusters.results)>0){
    clusters.results$p.value = 0
    cluster_FWE = 1 + clusters*0
    
    for (i in seq(ncol(cluster.stats))){
      if (sum(clusters.results$index == i) > 0){
        clusters.results$p.value[clusters.results$index == i] = 
          sapply(clusters.results[clusters.results$index == i, statistic], 
                 function(x) sum(x <= cluster.stats[, i])/NPERMS)

        # approximation
        # fit curve to cluster.stats[, i] and derive p value
        myclusters = clusters.results[clusters.results$index == i, , drop = F]
        myclusters.sig = myclusters$p.value <= alpha.FWE
        
        # eliminate non-significant clusters from cluster image
        for ( j in seq(nrow(myclusters))){
          cluster_FWE[clusters[, i] == j, i] = myclusters$p.value[j] 
        }
        
        if (sum(myclusters.sig)> 0){
          which.clusters = myclusters$Cluster[myclusters.sig]
          clusters[, i][ ! clusters[, i] %in% which.clusters ] = 0
          print(sum(abs(clusters)))
        } else {
          clusters[, i] = 0        
        }
      }    
    } # cluster.stats
    print(clusters.results)
    write.table(clusters.results, row.names = F, col.names = T, file = 'clusters.sum')
    
    #FDR adjustment
    #pvalues = results[, grep( "\\_p$", tags)]
    p.fdr = apply(pvalues, 2, function(x) p.adjust(x, method="fdr"))
    colnames(p.fdr) = paste(colnames(pvalues), "fdr", sep="_")
    colnames(clusters) = paste(colnames(pvalues), "cluster", sep="_")
    
    cluster_maps = pvalues*(clusters > 0) + (clusters == 0)
    colnames(cluster_maps) = paste(colnames(pvalues), "cluster_map", sep="_")
    
    colnames(cluster_FWE) = paste(colnames(pvalues), "cluster_FWE", sep="_")
    
    #if (sum(myclusters.sig)> 0){    
      results = cbind(results, p.fdr, clusters, cluster_maps, cluster_FWE)
    #}
    tags = colnames(results) 
    results[, grep( "\\_p$", tags)] = 1 - results[, grep( "\\_p$", tags)]
    results[, grep( "\\_fdr$", tags)] = 1 - results[, grep( "\\_fdr$", tags)]
    results[, grep( "\\_FWE$", tags)] = 1 - results[, grep( "\\_FWE$", tags)]
    results[, grep( "\\_map$", tags)] = 1 - results[, grep( "\\_map$", tags)]
    
    #write out results
    get_time()
    print(paste('Writing to:', OUTPUT_DIR))
    
    for (i in seq(length(tags))){
      save_result(results[, i], mask, tags[i], MASK_FILE, to_gifti)
      if (!is.null(upsample)){
        system(paste0('mri_vol2vol --mov ', tags[i], '.nii.gz --targ ', upsample, ' --o ', tags[i], '.nii.gz --regheader --nearest'),
               ignore.stdout = T, ignore.stderr = T)
      }
      
    }
  } 
  
  if(remove_outliers){
     save_result(num_outliers, mask, "num_outliers", MASK_FILE, to_gifti)
     save_result(dev, mask, "std", MASK_FILE, to_gifti)
   
     if (!is.null(upsample)){
        system(paste0('mri_vol2vol --mov std.nii.gz --targ ', upsample, ' --o std.nii.gz --regheader'))
        system(paste0('mri_vol2vol --mov num_outliers.nii.gz --targ ', upsample, ' --o num_outliers.nii.gz --regheader --nearest'))
     }
  }
  print('Done!')  
  print('------------------------------------------------------------')  
  return(list (pvalues = pvalues, imaging.mat = imaging.mat, stats = results, mask = mask, excluded = excluded))
  
}


# demean
demean_imaging = function(X, data){
  for (sub in unique(data$SUBJECT)){
    X[data$SUBJECT == sub ] = 
      X[data$SUBJECT == sub] - mean(X[data$SUBJECT == sub], na.rm = T)
  }
  return(X)  
}

doPCA = function(X){
  comp = complete.cases(X)
  pca = princomp(scale(X[comp, ]))
  p = list()
  p$sdev = pca$sdev
  p$scores = X*NA
  p$scores[comp, ] = pca$scores
  p$loadings = pca$loadings
  return(list(p))
}

get_anis = function(x){
  s = x$sdev**2
  anis = (max(s) - min(s)) /sum(s)
}

get_anis_2 = function(x){
  s = x$sdev**2
  anis = (s[1] - s[2]) /sum(s)
}

save_fig = function(figname=NULL, width=6.5, height=6.5, res=600, jpg=F){
  #size in inches
  if (dev.cur()!=1) dev.off()
  if (is.null(figname)) {
    figname = paste0('plot', fig_count)
    
    fig_count <<- fig_count+1
  } 
  
  print(paste("Generating figure: ", figname))
  if (jpg){
    figname = paste0(FIGS_DIR, figname, '.jpg')
    jpeg(figname, width = floor(res/2.54*width), height = floor(res/2.54*height), pointsize=POINTSIZE)
  } else {
    figname = paste0(FIGS_DIR, figname, '.png')
    png(figname, width = floor(res/2.54*width), height = floor(res/2.54*height), pointsize=POINTSIZE)
  }
  
  
}

### NOT USED
save_index = function(x, mask, tag, REF_FILE, flip = F){
  
  if (is.character(mask) & length(mask) == 1) mask = readNIfTI(mask)
  print(paste('Saving file', tag))
  
  n = nrow(x)
  if (is.null(n)) {
    n = 1
    l = length(x)
  } else l = ncol(x)
  
  d = dim(mask)
  if (is.na(d[2])) d[2] = 1
  if (is.na(d[3])) d[3] = 1
  
  data.out = nifti(array(0, d)) 
  pixdim(data.out) <- pixdim(mask)
  data.out@pixdim[5] = 1  
  datatype(data.out) <- 16
  bitpix(data.out) <- 32
  data.out@xyzt_units <- mask@xyzt_units
  aux = mask
  indices = seq(l)
  
  # if (n > 1){
  #   for (i in seq(n)){
  #     print(i/n*100)
  #     aux[mask > 0] <- x[i, ]
  #     data.out@.Data[, , , i] <- aux    
  #   }
  # }
  #  else {
  data.out[mask > 0] <- indices
  data.out@qform_code <- mask@qform_code
  data.out@sform_code <- mask@sform_code
  
  #  }
  
  suppressWarnings(writeNIfTI(data.out, filename = tag))
  system(paste("fslcpgeom", REF_FILE, tag))
  
  if (flip) system(paste("fslswapdim", tag, "-x y z", tag))
  
  if (to_gifti != ''){ 
    command = paste0(gifti_convert, tag, '.nii.gz ', to_gifti, ' ', tag, '.func.gii') 
    system(command)
    unlink(paste0(tag, '.nii.gz'))
  }
  
}  


save_result_oro = function(x, mask, tag, REF_FILE, to_gifti = '', flip = F){
  
  if (is.character(mask) & length(mask) == 1) mask = readNIfTI(mask)
  print(paste('Saving file', tag))
  
  n = nrow(x)
  if (is.null(n)) n = 1
  d = dim(mask)
  if (is.na(d[2])) d[2] = 1
  if (is.na(d[3])) d[3] = 1
  
  d = c(d[1:3], n)    
  
  data.out = nifti(array(0, d)) 
  pixdim(data.out) <- pixdim(mask)
  data.out@pixdim[5] = 1  
  datatype(data.out) <- 16
  bitpix(data.out) <- 32
  data.out@xyzt_units <- mask@xyzt_units
  aux = mask
  
  if (n > 1){
    for (i in seq(n)){
      print(i/n*100)
      aux[mask > 0] <- x[i, ]
      data.out@.Data[, , , i] <- aux    
    }
  }
  else {
    data.out[mask > 0] <- x
    data.out@qform_code <- mask@qform_code
    data.out@sform_code <- mask@sform_code
    
  }
  
  suppressWarnings(writeNIfTI(data.out, filename = tag))
  system(paste("fslcpgeom", REF_FILE, tag))
  
  if (flip) system(paste("fslswapdim", tag, "-x y z", tag))
  
  if (to_gifti != ''){ 
    command = paste0(gifti_convert, tag, '.nii.gz ', to_gifti, ' ', tag, '.func.gii') 
    system(command)
    unlink(paste0(tag, '.nii.gz'))
  }
  
}  
