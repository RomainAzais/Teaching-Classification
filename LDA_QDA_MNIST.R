# Analyse discriminante lineaire (LDA)
# et analyse discriminante quadratique (QDA)
# sur les donnees MNIST

y_traites = c(1,8)

# Estimation des parametres par maximum de vraisemblance
emv_gaussien = function(lab){
  x_train_this_lab = x_train[y_train == lab,]
  m = apply(x_train_this_lab,2,mean) # mean
  visu_mnist(m,lab)
  
  cov = t((x_train_this_lab-m))%*%(x_train_this_lab-m)/dim(x_train_this_lab)[1]
  det_cov = determinant(cov,logarithm=TRUE)$modulus[1]
  inv_cov = solve(cov)
  
  list('freq'=dim(x_train_this_lab)[1]/n_train, 'mean'=m, 'cov'=cov, 'det'=det_cov, 'inv'=inv_cov)
}

estimA = emv_gaussien(y_traites[1])
estimB = emv_gaussien(y_traites[2])

# Estimation pour la LDA (covariance commune)
emv_gaussien2 = function(lab){
  x_train_this_lab = x_train[y_train == lab,]
  m = apply(x_train_this_lab,2,mean) # mean
  visu_mnist(m,lab)
  
  cov = t((x_train_this_lab-m))%*%(x_train_this_lab-m)/dim(x_train_this_lab)[1]
  list('freq'=dim(x_train_this_lab)[1]/n_train, 'mean'=m, 'cov'=cov)
}

estims = list()
common_cov = matrix(0,nrow=784,ncol=784)
for (lab in y_traites){
  est = emv_gaussien2(lab)
  estims[[as.character(lab)]] = est
  common_cov = common_cov + est$freq*est$cov
}
det_common_cov = determinant(common_cov,logarithm=TRUE)$modulus[1]
inv_common_cov = solve(common_cov)

# LDA
lda = function(x){
  crits = c()
  for (lab in y_traites){
    est = estims[[as.character(lab)]]
    f = est$freq
    m = est$mean
    
    d = as.matrix((x-m))
    d = t(d)
    s = -d %*% inv_common_cov %*% (t(d)) + 2*log(f)
    
    crits = c(crits , s)
  }
  list('pred'=y_traites[which.max(crits)],'values'=crits)
}

# Exemple
x = x_test[3,]
y = y_test[3]
lda(x)
y

# Evaluation de l'erreur de test
err_lda = 0

x_testAB = x_test[y_test == y_traites,]
y_testAB = y_test[y_test == y_traites]
n_testAB = length(y_testAB)

for (i in 1:n_testAB){
  x = x_testAB[i,] 
  y = y_testAB[i]
  l = lda(x)$pred

  if (y!=l){
    err_lda = err_lda + 1
  }
  print(err_lda/i)
}
print(paste("Estimated LDA : erreur de test =",err_lda/n_testAB))
