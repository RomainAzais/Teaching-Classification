# Analyse discriminante lineaire (LDA)
# et analyse discriminante quadratique (QDA)
# sur des donnees simulees en dimension 2

# Parametres d'un melange de gaussiennes en dimension 2 (covariances nulles)
m1 = c(-1,0)
sd1 = c(1,2)
c1 = matrix(c(sd1[1]^2,0,0,sd1[2]^2),nrow=2)
inv_c1 = solve(c1)
det_c1 = determinant(c1,logarithm=TRUE)$modulus[1]
f1 = 0.5

m2 = c(3,3)
sd2 = c(2,3)
c2 = matrix(c(sd2[1]^2,0,0,sd2[2]^2),nrow=2)
inv_c2 = solve(c2)
det_c2 = determinant(c2,logarithm=TRUE)$modulus[1]
f2 = 0.5

c = f1*c1 + f2*c2 # covariance commune (seulement pour LDA)
inv_c = solve(c)
det_c = determinant(c,logarithm=TRUE)$modulus[1]

# Classifieur bayesien (parametres connus) : QDA
exact_qda = function(x){
  d1 = as.matrix((x-m1))
  d1 = t(d1)
  s1 = -det_c1 - d1 %*% inv_c1 %*% (t(d1)) + 2*log(f1)
  
  d2 = as.matrix((x-m2))
  d2 = t(d2)
  s2 = -det_c2 - d2 %*% inv_c2 %*% (t(d2)) + 2*log(f2)
    
  crits = c(s1,s2)
  list('pred'=which.max(crits),'values'=crits)
}

# Classifieur bayesien (parametres connus) : LDA
exact_lda = function(x){
  d1 = as.matrix((x-m1))
  d1 = t(d1)
  s1 = - d1 %*% solve(c) %*% (t(d1)) + 2*log(f1)
  
  d2 = as.matrix((x-m2))
  d2 = t(d2)
  s2 = - d2 %*% solve(c) %*% (t(d2)) + 2*log(f2)

  crits = c(s1,s2)
  list('pred'=which.max(crits),'values'=crits)
}

# Visualisation de la QDA
xgrid = seq(-4,10,by=0.2)
ygrid = xgrid
par(mfrow=c(1,2))

plot(-100,-100,xlim=c(-4,10),ylim=c(-4,10),main = "Exact QDA",xlab='X1',ylab='X2')
for (x in xgrid){
  for (y in ygrid){
    points(x,y,col=exact_qda(c(x,y))$pred)
  }
}

# Visualisation de la LDA
plot(-100,-100,xlim=c(-4,10),ylim=c(-4,10),main = "Exact LDA",xlab='X1',ylab='X2')
for (x in xgrid){
  for (y in ygrid){
    points(x,y,col=exact_lda(c(x,y))$pred)
  }
}

# Simulation de donnees d'apprentissage et de test
n_train = 100
n_test = 100

n_train1 = floor(n_train*f1)
n_train2 = n_train - n_train1
n_test1 = floor(n_test*f1)
n_test2 = n_test - n_test1

x_train1 = cbind( rnorm(n_train1,mean=m1[1],sd=sd1[1]) , rnorm(n_train1,mean=m1[2],sd=sd1[2]))
x_train2 = cbind( rnorm(n_train2,mean=m2[1],sd=sd2[1]) , rnorm(n_train2,mean=m2[2],sd=sd2[2]))
y_train1 = rep(1,floor(n_train1))
y_train2 = rep(2,n_train2)

x_train = rbind(x_train1, x_train2)
y_train = c(y_train1,y_train2)

x_test1 = cbind( rnorm(n_test1,mean=m1[1],sd=sd1[1]) , rnorm(n_test1,mean=m1[2],sd=sd1[2]))
x_test2 = cbind( rnorm(n_test2,mean=m2[1],sd=sd2[1]) , rnorm(n_test2,mean=m2[2],sd=sd2[2]))
y_test1 = rep(1,n_test1)
y_test2 = rep(2,n_test2)

x_test = rbind(x_test1, x_test2)
y_test = c(y_test1,y_test2)

plot(x_train,col=y_train,xlab='X1',ylab='X2',main = "Train")
plot(x_test,col=y_test,xlab='X1',ylab='X2',main="Test")

# Estimation des parametres par maximum de vraisemblance
emv_gaussien = function(lab){
  x_train_this_lab = x_train[y_train == lab,]
  m = apply(x_train_this_lab,2,mean) # mean
  cov = t((x_train_this_lab-m))%*%(x_train_this_lab-m)/dim(x_train_this_lab)[1]
  det_cov = determinant(cov,logarithm=TRUE)$modulus[1]
  inv_cov = solve(cov)
    
  list('freq'=dim(x_train_this_lab)[1]/n_train, 'mean'=m, 'cov'=cov, 'det'=det_cov, 'inv'=inv_cov)
}

# Estimation pour la QDA
estim1 = emv_gaussien(1)
estim2 = emv_gaussien(2)

# Estimation pour la LDA : matrice de covariance commune
estims = list()
common_cov = matrix(0,nrow=2,ncol=2)
for (lab in 1:2){
  est = emv_gaussien(lab)
  estims[[as.character(lab)]] = est
  common_cov = common_cov + est$freq*est$cov
}
det_common_cov = determinant(common_cov,logarithm=TRUE)$modulus[1]
inv_common_cov = solve(common_cov)

# QDA avec parametres estimes
qda = function(x){
  crits = c()
  for (lab in 1:2){
    est = estims[[as.character(lab)]]
    f = est$freq
    m = est$mean
    det_c = est$det
    inv_c = est$inv
    
    d = as.matrix((x-m))
    d = t(d)
    s = -det_c - d %*% inv_c %*% (t(d)) + 2*log(f)
    
    crits = c(crits , s)
  }
  list('pred'=which.max(crits),'values'=crits)
}

# LDA avec parametres estimes
lda = function(x){
  crits = c()
  for (lab in 1:2){
    est = estims[[as.character(lab)]]
    f = est$freq
    m = est$mean
    
    d = as.matrix((x-m))
    d = t(d)
    s = - d %*% inv_common_cov %*% (t(d)) + 2*log(f)
    
    crits = c(crits , s)
  }
  list('pred'=which.max(crits),'values'=crits)
}

# Visualisation de la QDA (parametres estimes)
plot(-100,-100,xlim=c(-4,10),ylim=c(-4,10),main = "Estimated QDA",xlab='X1',ylab='X2')
for (x in xgrid){
  for (y in ygrid){
    points(x,y,col=qda(c(x,y))$pred)
  }
}

# Visualisation de la LDA (parametres estimes)
plot(-100,-100,xlim=c(-4,10),ylim=c(-4,10),main = "Estimated LDA",xlab='X1',ylab='X2')
for (x in xgrid){
  for (y in ygrid){
    points(x,y,col=lda(c(x,y))$pred)
  }
}

# Comparaison : parametres connus vs. parametres estimes

# QDA
plot(-100,-100,xlim=c(-4,10),ylim=c(-4,10),main = "Exact QDA",xlab='X1',ylab='X2')
for (x in xgrid){
  for (y in ygrid){
    points(x,y,col=exact_qda(c(x,y))$pred)
  }
}
plot(-100,-100,xlim=c(-4,10),ylim=c(-4,10),main = "Estimated QDA",xlab='X1',ylab='X2')
for (x in xgrid){
  for (y in ygrid){
    points(x,y,col=qda(c(x,y))$pred)
  }
}

# LDA
plot(-100,-100,xlim=c(-4,10),ylim=c(-4,10),main = "Exact LDA",xlab='X1',ylab='X2')
for (x in xgrid){
  for (y in ygrid){
    points(x,y,col=exact_lda(c(x,y))$pred)
  }
}

plot(-100,-100,xlim=c(-4,10),ylim=c(-4,10),main = "Estimated LDA",xlab='X1',ylab='X2')
for (x in xgrid){
  for (y in ygrid){
    points(x,y,col=lda(c(x,y))$pred)
  }
}

# Evaluation des erreurs de test
err_exact_qda = 0
err_exact_lda = 0
err_qda = 0
err_lda = 0

for (i in 1:n_test){
  print(i/n_test)
  
  x = x_test[i,] 
  y = y_test[i]
  
  eq = exact_qda(x)$pred
  el = exact_lda(x)$pred
  q = qda(x)$pred
  l = lda(x)$pred
  
  if (y!=eq){
    err_exact_qda = err_exact_qda + 1
  }
  if (y!=el){
    err_exact_lda = err_exact_lda + 1
  }
  if (y!=q){
    err_qda = err_qda + 1
  }
  if (y!=l){
    err_lda = err_lda + 1
  }
}
print(paste("Exact QDA : erreur de test =",err_exact_qda/n_test))
print(paste("Exact LDA : erreur de test =",err_exact_lda/n_test))
print(paste("Estimated QDA : erreur de test =",err_qda/n_test))
print(paste("Estimated LDA : erreur de test =",err_lda/n_test))
