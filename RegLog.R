# Regression logistique
# sur des donnees simulees en dimension 2

# Simulation du modele en dimension 2 (X gaussien)
sigma = function(z){
  1/(1+exp(-z))
}

simulation = function(n, par){
  d = length(par)-1
  X = cbind(1 , matrix(rnorm(d*n,mean=0,sd=1),nrow=n,ncol=d))
  Y = rep(-1,n)
  for (i in 1:n){
    u = runif(1)
    p = sigma(par%*%X[i,])
    if (u<p){
      Y[i] = 1
    }
  }
  list('X'=X,'Y'=Y)
}

# Classifieur bayesien
classifieur = function(x, par){
  y = -1
  if (sigma(par %*% x)>0.5){
    y = 1
  }
  y
}

# Estimation des parametres : descente de gradient
gradient = function(X, Y, par){
  s = apply(X[Y==1,],2,sum)
  for (i in 1:length(Y)){
    s = s - as.numeric(sigma(par%*%X[i,]))*X[i,]
  }
  s
}

hessienne = function(X, Y, par){
  mat = diag( as.vector( sigma(par%*%t(X))*(1-sigma(par%*%t(X))) ) , nrow = length(Y))
  -t(X) %*% mat %*% X
}

estim = function(X, Y, ini, horizon){
  par = ini
  for (i in 1:horizon){
    print(par)
    par = par - solve(hessienne(X, Y, par)) %*% gradient(X, Y, par)
    par = as.vector(par)
  }
  par
}

# Exemple
par = c(1,-0.5,2)

# Visualisation du classifieur bayesien (parametres connus)
xgrid = c()
y = c()
for (x1 in seq(-4,4,by=0.1)){
  for (x2 in seq(-4,4,by=0.1)){
    x = c(1,x1,x2)
    xgrid = rbind(xgrid , x)
    y = c(y, classifieur(x,par))
  }
}
par(mfrow=c(1,1))
plot(xgrid[,2],xgrid[,3],col=(y+1)/2+1,xlab='X1',ylab='X2')

# Risque bayesien
ntest = 1000
test = simulation(ntest,par)
Xtest = test$X
Ytest = test$Y

err_test = 0
for (i in 1:ntest){
  y_hat = classifieur(Xtest[i,] , par)
  if (y_hat != Ytest[i]){
    err_test = err_test + 1
  }
}
err_test = err_test / ntest
print(paste("Regression logistique : risque bayesien =",err_test))

# Estimation des parametres

ntrain = 25
train = simulation(ntrain,par)
Xtrain = train$X
Ytrain = train$Y

par_hat = estim(Xtrain, Ytrain, c(0,0,0), 200)

# Visualisation de la regle apprise
y2 = c()
for (x1 in seq(-4,4,by=0.1)){
  for (x2 in seq(-4,4,by=0.1)){
    x = c(1,x1,x2)
    y2 = c(y2, classifieur(x,par_hat))
  }
}

par(mfrow=c(1,2))
plot(xgrid[,2],xgrid[,3],col=(y+1)/2+1,xlab='X1',ylab='X2',main='Vrai parametre')
plot(xgrid[,2],xgrid[,3],col=(y2+1)/2+1,xlab='X1',ylab='X2',main='Parametre estime')


# Erreur du classifieur entraine
err_test_hat = 0
for (i in 1:ntest){
  y_hat = classifieur(Xtest[i,] , par_hat)
  if (y_hat != Ytest[i]){
    err_test_hat = err_test_hat + 1
  }
}
err_test_hat = err_test_hat / ntest
print(paste("Regression logistique : risque bayesien (plug-in) =",err_test_hat))
print(paste("Regression logistique : risque bayesien =",err_test))
