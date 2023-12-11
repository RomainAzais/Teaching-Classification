# Chargement et visualisation des donnees MNIST

library(dslabs)
dat = read_mnist()

# Donnees d'apprentissage
x_train = dat$train$images
y_train = dat$train$labels
n_train = length(y_train)

# Donnees de test
x_test = dat$test$images
y_test = dat$test$labels
n_test = length(y_test)

# Exemple
i=5
x = x_train[i,]
y = y_train[i]

visu_mnist = function(x,y=NULL){
  im = matrix(x,nrow=28)[,28:1]
  if (is.null(y)){
    lab_ = 'unknown'
  } else {
    lab_ = y
  }
  image(1:28, 1:28, im, col = gray(seq(0, 1, 0.05)), xlab = "", ylab="",main=paste("label=",lab_,sep=''))
}

par(mfrow=c(1,1))
visu_mnist(x,y)
