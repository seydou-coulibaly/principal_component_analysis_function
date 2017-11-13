# pré-traitement des données
x = read.table(file.choose())
View(x)
y = x[,c(3,6,7,8,9,10,11,12,13)]
y = as.matrix(y)
colnames(y) = y[1,]
View(y)
y = y[-1,]
y = apply(y,c(1,2),as.numeric)
View(y)

# Description univariée
summary(y)
round (apply(y,2,var),3)
round (apply(y,2,sd),3)
# Paramètre de formes
install.packages("e1071")
library(e1071)
for (i in 1:9 )
{
  print(round(skewness(y[, i],na.rm=T),3))
  print(round(kurtosis(y[,i],na.rm=T),3))
}
# Histogramme
layout(matrix(c(1:9),3,3))
for(i in 1:9) { hist(y[,i],main=colnames(y)[i],xlab="")}
layout(1)
# Densiste
layout(matrix(c(1:9),3,3))
for(i in 1:9) {plot(density(y[,i], na.rm = TRUE), main=colnames(y)[i])}
layout(1)
# Description bivariée
summary(x[,15])
cor(y)
pairs(y)
#pie(summary(as.factor(z[,2])), label=z[,2], main=titre)








# Réalisation de la fonction ACP
acpFunction = function(X){
  #nombre de colonnes
  m = ncol(X)
  #nobmbre de lignes
  n = nrow(X)
  #la matrice centrée reduite
  Y = scale(X, center=TRUE, scale=TRUE) * sqrt((n)/(n-1))
  #Y
  #matrice de correlation ou de variance-covariance vu que la matrice est centrée reduite
  R = cor(Y)
  #R
  #valeurs propores et vecteur propres
  eig = eigen(R)
  Lambda = eig$values
  #Lambda
  #vecteurs propores
  U = eig$vectors
  #U
  #inertie
  inertie = sum(diag(R))
  #inertie
  #afficher pourcentage inertie
  pourcentage = (Lambda / inertie) * 100
  #on retiendra les axes qui repondent aux critères de Kaiser
  #c'est à dire ceux dont la valeur propre est supérieur ou égale à 1
  Axes = which(Lambda >= 1)
  #Axes
  #fonction anonyme
  somme = 0
  for (i in 1:length(Axes)) {somme = somme + pourcentage[i]}
  somme
  #print
  #composantes principales
  F = Y %*% U
  #F
  #Qualité de répresentations des individus sur les axes
  Q = (F^2) / (apply((F^2),1,sum))
  #Q
  #Q * 100
  #contribution des individus à la construction des axes
  Ctr = (t(F)^2)/(n * Lambda)
  Ctr = t(Ctr)
  #Ctr
  #Ctr * 100
  #coordonées des variables sur les différents axes
  G = U %*% diag(Lambda^(0.5))
  #G
  #contribution des variables aux axes
  #t(cor(F,Y)^2 / Lambda)
  Ctrv = t(t(G^2)/Lambda)
  #Ctrv
  #Qualité de répresentation des variables suivant les axes
  Qv = G^2
  #Qv
  #Qv * 100
  
  
  #retourner une liste
  liste = list(Lambda,U,length(Axes),somme,F,Q,Ctr,G,Qv,Ctrv)
  #valeurs propres, vecteurs propres,nbre axe retenu,pourcentage inertie cumulé des axes retenu,coordonnées F
  #Qualité Q, contribution Ctr, coordonnées G des varibales, Qualité variables et contribution variables
  return(liste)
}




# Execution de la fonction 
acp = acpFunction(y)
library(ade4)

# Répresentation cercle de corrélation
# G = acp[8]
s.corcircle(acp[8],xax=2,yax=3)


# Répresentations des individus
#avec F
F.range = apply(F, 2, range)
F.range
xlim = c(F.range[1,2], F.range[2,2])
ylim = c(F.range[1,3], F.range[2,3])
plot(F[,2], F[,3], xlim=xlim, ylim=ylim, asp=1, xlab="Axe 2", ylab="Axe 3")
abline(h=0,v=0)
text(F[,2],F[,3],rownames(F))

#arrows(x0=0,y0=0,U[,1],U[,2])

ls()

# L'ACP de la librairie Ade4
z<- dudi.pca(y, center = T, scale = T)
z$eig
#3.0120696 1.8173593 1.1386658 0.8891271 0.7946599 0.5870253 0.3389612 0.2384442 0.1836876
z$eig/sum(z$eig)*100
inertierow<-inertia.dudi(z,row.inertia=TRUE)
inertiecol<-inertia.dudi(z,col.inertia=TRUE)
names(inertierow)

#individus
# contribution =  inertierow$row.abs
# qualite = inertierow$row.rel

# variables
# contribution =  inertiecol$col.abs
# qualite = inertiecol$col.rel

#cercle correlation variables
names(z)
s.corcircle(z$co,xax=1,yax=2)
s.corcircle(z$co,xax=1,yax=3)
s.corcircle(z$co,xax=2,yax=3)


s.label(z$li,xax=2,yax=3)
#s.class(dfxy=z$li,fac=race,col=col,xax=1,yax=2)

#Les variables supplemenetaires
View(x)
x[,15]
x[,14]
x[,2]
supcol(z,x[-1,15])
  #Erreur dans t(as.matrix(Xsup)) %*% (as.matrix(x$l1) * x$lw) : 
  #nécessite des arguments numériques/complexes matrice/vecteur
supcol(z,x[-1,14])
supcol(z,x[-1,2])
