#réalisation de la fonction
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