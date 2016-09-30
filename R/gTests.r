
## main function
g.tests = function(E, sample1ID, sample2ID, test.type="all", perm=0){
  temp = getR1R2(E, sample1ID)
  R1 = temp$R1
  R2 = temp$R2
  n = length(sample1ID)
  m = length(sample2ID)
  N = n+m
  Ebynode = vector("list", N)
  for (i in 1:N) Ebynode[[i]] = rep(0,0)
  for (i in 1:nrow(E)){
    Ebynode[[E[i,1]]] = c(Ebynode[[E[i,1]]], E[i,2])
    Ebynode[[E[i,2]]] = c(Ebynode[[E[i,2]]], E[i,1])
  }
  nE = nrow(E)
  nodedeg = rep(0,N)
  for (i in 1:N) nodedeg[i] = length(Ebynode[[i]])
  nEi = sum(nodedeg*(nodedeg-1))  # pair of nodes sharing a node * 2
  mu0 = nE*2*n*m/N/(N-1)
  mu1 = nE*n*(n-1)/N/(N-1)
  mu2 = nE*m*(m-1)/N/(N-1)
  V0 = nEi * n*m/N/(N-1) + (nE*(nE-1)-nEi) * 4*n*m*(n-1)*(m-1)/N/(N-1)/(N-2)/(N-3) + mu0 - mu0^2
  V1 = nEi * n*(n-1)*(n-2)/N/(N-1)/(N-2) + (nE*(nE-1)-nEi) * n*(n-1)*(n-2)*(n-3)/N/(N-1)/(N-2)/(N-3) + mu1 - mu1^2
  V2 = nEi * m*(m-1)*(m-2)/N/(N-1)/(N-2) + (nE*(nE-1)-nEi) * m*(m-1)*(m-2)*(m-3)/N/(N-1)/(N-2)/(N-3) + mu2 - mu2^2
  V12 = (nE*(nE-1)-nEi) * m*n*(m-1)*(n-1)/N/(N-1)/(N-2)/(N-3) - mu1*mu2
  S = matrix(c(V1,V12,V12,V2), nrow=2)
  if (is.na(match(test.type,c("all","original","o","generalized","g","weighted","w")))){
    cat("Wrong test.type input! All tests are performed!\n")
    test.type="all"
  }
  if (test.type=="all" || test.type=="original" || test.type=="o"){
    Zo = (nE-R1-R2-mu0)/sqrt(V0)
    po.approx = pnorm(Zo)
    ro = list(test.statistic=Zo, pval.approx=po.approx)
  }
  if (test.type=="all" || test.type=="generalized" || test.type=="g"){
    Sinv = solve(S)
    Rmv = c(R1-mu1, R2-mu2)
    Zg = Rmv %*% Sinv %*% Rmv
    pg.approx = pchisq(Zg, df=2, lower.tail=F)
    rg = list(test.statistic=Zg, pval.approx=pg.approx)
  }
  if (test.type=="all" || test.type=="weighted" || test.type=="w"){
    Zw = ((m)*(R1-mu1)+(n)*(R2-mu2))/sqrt((m)^2*V1+(n)^2*V2+2*(m)*(n)*V12)
    pw.approx = pnorm(-Zw)
    rw = list(test.statistic=Zw, pval.approx=pw.approx)
  }
  if (perm>0){
    Zov = Zgv = Zwv = rep(0,perm)
    for (k in 1:perm){
      g = sample(c(sample1ID, sample2ID), n)
      temp.p = getR1R2(E,g)
      R1.p = temp.p$R1
      R2.p = temp.p$R2
      if (test.type=="all" || test.type=="original" || test.type=="o"){
        Zov[k] = (nE-R1.p-R2.p-mu0)/sqrt(V0)
      }
      if (test.type=="all" || test.type=="generalized" || test.type=="g"){
        Rmv.p = c(R1.p-mu1, R2.p-mu2)
        Zgv[k] = Rmv.p %*% Sinv %*% Rmv.p
      }
      if (test.type=="all" || test.type=="weighted" || test.type=="w"){
        Zwv[k] = ((m)*(R1.p-mu1)+(n)*(R2.p-mu2))/sqrt((m)^2*V1+(n)^2*V2+2*(m)*(n)*V12)
      }
    }
    if (test.type=="all" || test.type=="original" || test.type=="o"){
      po.perm = length(which(Zov<=Zo))/perm
      ro = c(ro, list(pval.perm=po.perm))
    }
    if (test.type=="all" || test.type=="generalized" || test.type=="g"){
      pg.perm = length(which(Zgv>=Zg[1]))/perm
      rg = c(rg, list(pval.perm=pg.perm))
    }
    if (test.type=="all" || test.type=="weighted" || test.type=="w"){
      pw.perm = length(which(Zwv>=Zw))/perm
      rw = c(rw, list(pval.perm=pw.perm))
    }
  }
  r = list()
  if (test.type=="all" || test.type=="original" || test.type=="o"){
    r = c(r,list(original=ro))
  }
  if (test.type=="all" || test.type=="generalized" || test.type=="g"){
    r = c(r,list(generalized=rg))
  }
  if (test.type=="all" || test.type=="weighted" || test.type=="w"){
    r = c(r,list(weighted=rw))
  }
  return(r)
}


## supporting function
getR1R2 = function(E, G1){
  R1 = R2 = 0
  for (i in 1:nrow(E)){
    e1 = is.na(match(E[i,1],G1))
    e2 = is.na(match(E[i,2],G1))
    if ((!e1) && (!e2))  R1 = R1 + 1
    if (e1 && e2)  R2 = R2 + 1
  }
  return(list(R1=R1, R2=R2))
}
