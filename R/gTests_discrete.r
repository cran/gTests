## main function
g.tests_discrete = function(E, counts, test.type = "all", maxtype.kappa = 1.14, perm = 0){ 	
	v = counts[,1] + counts[,2]
	ids = which(v!=0)
	v1 = counts[ids,1]
	v2 = counts[ids,2]
	vmat = counts[ids,]
	colnames(vmat) = NULL
	
	temp = getMV_discrete(E,vmat)
	mu1_a = temp$mu1_a
	mu2_a = temp$mu2_a
	var1_a = temp$var1_a
	var2_a = temp$var2_a
	var12_a = temp$var12_a
	mu1_u = temp$mu1_u
	mu2_u = temp$mu2_u
	var1_u = temp$var1_u
	var2_u = temp$var2_u
	var12_u = temp$var12_u

	p0 = q0 = 0.5
	qw = (sum(v2)-1)/(sum(v1+v2)-2)
	pw = 1-qw
		
	# mu and variance of the original edge-count test (average method)	
	muo_a = q0*mu1_a + p0*mu2_a 
	varo_a = q0^2*var1_a + p0^2*var2_a + 2*q0*p0*var12_a
	# mu and variance of the weighted edge-count test (average method)
	muw_a = qw*mu1_a + pw*mu2_a 
	varw_a = qw^2*var1_a + pw^2*var2_a + 2*qw*pw*var12_a
	# mu and variance of the difference of two with-in group edge-counts (average method)
	mud_a = mu1_a - mu2_a 
	vard_a = var1_a + var2_a - 2*var12_a
	
	# mu and variance of the original edge-count test (union method)
	muo_u = q0*mu1_u + p0*mu2_u
	varo_u = q0^2*var1_u + p0^2*var2_u + 2*q0*p0*var12_u
	# mu and variance of the weighted edge-count test (union method)
	muw_u = qw*mu1_u + pw*mu2_u 
	varw_u = qw^2*var1_u + pw^2*var2_u + 2*qw*pw*var12_u
	# mu and variance of the difference of two with-in group edge-counts (union method)
	mud_u = mu1_u - mu2_u 
	vard_u = var1_u + var2_u - 2*var12_u
	
	temp = getR1R2_discrete(E,vmat)
	R1_a = temp$R1_a
	R2_a = temp$R2_a
	R1_u = temp$R1_u
	R2_u = temp$R2_u
	
	Ro_a = q0*R1_a + p0*R2_a
	Rw_a = qw*R1_a + pw*R2_a
	Rd_a = R1_a - R2_a	
	Zw_a = (Rw_a-muw_a)/sqrt(varw_a)
	Zd_a = (Rd_a-mud_a)/sqrt(vard_a)	
	
	Ro_u = q0*R1_u + p0*R2_u
	Rw_u = qw*R1_u + pw*R2_u
	Rd_u = R1_u - R2_u	
	Zw_u = (Rw_u-muw_u)/sqrt(varw_u)
	Zd_u = (Rd_u-mud_u)/sqrt(vard_u)	

	if (is.na(match(test.type,c("all","original","o","generalized","g","weighted","w","maxtype","m")))){
		cat("Wrong test.type input! All tests are performed!\n")
		test.type="all"
	}
	if (test.type=="all" || test.type=="original" || test.type=="o"){
		Zo_a = -(Ro_a-muo_a)/sqrt(varo_a)
		Zo_u = -(Ro_u-muo_u)/sqrt(varo_u)
		po.approx_a = pnorm(Zo_a)
		po.approx_u = pnorm(Zo_u)
		ro = list(test.statistic_a=Zo_a, pval.approx_a=po.approx_a, test.statistic_u=Zo_u, pval.approx_u=po.approx_u)
	}
	if (test.type=="all" || test.type=="generalized" || test.type=="g"){
		S_a = Zw_a^2 + Zd_a^2
		S_u = Zw_u^2 + Zd_u^2
		pg.approx_a = pchisq(S_a, df=2, lower.tail=F)
		pg.approx_u = pchisq(S_u, df=2, lower.tail=F)
		rg = list(test.statistic_a=S_a, pval.approx_a=pg.approx_a, test.statistic_u=S_u, pval.approx_u=pg.approx_u)
	}
	if (test.type=="all" || test.type=="weighted" || test.type=="w"){
		pw.approx_a = pnorm(-Zw_a)
		pw.approx_u = pnorm(-Zw_u)
		rw = list(test.statistic_a=Zw_a, pval.approx_a=pw.approx_a, test.statistic_u=Zw_u, pval.approx_u=pw.approx_u)
	}
	if (test.type=="all" || test.type=="maxtype" || test.type=="m"){
		M_a = max(maxtype.kappa*Zw_a,abs(Zd_a))
		M_u = max(maxtype.kappa*Zw_u,abs(Zd_u))
		pm.approx_a = 1-pnorm(M_a/maxtype.kappa)*(2*pnorm(M_a)-1)
		pm.approx_u = 1-pnorm(M_u/maxtype.kappa)*(2*pnorm(M_u)-1)
		rmax = list(test.statistic_a=M_a, pval.approx_a=pm.approx_a, test.statistic_u=M_u, pval.approx_u=pm.approx_u)
	}

	if (perm > 0){
		Zov_a = Sv_a = Zwv_a = Mv_a = rep(0,perm)
		Zov_u = Sv_u = Zwv_u = Mv_u = rep(0,perm)
		for (k in 1:perm){
			vmatnew = permute_discrete(vmat)
			colnames(vmatnew) = NULL
			
			temp.p = getR1R2_discrete(E,vmatnew)
			R1.p_a = temp.p$R1_a
			R2.p_a = temp.p$R2_a
			R1.p_u = temp.p$R1_u
			R2.p_u = temp.p$R2_u
			
			Ro.p_a = q0*R1.p_a + p0*R2.p_a
			Rw.p_a = qw*R1.p_a + pw*R2.p_a
			Rd.p_a = R1.p_a - R2.p_a	
			Zwv_a[k] = (Rw.p_a-muw_a)/sqrt(varw_a)
			Zd.p_a = (Rd.p_a-mud_a)/sqrt(vard_a)	
	
			Ro.p_u = q0*R1.p_u + p0*R2.p_u
			Rw.p_u = qw*R1.p_u + pw*R2.p_u
			Rd.p_u = R1.p_u - R2.p_u	
			Zwv_u[k] = (Rw.p_u-muw_u)/sqrt(varw_u)
			Zd.p_u = (Rd.p_u-mud_u)/sqrt(vard_u)	
	
			if (test.type=="all" || test.type=="original" || test.type=="o"){
				Zov_a[k] = -(Ro.p_a-muo_a)/sqrt(varo_a)
				Zov_u[k] = -(Ro.p_u-muo_u)/sqrt(varo_u)
			}
			if (test.type=="all" || test.type=="generalized" || test.type=="g"){
				Sv_a[k] = Zwv_a[k]^2 + Zd.p_a^2
				Sv_u[k] = Zwv_u[k]^2 + Zd.p_u^2
			}
			if (test.type=="all" || test.type=="maxtype" || test.type=="m"){
				Mv_a[k] = max(maxtype.kappa*Zwv_a[k],abs(Zd.p_a))
				Mv_u[k] = max(maxtype.kappa*Zwv_u[k],abs(Zd.p_u))
			}
		}
		if (test.type=="all" || test.type=="original" || test.type=="o"){
			po.perm_a = (length(which(Zov_a<=Zo_a)) + length(which(Zov_a<Zo_a)))/2/perm
			po.perm_u = (length(which(Zov_u<=Zo_u)) + length(which(Zov_u<Zo_u)))/2/perm
			ro = c(ro, list(pval.perm_a=po.perm_a, pval.perm_u=po.perm_u))
		}
		if (test.type=="all" || test.type=="generalized" || test.type=="g"){
			pg.perm_a = (length(which(Sv_a>=S_a)) + length(which(Sv_a>S_a)))/2/perm
			pg.perm_u = (length(which(Sv_u>=S_u)) + length(which(Sv_u>S_u)))/2/perm
			rg = c(rg, list(pval.perm_a=pg.perm_a, pval.perm_u=pg.perm_u))
		}
		if (test.type=="all" || test.type=="weighted" || test.type=="w"){
			pw.perm_a = (length(which(Zwv_a>=Zw_a)) + length(which(Zwv_a>Zw_a)))/2/perm
			pw.perm_u = (length(which(Zwv_u>=Zw_u)) + length(which(Zwv_u>Zw_u)))/2/perm
			rw = c(rw, list(pval.perm_a=pw.perm_a, pval.perm_u=pw.perm_u))
		}
		if (test.type=="all" || test.type=="maxtype" || test.type=="m"){
			pm.perm_a = (length(which(Mv_a>=M_a)) + length(which(Mv_a>M_a)))/2/perm
			pm.perm_u = (length(which(Mv_u>=M_u)) + length(which(Mv_u>M_u)))/2/perm
			rmax = c(rmax, list(pval.perm_a=pm.perm_a, pval.perm_u=pm.perm_u))
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
	if (test.type=="all" || test.type=="maxtype" || test.type=="m"){
		r = c(r,list(maxtype=rmax))
	}
	return(r)
}

## generate a permutation 
permute_discrete = function(vmat){
	m = vmat[,1] + vmat[,2]
	N = sum(m)
	X = rep(0, sum(m))
	id = sample(1:N, sum(vmat[,1]))
	X[id] = 1
	K = dim(vmat)[1]
	newvec1 = rep(0,K)
	newvec2 = rep(0,K)
	for (i in 1:K){
		j1 = sum(m[0:(i-1)])+1
		j2 = sum(m[0:i])
		newvec1[i] = sum(X[j1:j2])
		newvec2[i] = sum(1-X[j1:j2])
	}
	return(cbind(newvec1, newvec2))
}

## supporting function
getR1R2_discrete = function(E,vmat){
	SA_a = SB_a = SA_u = SB_u = rep(0,2)
	catenum = vmat[,1] + vmat[,2]
	len = dim(vmat)[1]
	SA_a = colSums(vmat*(vmat-1)/matrix(rep(catenum,2),nrow=len)) # 1*2
	SA_u = colSums(vmat*(vmat-1)/2) # 1*2
	prob_vmat = vmat/matrix(rep(catenum,2),nrow=len)
	
	for (e in 1:nrow(E)){
		SB_a = SB_a + prob_vmat[E[e,1],]*prob_vmat[E[e,2],]
		SB_u = SB_u + vmat[E[e,1],]*vmat[E[e,2],]
	}

	R1_a = SA_a[1] + SB_a[1]
	R2_a = SA_a[2] + SB_a[2]
	R1_u = SA_u[1] + SB_u[1]
	R2_u = SA_u[2] + SB_u[2]
	
	return(list(R1_a = R1_a, R2_a = R2_a, R1_u = R1_u, R2_u = R2_u))
}

getMV_discrete = function(E,vmat){
	K = dim(vmat)[1]
	na = sum(vmat[,1])
	nb = sum(vmat[,2])
	m = vmat[,1] + vmat[,2]
	N = na + nb
	p1 = na*(na-1)/N/(N-1)
	p2 = p1*(na-2)/(N-2)
	p3 = p2*(na-3)/(N-3)
	q1 = nb*(nb-1)/N/(N-1)
	q2 = q1*(nb-2)/(N-2)
	q3 = q2*(nb-3)/(N-3)
	f1 = p1*nb*(nb-1)/(N-2)/(N-3)
	nodedeg_a = nodedeg_u = rep(0,K)
	sumE = dim(E)[1]
	quan = quan_u = 0
	mu1_a = (N-K+sumE)*p1
	mu2_a = (N-K+sumE)*q1

	for (i in 1:sumE){
		e1 = E[i,1]
		e2 = E[i,2]
		nodedeg_a[e1] = nodedeg_a[e1]+1
		nodedeg_a[e2] = nodedeg_a[e2]+1
		nodedeg_u[e1] = nodedeg_u[e1] + m[e2]
		nodedeg_u[e2] = nodedeg_u[e2] + m[e1]
		quan= quan + 1/m[e1]/m[e2]
		quan_u = quan_u + m[e1]*m[e2]
	}
	temp1 = N-K+2*sumE+sum(nodedeg_a^2/4/m)-sum(nodedeg_a/m)
	temp2 = K-sum(1/m)
	temp3 = quan
	temp4 = (N-K+sumE)^2
	Gedge = sum(m*(m-1))/2 + quan_u
	nodedeg_u = nodedeg_u + m - 1

	var1_a = 4*(p2-p3)*temp1+2*(p1-4*p2+3*p3)*temp2+(p1-2*p2+p3)*temp3+(p3-p1^2)*temp4
	var2_a = 4*(q2-q3)*temp1+2*(q1-4*q2+3*q3)*temp2+(q1-2*q2+q3)*temp3+(q3-q1^2)*temp4
	var12_a = f1*(-4*temp1+6*temp2+temp3)+(f1-p1*q1)*temp4

	mu1_u = Gedge*p1
	mu2_u = Gedge*q1
	var1_u = (p1-p3)*Gedge + (p2-p3)*sum(m*nodedeg_u*(nodedeg_u-1)) + (p3-p1^2)*Gedge^2
	var2_u = (q1-q3)*Gedge + (q2-q3)*sum(m*nodedeg_u*(nodedeg_u-1)) + (q3-q1^2)*Gedge^2
	var12_u = f1*(Gedge^2-Gedge-sum(m*nodedeg_u*(nodedeg_u-1))) - p1*q1*Gedge^2
	return(list(mu1_a=mu1_a,mu2_a=mu2_a,var1_a=var1_a,var2_a=var2_a,var12_a=var12_a,
			mu1_u=mu1_u,mu2_u=mu2_u,var1_u=var1_u,var2_u=var2_u,var12_u=var12_u))
}

