getPsi <-
function(mat,mat2=NULL,...){
	#R	
	dots<-list(...);
	doCheck = TRUE;
	if (!is.null(dots$doCheck)){doCheck = dots$doCheck;}
	standardChecks <- .doStandardChecks(mat1=mat,mat2=mat2,doCheck=doCheck);
	mat <-standardChecks$mat1
	mat2<-standardChecks$mat2
	if (length(mat2)==0){
		psi<-.Call("getPsi202",mat)
		return(psi)
	} else{
		psi1<-.Call("getPsi202",mat)
		psi2<-.Call("getPsi202",mat2)
		return(c(psi1,psi2))
	}
}

.f1<-function(b){
	.f2<-function(i){
		if (any(i>50)){stop("Not save to calculate values i>50");}	
		if (length(i)>1){return(sapply(i,.f2))}
		1+sum((-1)^(0:i)*4/((0:i)+1)*choose(i-1,c(0:i)-1))
	}
	if (length(b)>1){return(sapply(b,f1))}
	9+sum(.f2(0:(b-2)))
}
getVar<-function(mat,...){
	#R
	dots<-list(...);
	doCheck = TRUE;
	if (!is.null(dots$doCheck)){doCheck = dots$doCheck;}
	standardChecks <- .doStandardChecks(mat1=mat,mat2=NULL,doCheck=doCheck);
	mat <-standardChecks$mat1

	bn<-rowSums(is.finite(mat))
	ub<-table(bn)
	T<-sum(bn)

	v<-0
	for (i in c(1:length(ub))){
		Bi<-as.numeric(names(ub)[i])
		cii<-(Bi+3)*(T+1)/(45*Bi*(Bi-1)*(T-Bi)) *(Bi/T)^2
		v<-v+cii*ub[i]
		for (j in c(1:length(ub))){
			Bj<-as.numeric(names(ub)[j])
			cjj<-(Bj+3)*(T+1)/(45*Bj*(Bj-1)*(T-Bj))*(Bj/T)^2
			v = v-(ub[j]-(i==j))*ub[i]*sqrt(cii*cjj*Bi*Bj*Bi*Bj/((T-Bi)*(T-Bj)*.f1(Bi)*.f1(Bj)))
		}	
	}
	return(v)
}


 






