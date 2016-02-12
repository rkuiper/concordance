concordance.test <-
function(mat,mat2=NULL,alternative=NULL,exact=FALSE){
	#R	
	if (is.null(alternative) & is.null(mat2)){
		alternative=match.arg(alternative,c('greater','less','two.sided')) #Null hypothesis psi<=2/3
	} else {	
		alternative=match.arg(alternative,c('two.sided','greater','less')) #Null hypothesis |psi1-psi2|==0
	}
	cl <- match.call()
	standardChecks <- .doStandardChecks(mat=mat,mat2=mat2,doCheck=TRUE);
	mat <- standardChecks$mat1
	mat2 <- standardChecks$mat2
	
	
	if (exact & !is.null(mat2)){ warning('No exact distribution is calculated for difference between concordances'); }

	obs<-getPsi(mat,mat2,doCheck=FALSE)

	p<-NA
	if (exact){	
		exactDist<-list();
		bn<- rowSums(is.finite(mat))
		sumN<-length(bn)
		exactDist<-exactProb(bn); 
		
		if (alternative=='less'){
			exact_pvalue<-sum(exactDist[((exactDist[,1])<=obs),2])
		} else if (alternative=='two.sided'){
			exact_pvalue<- sum(exactDist[ ! ((exactDist[,1] > 2/3-abs(obs-2/3)) &  (exactDist[,1] < 2/3+abs(obs-2/3))) ,2])

		} else if (alternative=='greater'){
			exact_pvalue<-sum(exactDist[ ((exactDist[,1])>=obs),2])
		}
		p<-exact_pvalue
	} else {
		if (!is.null(mat2)){
			obsvar<-var(doBootstrap(mat1=mat,mat2=mat2,doCheck=FALSE));
			delta<-diff(obs)
			p.norm<-NA
			if (alternative=='less') p.norm<-pnorm(delta,mean=0,sd=sqrt(obsvar) ,lower.tail=T)
			else if (alternative=='greater') p.norm<-pnorm(delta,mean=0,sd=sqrt(obsvar) ,lower.tail=F)  
			else if (alternative=='two.sided')  p.norm<-2*pnorm(abs(delta),mean=0,sd=sqrt(obsvar) ,lower.tail=F)
			p<-p.norm;

		}
		if (is.null(mat2)){
			obsvar<-getVar(mat,doCheck=FALSE)
			bn<-rowSums(is.finite(mat))
			meanB<-mean(bn)
			x<-(6*exp(1)*meanB-10*exp(1))/(9*exp(1)*meanB-15*exp(1)+9)
			alpha<--1 * ((2/3-x)*obsvar + (2/3-x)^3 - (2/3-x)^2)/obsvar
			beta<--1 * (alpha * (3*x+1))/(3*x-2)

	
			if (alternative=='less') p.beta<-pbeta(obs-x,shape1=alpha,shape2=beta,lower.tail=T)
			else if (alternative=='greater') p.beta<-pbeta(obs-x,shape1=alpha,shape2=beta,lower.tail=F)
			else if (alternative=='two.sided') {
				p.beta.left<-pbeta(alpha/(alpha+beta)-(abs((obs-x)-alpha/(alpha+beta))) ,shape1=alpha,shape2=beta,lower.tail=T)
				p.beta.right<-pbeta(alpha/(alpha+beta)+(abs((obs-x)-alpha/(alpha+beta))) ,shape1=alpha,shape2=beta,lower.tail=F)
				p.beta<-p.beta.left+ p.beta.right;
			}
			p<-p.beta	
		}	
	}
	
	result<-list(p.value = p,statistic = obs,exact=exact,alternative=alternative,call=cl)
	class(result)<-"concordance.test"
	result
}

print.concordance.test<-function(x){
	#R
	if (length(x$statistic)==1){
		if (x$alternative=="greater"){
			cat("One sided concordance test\nAlternative hypothesis: psi > 2/3\n")
		}
		else if (x$alternative=="less"){
			cat("One sided concordance test\nAlternative hypothesis: psi < 2/3\n")
		}
		else if (x$alternative=="two.sided"){
			cat("Two sided concordance test\nAlternative hypothesis: |psi| > 2/3\n")
		}
		cat("Statistic = " , signif(x$statistic,3),"\n")
		cat("Pvalue = " , signif(x$p.value,3),"\n")
	} else if (length(x$statistic)==2){
		if (x$alternative=="greater"){
			cat("One sided concordance test for difference between two concordances\nAlternative hypothesis: psi1 > psi2\n")
		}
		else if (x$alternative=="less"){
			cat("One sided concordance test for difference between two concordances\nAlternative hypothesis: psi1 < 2/3\n")
		}
		else if (x$alternative=="two.sided"){
			cat("Two sided concordance test for difference between two concordances\nAlternative hypothesis: |psi1 - psi2| > 0\n")
		}
		cat("Statistics = " , signif(x$statistic[1],3),"(psi1); ", signif(x$statistic[2],3),"(psi2); \n")
		cat("Pvalue = " , signif(x$p.value,3),"\n")
	} else { stop("Not implemented situation");
	
	}
}
