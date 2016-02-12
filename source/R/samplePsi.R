samplePsi<-function(b,n=NULL,r=0){
	bn<-b;
	if  (is.null(n)){
		if ( length(b)<2 ) { stop("If 'n' is not given 'b' should be of length >1");}
		bn<-b;
		n<-length(bn);
	} else {
		if ( length(b)!=1 ) { stop("If 'n' is given 'b' should be of length 1");}
		bn<-rep(b,n);
	}
	chol<-matrix(r,ncol=max(bn),nrow=max(bn));diag(chol)<-1;chol<-chol(chol)
	psi<-.Call("samplePsi",chol,as.integer(bn),NULL,as.integer(options("concordance.nDraws")$concordance.nDraws),as.integer(options("concordance.nCPU")$concordance.nCPU),as.integer(options("concordance.seed")$concordance.seed));
	psi
}


sampleVariance<-function(mat1,mat2=NULL,...){
	dots<-list(...);
	doCheck = TRUE;
	if (!is.null(dots$doCheck)){doCheck = dots$doCheck;}
	standardChecks <- .doStandardChecks(mat1=mat1,mat2=mat2,doCheck=doCheck);
	mat1 <-standardChecks$mat1
	mat2<-standardChecks$mat2

	mat<-mat1;
	if (!is.null(mat2)){ mat<-cbind(mat,mat2);}
	chol<-chol(cor(mat,method='spearman',use='pairw'))
	
	missingmat1<-is.finite(mat1)
	missingmat2<-is.finite(mat2)
	psi<-.Call("samplePsi",chol,missingmat1,missingmat2,as.integer(options("concordance.nDraws")$concordance.nDraws),as.integer(options("concordance.nCPU")$concordance.nCPU),as.integer(options("concordance.seed")$concordance.seed));
	psi;
}

doBootstrap<-function(mat1,mat2,...){
	dots<-list(...);
	doCheck = TRUE;

	if (!is.null(dots$doCheck)){doCheck = dots$doCheck;}
	standardChecks <- .doStandardChecks(mat1=mat1,mat2=mat2,doCheck=doCheck);
	mat1 <-standardChecks$mat1
	mat2 <-standardChecks$mat2


	psi<-.Call("doBootstrap203",mat1,mat2,as.integer(options("concordance.nDraws")$concordance.nDraws),as.integer(options("concordance.nCPU")$concordance.nCPU),as.integer(options("concordance.seed")$concordance.seed));
	psi;
}



