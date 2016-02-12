exactProb<-function(bn, version=c("2.02")){
	version=match.arg(version)
	storage.mode(bn) <- "integer"
	if(version=="2.02"){	return(.Call('exactDistr202',bn,as.logical(options("concordance.verbose")$concordance.verbose)))}
}
 

