doPerm <-
function(n,b,nperm=10000,seed=1){
if (missing(n)){n<-NULL;}
ncpu = options("nCPU.concordance")$nCPU.concordance;
if (nperm<(1e6+1)){ncpu=1;}

if (!require("snow")){stop('This function requires packages "snow" and "rlecuyer"')}
if (ncpu>1){
	cl <- makeSOCKcluster(rep("localhost",ncpu),prngkind="default")
	if (!is.null(seed)){
	clusterSetupRNG(cl, seed=rep(seed,6))
}
clusterExport(cl, "doPerm")
res<-unlist(clusterCall(cl, doPerm,b=b,n=n,nperm=ceiling(nperm/ncpu),seed=seed,ncpu=1))
stopCluster(cl)
return(res)
}
orgSeed<- get(".Random.seed", .GlobalEnv)
set.seed(seed)

if (is.null(seed)){seed<- -1;}
if (length(b)==1){b<-rep(b,n);}
tab<-table(b);
ub<-as.integer(names(tab));
nb<-as.integer(tab);
nub<-as.integer(length(ub));
nperm<-as.integer(nperm);
psi<-as.double(rep(0,nperm))
seed<-as.integer(seed)
assign(".Random.seed", orgSeed, .GlobalEnv)

.C('doPerm',ub,nb,nub,nperm,seed,psi)[[6]]
}
