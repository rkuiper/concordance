.onLoad<-function(libname, pkgname){
	library.dynam("concordance",pkgname,libname)
}
.onunLoad<-function(libname, pkgname){
	library.dynam.unload("concordance",pkgname,libname)
}

setGeneric("summary", function (object, ...) {base::summary(object...);});

.doStandardChecks<-function(mat1,mat2=NULL,doCheck=TRUE){
	#R

	if (doCheck == TRUE) {
		
		 
		if(is.data.frame(mat1)){
			if (!all(unlist(lapply(mat1,class))=="numeric")) stop("mat1 contains non numeric values")
		}	
		if(is.data.frame(mat2)){
			if (!all(unlist(lapply(mat2,class))=="numeric")) stop("mat2 contains non numeric values")
		}
		mat1<-as.matrix(mat1)

		if (length(mat2)==0){mat2<-NULL;}
		else {mat2<-as.matrix(mat2);}

		if(is.null(rownames(mat1))){rownames(mat1)<-c(1:nrow(mat1))}
		if(is.null(colnames(mat1))){colnames(mat1)<-c(1:ncol(mat1))}

		if (!is.null(mat2)){
			if(is.null(rownames(mat2))){rownames(mat2)<-c(1:nrow(mat2))}
			if(is.null(colnames(mat2))){colnames(mat2)<-c(1:ncol(mat2))}
		}

		mat1[!is.finite(mat1)]<-NA

		##Check for ties
		t1<-table(mat1)
		t1.s<-sum(t1[which(t1>1)])
		if (t1.s>0) {warning("Found ",t1.s," ties in mat1")}


		#Check mat2 compatability with mat1 
		if (!is.null(mat2) ){
			mat2[!is.finite(mat2)]<-NA

			##Check for ties
			t2<-table(mat2)
			t2.s<-sum(t2[which(t2>1)])
			if (t2.s>0) {warning("Found ",t2.s," ties in mat2")}

			if(all(dim(mat1)==dim(mat2))){
				rnames1<-rownames(mat1)
				rnames2<-rownames(mat2)
				if (any(rnames1!=rnames2) ) {stop("Rownames of matrices are not in the same order!")}
			}
			else { stop("Matrix dimensions are not comparable!")}

			#Check similarity columnnames matrices
			cnames1<-colnames(mat1)
			cnames2<-colnames(mat2)
			l<-length(intersect(cnames1,cnames2))
			if (l!=length(cnames1) | l!=length(cnames2) ){warning("Colnames of matrices are not comparable!")}
			else if (l==ncol(mat1) & l==ncol(mat2)){
				mat2<-mat2[,cnames1,drop=FALSE];
			}
			if (any(is.finite(mat1)!=is.finite(mat2))){
				warning("Excluding some measurements to make matrices comparable!")
				mat2[!is.finite(mat1)]<-NA
				mat1[!is.finite(mat2)]<-NA
			}
		}
	

		##Exclude observations with less than two replicates
		isOK1<-rowSums(is.finite(mat1))>1
		if (!all(isOK1)){
			mat1<-mat1[isOK1,,drop=FALSE]
			if (!is.null(mat2)){
				mat2<-mat2[isOK1,,drop=FALSE]
			}
		}
		isOK2<-TRUE;
		if (!is.null(mat2)){
			isOK2<-rowSums(is.finite(mat2))>1
			if (!all(isOK2)){
				mat2<-mat2[isOK2,,drop=FALSE]
				if (!is.null(mat1)){
					mat1<-mat1[isOK2,,drop=FALSE]
				}
			}
		}

		if ( (sum(!isOK1)+sum(!isOK2))>0){
			warning(paste("Excluding observations with less than two known replicates (n=",sum(!isOK1)+sum(!isOK2),")",sep=""))
		}	
	}
	mat1[!is.finite(mat1)]<-Inf
	if (!is.null(mat2)) mat2[!is.finite(mat2)]<-Inf
	return(list(mat1=mat1,mat2=mat2))

}

