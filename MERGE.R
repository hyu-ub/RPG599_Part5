# Parameter 'x_norm' is standardized matrix of the gene expression values (dimensions #genes x #samples)
# Parameter 'y_norm' is standardized matrix of the drug sensitivity values (dimensions #drugs x #samples)
# Parameter 'Fnorm' is standardized matrix of the driver feature values (dimensions #genes x #features)
# Parameter 'lambda0' is the tuning parameter for regularization
# Parameter 'init' defines how to initialize the driver features weight vector v (of length #features); options for the parameter are 'zero' and 'rnd'
  
  MERGE = function (x_norm, y_norm, Fnorm, lambda0, init='zero') {
	ngenes = nrow(x_norm)
	nsamples = ncol(x_norm)
	ndrugs = nrow(y_norm)

	xy = x_norm %*% t(y_norm)
	if (init == 'zero') {
		vinit = matrix(0, ncol=1, nrow=ncol(Fnorm))
	} else if (init == 'rnd') {
		vinit = matrix(rnorm(ncol(Fnorm)), ncol=1)
	}
	v = vinit
	
	lambda = lambda0 - Fnorm %*% v

	for (iter in 1:10) {
		W = xy / ( nsamples - 1 + lambda[,rep(1,ndrugs)])
		sumW2 = matrix(rowSums(W^2), nrow=1)

		obj0 = 0
		for (i in 1:ngenes) {
			obj0 = obj0 + sum((t(W[i,,drop=F]) %*% x_norm[i,,drop=F] - y_norm)^2)
		}	

		gradinit = -t(Fnorm) %*% t(sumW2)
		for (it in 1:10) {
			tmp = rowSums(Fnorm * t(as.matrix(v)[,rep(1,ngenes)]))
			grad = gradinit + ndrugs * colSums(Fnorm / (as.matrix(lambda0 - tmp)[,rep(1,length(v))]))
			grad2 = sum(grad^2)

			tt = 0.002
			alphaa = 0.3
			betaa = 0.8

			preobj = obj0 + sumW2 %*% (lambda0 - Fnorm %*% v) - ndrugs * sum(log(lambda0 - tmp))
			for (biter in 1:200) {
				newv = v - tt * grad
				if(sum(Fnorm %*% newv <= lambda0) == ngenes) {
					tmp = rowSums(Fnorm * t(as.matrix(newv)[,rep(1,ngenes)]))
					postobj = obj0 + sumW2 %*% (lambda0 - Fnorm %*% newv) - ndrugs * sum(log(lambda0 - tmp))
					if(is.na(preobj) || postobj < preobj - alphaa * tt * grad2) { # duzgun bir v buldugun anda cik
						break
					} else {
						tt = betaa * tt
					}
				} else {
					tt = betaa * tt
				}
			}
			v = newv
		}
		lambda = lambda0 - Fnorm %*% v

		print(paste0('iter = ', iter, ' min = ', min(lambda), '  max = ', max(lambda)))
	}
	return(list('v'=v, 'W'=W, 'lambda'=lambda, 'sumW2'=sumW2, 'vinit'=vinit))
}
