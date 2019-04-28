ProfLike <- function(
  theta, 
  zt,
  X,
  dist.hydro,
  weight,
  net.zero,
  a.mat,
  b.mat,
  x.dat,
  y.dat,
  CorModels,
  use.nugget,
  use.anisotropy,
  EstMeth,
  REs, 
  scale,
  maxrang = NULL
){

	LL <- -2*m2LL.stream1(theta, m2LLdata = zt,
		X = X, dist.hydro = dist.hydro, weight =
		weight, net.zero = net.zero,
		a.mat = a.mat, b.mat = b.mat,
		x.dat = x.dat, y.dat = y.dat,
		CorModels = CorModels, 
		use.nugget = use.nugget, 
		use.anisotropy = use.anisotropy,
		EstMeth = EstMeth, 
		REs = REs, 
		scale = scale,
		maxrang = maxrang
	)
						
	return(LL)

}

