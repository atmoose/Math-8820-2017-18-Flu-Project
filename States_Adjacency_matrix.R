library(maptools)
library(spdep)
library(maps)


usa.state = map(database="state", fill=TRUE, plot=FALSE)
state.ID <- sapply(strsplit(usa.state$names, ":"), function(x) x[1])
usa.poly = map2SpatialPolygons(usa.state, IDs=state.ID)
usa.nb = poly2nb(usa.poly)
usa.adj.mat = nb2mat(usa.nb, style="B")

##Write the 0-1 adjacency matrix
W.raw = usa.adj.mat
#Remove D.C. and Florida
W=W.raw[c(1:7,10:49),c(1:7,10:49)]
colnames(W) <- rownames(W)
D = diag(colSums(W))