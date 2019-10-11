#spacetree
library(plot3D);library(sp);library(ape)

nsamples <- 5;dispersal <- 20;maxgen <- 40
x <- runif(nsamples,0,100)
y <- runif(nsamples,0,100)
z <- rep(0,nsamples)
samples <- data.frame(x,y,z,sample=letters[1:nsamples],stringsAsFactors = F)
testsamples <-  samples

s2 <- samples
coalsamples <- c()
node_ages <- rep(0,nsamples)
while(length(unique(s2$sample))>1 & max(samples$z)<=maxgen){
  #dispersal
  a <- data.frame(x=NA,y=NA,z=NA)[0,]
  for(s in unique(s2$sample)){
    e <- subset(s2,sample==s)
    newx <- e$x[1]+rnorm(1,mean=0,sd=dispersal)
    if(newx<0) newx <- abs(newx) else if(newx>100) newx <- 100-(newx-100) #reflecting boundaries
    newy <- e$y[1]+rnorm(1,mean=0,sd=dispersal)
    if(newy<0) newy <- abs(newx) else if(newy>100) newy <- 100-(newy-100)
    newz <- e$z[1]+1
    a <- rbind(a,data.frame(x=rep(newx,nrow(e)),y=rep(newy,nrow(e)),z=rep(newz,nrow(e))))
    if(nrow(e)>1) coalsamples <- e$sample
  }
  
  #mate selection
  a$sample <- s2$sample
  dists <- spDists(as.matrix(a[,c("x","y")]))
  if(any(dists[dists>0]<dispersal)){ #one pair of lineages within `dispersal` coalesce (could set probability instead but would take too long for plotting)
    pairind <- sample(which(dists>0 & dists<dispersal),1)
    mate1 <- a$sample[pairind%%nrow(a)]
    if(pairind%%nrow(a)==0) mate1 <- a$sample[nrow(a)]
    mate2 <- a$sample[ceiling(pairind/nrow(a))]
    bl1 <- a$z[1]-node_ages[a$sample==mate1]
    bl2 <- a$z[1]-node_ages[a$sample==mate2]
    node_ages[a$sample==mate1] <- a$z[1]
    node_ages <- node_ages[-which(a$sample==mate2)] 
    a$sample[a$sample==mate1] <- paste0("(",mate1,":",bl1,",",
                                            mate2,":",bl2,")")
    a <- a[!a$sample==mate2,]
  }
  s2 <- a
  samples <- rbind(samples,a)
  # print(samples)
  # print(node_ages)
}
for(i in 1:5){ #extra generations to show drift in the ancestral lineage
  #dispersal
  a <- data.frame(x=NA,y=NA,z=NA)[0,]
  for(s in unique(s2$sample)){
    e <- subset(s2,sample==s)
    newx <- e$x[1]+rnorm(1,mean=0,sd=dispersal)
    if(newx<0) newx <- abs(newx) else if(newx>100) newx <- 100-(newx-100) #reflecting boundaries
    newy <- e$y[1]+rnorm(1,mean=0,sd=dispersal)
    if(newy<0) newy <- abs(newx) else if(newy>100) newy <- 100-(newy-100)
    newz <- e$z[1]+1
    a <- rbind(a,data.frame(x=rep(newx,nrow(e)),y=rep(newy,nrow(e)),z=rep(newz,nrow(e))))
    if(nrow(e)>1) coalsamples <- e$sample
  }
  
  #mate selection
  a$sample <- s2$sample
  dists <- spDists(as.matrix(a[,c("x","y")]))
  if(any(dists[dists>0]<dispersal)){ #one pair of lineages within `dispersal` coalesce (could set probability instead but would take too long for plotting)
    pairind <- sample(which(dists>0 & dists<dispersal),1)
    mate1 <- a$sample[pairind%%nrow(a)]
    if(pairind%%nrow(a)==0) mate1 <- a$sample[nrow(a)]
    mate2 <- a$sample[ceiling(pairind/nrow(a))]
    bl1 <- a$z[1]-node_ages[a$sample==mate1]
    bl2 <- a$z[1]-node_ages[a$sample==mate2]
    node_ages[a$sample==mate1] <- a$z[1]
    node_ages <- node_ages[-which(a$sample==mate2)] 
    a$sample[a$sample==mate1] <- paste0("(",mate1,":",bl1,",",
                                        mate2,":",bl2,")")
    a <- a[!a$sample==mate2,]
  }
  s2 <- a
  samples <- rbind(samples,a)
}
tail(samples,n=1)
tree <- paste0(samples$sample[nrow(samples)],";")
tree <- read.tree(text=tree)


max(samples$z)
layout(matrix(c(1, 1, 2, 2, 2,
                1, 1, 2, 2, 2), nrow=2, byrow=TRUE))
par(xpd=T)
plot(tree,cex=1.5)
segments3D(x0=c(0,100),y0=c(0,100),z0=c(0,max(samples$z)),bty="b2",zlab="T",phi = 1)
points3D(x=samples$x[samples$z==0],y=samples$y[samples$z==0],z=samples$z[samples$z==0],pch=16,col="black",cex=0.75,add=T)
text3D(x=samples$x[samples$z==0]+5,y=samples$y[samples$z==0],z=samples$z[samples$z==0],labels=letters[1:nsamples],pch=16,col="black",cex=1.25,add=T)
for(i in 0:(max(samples$z))){
  for(j in letters[1:nsamples]){
    b <- samples[grepl(as.character(j),samples$sample),]
    b <- b[b$z %in% c(i,i+1),]
    segments3D(x0=b$x[1],x1=b$x[2],y0=b$y[1],y1=b$y[2],z0=b$z[1],z1=b$z[2],add=T,col="steelblue3")
    points3D(x=b$x[2],y=b$y[2],z=b$z[2],add=T,cex=0.5)
  }
}


