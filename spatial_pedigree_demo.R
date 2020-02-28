library(ggplot2);library(plot3D);library(data.table)
setwd("~/Dropbox/job_talk/")

reflect <- function(x,y,xlim,ylim){
  x[x>xlim[2]] <- xlim[2]-(x[x>xlim[2]]-xlim[2])
  y[y>ylim[2]] <- ylim[2]-(y[y>ylim[2]]-ylim[2])
  x[x<xlim[1]] <- xlim[1]-x[x<xlim[1]]
  y[y<ylim[1]] <- ylim[1]-y[y<ylim[1]]
  return(c(x,y))
}

#paths of a pedigree through time
pdf("spatial_pedigree/5sample.pdf",width=5,height=5,useDingbats = F)
segments3D(x0=c(0,20),y0=c(0,20),z0=c(0,5),bty="b2",zlab="T",phi = 1)
x <- c(runif(1,0,20));y <- c(runif(1,0,20))
old <- data.frame(x,y)
samples <- data.frame(x,y,z=0)
for(i in 1:5){
  nancestors <- 2^i
  newx <- rnorm(nancestors,old$x,2.5)
  newy <- rnorm(nancestors,old$y,2.5)
  newx[newx>20] <- 20-(newx[newx>20]-20)
  newy[newy>20] <- 20-(newy[newy>20]-20)
  newx[newx<0] <- -newx[newx<0]
  newy[newy<0] <- -newy[newy<0]
  new <- data.frame(x=newx,y=newy)
  for(j in 1:nrow(old)){
    segments3D(x0=old$x[j],x1=new$x[(j*2)-1],
               y0=old$y[j],y1=new$y[(j*2)-1],
               z0=i-1,z1=i,add=T,col="forestgreen",lwd=0.5)
    segments3D(x0=old$x[j],x1=new$x[(j*2)],
               y0=old$y[j],y1=new$y[(j*2)],
               z0=i-1,z1=i,add=T,col="forestgreen",lwd=0.5)
  }
  old <- new
  new$z <- i
  samples <- rbind(samples,new)
  points3D(x=samples$x,y=samples$y,z=samples$z,pch=16,col="black",cex=0.5,add=T,bty="b2",zlab="T",phi = 1)
}
dev.off()

#spatial pedigrees with the highlighted path of a locus
pdf("spatial_pedigree/pedigree_w_lineage_path.pdf",width=5,height=5,useDingbats = F)
segments3D(x0=c(0,20),y0=c(0,20),z0=c(0,5),bty="b2",zlab="T",phi = 1)
x <- c(runif(1,0,20));y <- c(runif(1,0,20))
old <- data.frame(x,y)
samples <- data.frame(x,y,z=0)
lineage <- data.frame(x=0,y=0,z=0)[0,]
for(i in 1:5){
  nancestors <- 2^i
  meiosis <- sample(1:2,1)
  newx <- rnorm(nancestors,old$x,2.5)
  newy <- rnorm(nancestors,old$y,2.5)
  newx[newx>20] <- 20-(newx[newx>20]-20)
  newy[newy>20] <- 20-(newy[newy>20]-20)
  newx[newx<0] <- -newx[newx<0]
  newy[newy<0] <- -newy[newy<0]
  new <- data.frame(x=newx,y=newy)
  for(j in 1:nrow(old)){
    segments3D(x0=old$x[j],x1=new$x[(j*2)-1],
               y0=old$y[j],y1=new$y[(j*2)-1],
               z0=i-1,z1=i,add=T,col="cornflowerblue",lwd=0.5)
    segments3D(x0=old$x[j],x1=new$x[(j*2)],
               y0=old$y[j],y1=new$y[(j*2)],
               z0=i-1,z1=i,add=T,col="cornflowerblue",lwd=0.5)
  }
  old <- new
  new$z <- i
  samples <- rbind(samples,new)
  points3D(x=samples$x,y=samples$y,z=samples$z,pch=16,col="black",cex=0.4,add=T,bty="b2",zlab="T",phi = 1)
}
meioses <- sapply(1:6,function(e) sample(1:2,1))
anc <- samples[0,]
for(i in 1:6){
  gen <- subset(samples,z==unique(samples$z)[i])
  anc <- rbind(anc,gen[meioses[i],])
  if(i>1){
    segments3D(x0=anc$x[i-1],x1=anc$x[i],
               y0=anc$y[i-1],y1=anc$y[i],
               z0=anc$z[i-1],z1=anc$z[i],
               col="orange",lwd=2,add=T)
  }
}
dev.off()


#path of a single locus through time
x <- c(0);y <- c(0)
for(i in 1:100){
  x[i+1] <- rnorm(1,x[i],1)
  y[i+1] <- rnorm(1,y[i],1)
}
samples <- data.frame(x,y)
samples$z <- 1:nrow(samples)

ggplot(data=samples,aes(x=x,y=y))+
  geom_point(size=6)+
  theme(panel.grid.major=element_line(color="grey"))+
  transition_manual(z)
anim_save("dispersal.gif")


segments3D(x0=c(-20,20),y0=c(-20,20),z0=c(0,max(samples$z)),bty="b2",zlab="T",phi = 1)
for(i in 1:(max(samples$z)-1)){
  segments3D(x0=samples$x[i],x1=samples$x[i+1],
             y0=samples$y[i],y1=samples$y[i+1],
             z0=samples$z[i],z1=samples$z[i+1],add=T,col='cornflowerblue')
}
points3D(x=samples$x,y=samples$y,z=samples$z,pch=16,col="black",cex=0.3,add=T,bty="b2",zlab="T",phi = 1)

segments3D(x0=c(0,100),y0=c(0,100),z0=c(0,max(samples$z)),bty="b2",zlab="T",phi = 1)
for(i in 0:(max(samples$z))){
  for(j in letters[1:nsamples]){
    b <- samples[grepl(as.character(j),samples$sample),]
    b <- b[b$z %in% c(i,i+1),]
    segments3D(x0=b$x[1],x1=b$x[2],y0=b$y[1],y1=b$y[2],z0=b$z[1],z1=b$z[2],add=T,col="steelblue3")
    points3D(x=b$x[2],y=b$y[2],z=b$z[2],add=T,cex=0.5)
  }
}
