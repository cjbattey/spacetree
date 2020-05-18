library(ggplot2);library(plot3D);library(data.table);library(plotly)
setwd("~/sptree/")

dispersal_rate <- 3
bias <- 0

#spatial pedigrees with the highlighted path of a locus
pdf("fig/pedigree_w_lineage_path.pdf",width=6.5,height=2.5,useDingbats = F)
par(mfrow=c(1,2))
par(mar=rep(0.4,4))
segments3D(x0=c(0,50),y0=c(0,50),z0=c(0,7),bty="b2",zlab="T",phi = 1)
x <- c(runif(1,0,50));y <- c(runif(1,0,50))
x <- 10; y <- 10
samples <- data.frame(x,y,z=0,id=0,parent=-1)
old <- samples
lineage <- data.frame(x=0,y=0,z=0)[0,]
for(i in 1){
  nancestors <- 2^i
  newx <- rnorm(nancestors,old$x+bias,dispersal_rate)
  newy <- rnorm(nancestors,old$y+bias,dispersal_rate)
  newx[newx>50] <- 50-(newx[newx>50]-50)
  newy[newy>50] <- 50-(newy[newy>50]-50)
  newx[newx<0] <- -newx[newx<0]
  newy[newy<0] <- -newy[newy<0]
  new <- data.frame(x=newx,y=newy)
  for(j in 1:nrow(old)){
    segments3D(x0=old$x[j],x1=new$x[(j*2)-1],
               y0=old$y[j],y1=new$y[(j*2)-1],
               z0=i-1,z1=i,add=T,col="steelblue2",lwd=0.75,alpha.col=0.5)
    segments3D(x0=old$x[j],x1=new$x[(j*2)],
               y0=old$y[j],y1=new$y[(j*2)],
               z0=i-1,z1=i,add=T,col="steelblue2",lwd=0.75,alpha.col=0.5)
  }
  new$z <- i
  new$id <- (max(samples$id)+1):(max(samples$id)+nancestors)
  new$parent <- unlist(lapply(old$id,function(e) rep(e,2)))
  samples <- rbind(samples,new)
  points3D(x=samples$x,y=samples$y,z=samples$z,pch=16,col="black",cex=0.4,add=T,bty="b2",zlab="T",phi = 1)
  old <- new
}
meioses <- sapply(1:8,function(e) sample(1:2,1))
anc <- samples[0,];parent_id <- 0;pids <- parent_id
for(i in 1){
  gen <- subset(samples,z==i & parent==parent_id)
  parent_id <- gen$id[meioses[i]]
  pids <- append(pids,parent_id)
  anc <- rbind(anc,gen[meioses[i],])
  if(i==1){
    segments3D(x0=x,x1=anc$x[i],
               y0=y,y1=anc$y[i],
               z0=0,z1=anc$z[i],
               col="orangered",lwd=2,add=T)
  }
  if(i>1){
    segments3D(x0=anc$x[i-1],x1=anc$x[i],
               y0=anc$y[i-1],y1=anc$y[i],
               z0=anc$z[i-1],z1=anc$z[i],
               col="orangered",lwd=2,add=T)
  }
}
par(mar=c(2,2,2,5),mgp=c(3,0.3,0))
scatter2D(x=samples$x,y=samples$y,xlim=c(0,50),ylim=c(0,50),
          colvar=samples$id %in% pids,col=c("steelblue2","orangered"),pch=1,lwd=1.5,
          cex=0.75,cex.axis=0.8,cex.main=0.8,
          main="Ancestors",colkey=list(labels=c("genealogical","genetic"),at=c(0.25,0.75),cex.axis=0.8),
          tck=-0.02)
grid(col="grey",lty=1)
scatter2D(x=samples$x[samples$id %in% pids],y=samples$y[samples$id %in% pids],add=T,col="orangered",pch=16)

dev.off()



# 
# fig <- plot_ly(data=samples,type='scatter3d',x=~x,y=~y,z=~z*4)
# lines <- list()
# 
# ids <- unique(samples$id[samples$z<5])
# for(i in 1:length(ids)){
#   anc <- subset(samples,parent==ids[i]|id==ids[i])
#   lines[[i]] <- list("x0"=anc$x[1],"x1"=anc$x[2],"y0"=anc$y[1],"y1"=anc$y[2],"z0"=anc$z[1],"z1"=anc$z[2])
#   lines[[i]] <- list("x0"=anc$x[1],"x1"=anc$x[3],"y0"=anc$y[1],"y1"=anc$y[3],"z0"=anc$z[1],"z1"=anc$z[3])
# }
# 
# fig <- layout(fig, shapes = lines)
# fig



