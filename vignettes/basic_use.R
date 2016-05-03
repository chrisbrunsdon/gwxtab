## ---- fig.show='hold'----------------------------------------------------
## Demonstrator for SpatialCrossTabs

# Packages and data
library(gwxtab)
library(sp)
library(RColorBrewer)
data(classification)

## ---- fig.show='hold'----------------------------------------------------
# Set up a crosstab function
dummy_xtab <- new_spxt(loc_confuse,'actual','svmpred')

# Define a portmanteau accuracy function
portmanteau <- function(x) data.frame(port=sum(diag(x))/sum(x))

# Define two sample grids
hg <- spsample(loc_confuse,5000,'hexagonal',offset=c(0.5,0.5))
st <- spsample(loc_confuse,5000,'stratified')

# Work out local Portmanteau accuracy at sample points (bw=20km)

hg_port <- gwxtab_sample(hg,dummy_xtab,fixed(15),portmanteau)
st_port <- gwxtab_sample(st,dummy_xtab,fixed(15),portmanteau)

## ---- fig.show='hold'----------------------------------------------------
# ... now draw it

par(mar=c(0,0,0,0)+0.4)
plot(loc_confuse,pch=16,col='darkred')
plot(hg_port,add=T,col='navy',pch=16,cex=hg_port$port/2)
plot(loc_confuse,pch=16,col='darkred')
plot(st_port,add=T,col='navy',pch=16,cex=st_port$port/3)


## ------------------------------------------------------------------------
# Now use adaptive bandwidths
hg_port_a <- gwxtab_sample(hg,dummy_xtab,adapt(7.5),portmanteau)
st_port_a <- gwxtab_sample(st,dummy_xtab,adapt(7.5),portmanteau)

## ---- fig.show='hold'----------------------------------------------------
# ... now draw it

par(mar=c(0,0,0,0)+0.4)
plot(loc_confuse,pch=16,col='darkgreen')
plot(hg_port_a,add=T,col='tomato',pch=16,cex=hg_port_a$port/2)
plot(loc_confuse,pch=16,col='darkgreen')
plot(st_port_a,add=T,col='tomato',pch=16,cex=st_port_a$port/3)


## ------------------------------------------------------------------------
# Define a majority vote function
maj <- function(x) {
  votes <- rowSums(x)
  idx <- which.max(votes)
  conf <- max(votes)
  data.frame(choice = names(votes)[idx], conf=conf)
}
# Make a sample spdf
hg_maj_a <- gwxtab_sample(hg,dummy_xtab,adapt(10.6),maj)
st_maj_a <- gwxtab_sample(st,dummy_xtab,adapt(10.6),maj)

## ---- fig.show='hold'----------------------------------------------------
# ... now draw it

col_list <- brewer.pal(8,'Dark2')

par(mar=c(0,0,0,0)+0.4)
plot(loc_confuse,pch=16,col='darkgreen')
plot(hg_maj_a,add=T,col=col_list[hg_maj_a$choice],pch=16,cex=hg_maj_a$conf*0.3)
plot(loc_confuse,pch=16,col='darkgreen')
plot(st_maj_a,add=T,col=col_list[st_maj_a$choice],pch=16,cex=st_maj_a$conf*0.3)

