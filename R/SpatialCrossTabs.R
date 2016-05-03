######################################################################
# Create the SpatialCrossTabs class
#
# This is used to represent an object field of contingency tables.
# The model is a 'melt-and-map' approach
# There is an mxn contingency table associated with each point (U,V)
# The contigency table can be summarised by a single number (eg odds ratio)
# or a 1-d vector (eg row proportions, column proportions) and converted to a
# SpatialPointsDataFrame - which may then be mapped.

# First, load up some c++ helpers
#

# Rcpp::sourceCpp('distkern.cpp')


# Worker tools,  not public
has_data <- function(sspdf) "data" %in% slotNames(sspdf)
guess_constructor <- function(x) {
  cl <- class(x)
  if (grepl("Points",cl)) return(SpatialPointsDataFrame)
  if (grepl("Polygons",cl)) return(SpatialPolygonsDataFrame)
  stop("Data must be SpatialPolygons or SpatialPoints")
}

#' gcr2km: Great circle radian angle to kilometers conversion
gcr2km <- 6371.008

#' gwxtab: A package for manipulating geographically weighted crosstabulations
#'
#' For more details see the help vignette:
#' \code{vignette("basic_use", package = "gwxtab")}
#'
#' @importClassesFrom sp Spatial CRS
#' @importMethodsFrom sp coordinates proj4string proj4string<- is.projected spTransform
#' @importFrom sp CRS
#' @importFrom RANN nn2
#' @import Rcpp
#' @docType package
#' @name gwxtab
#' @useDynLib gwxtab
#' @importClassesFrom  raster Raster
#' @export bisq bisq2 dist2
NULL

#' Spatially referenced land use data
#'
#' A SpatialPointsDataFrame named \code{loc_confuse} containing actual (\code{actual}) and predicted (\code{svmpred}) land use at
#' a number of data points. It has the following columns:
#' \tabular{rl}{
#'   \code{actual} \tab Actual land use class \cr
#'   \code{svmpred} \tab Predicted land use class using a support vector machine approach
#' }
#'
#' @name classification
#' @docType data
#' @keywords data
NULL

#' gcr2km: Great circle radian angle to kilometers conversion
#' @export
gcr2km <- 6371.008

#' gcr2mi: Great circle radian angle to miles conversion
#' @export
gcr2mi <- gcr2km * 0.621371

#' gcr2nm: Great circle radian angle to nautical miles conversion
#' @export
gcr2nm <- 10800 / pi

#' @include bandwidthclass.R

#' An S4 class to represent a geographically weighted cross-tabulation
#'
#' @slot locations A \code{c(n,2)} matrix giving coordinates of points at which cross-tabulation is evaluated
#' @slot xtabs A length \code{n} list of \code{c(l,m)} dimensional cross-tabulations
#' @slot proj4 A \code{CRS} object containing the map projection
#' @export
SpatialCrossTabs <- setClass(
  # Set the name for the class
  "SpatialCrossTabs",

  # Define the slots
  slots = c(
    locations = "matrix",
    xtabs = "list",
    proj4 = "CRS"
  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    locations = NULL,
    xtabs = NULL,
    proj4 = sp::CRS()
  ),

  # Make a function that can test to see if the data is consistent.
  # This is not called if you have an initialize function defined!
  validity=function(object)
  { if (is.null(object@locations) && is.null(object@xtabs)) return(TRUE)
    if (!is.matrix(object@locations)) return('Non numeric locations slot.')
    if (!is.list(object@xtabs)) return('Badly formed xtabs slot.')
    return(TRUE)
  }
)

## Set up a SpatialCrossTabs from a SpatialPointsDataFrame and a matrix

#' Create a new `SpatialCrossTabs` object
#'
#' @param spdf A spatial object or nx2 matrix of coordinates
#' @param var1 A 2-dimensional table object,  or a string specifying row variable
#' @param var2 A string specifying column variable (if \code{var1} specifies row variable)
#'
#' @name new_spxt
#' @docType methods
NULL

#' @rdname new_spxt
#' @export
setGeneric(name="new_spxt",
           def=function(spdf,var1,var2)
           {
             standardGeneric("new_spxt")
           }
)

#' @rdname new_spxt
#' @export
setMethod(f="new_spxt",
          signature=c("SpatialPointsDataFrame","table","missing"),
          definition=function(spdf,var1)
          {
            pts <- coordinates(spdf)
            tbl <- vector(mode='list',length=nrow(pts))
            for (i in 1:length(tbl)) tbl[[i]] <- var1
            theObject <- SpatialCrossTabs(locations=pts,xtabs=tbl)
            proj4string(theObject) <- CRS(proj4string(spdf))
            return(theObject)
          }
)

#' @rdname new_spxt
#' @export
setMethod(f="new_spxt",
          signature=c("matrix","table","missing"),
          definition=function(spdf,var1)
          {
            tbl <- vector(mode='list',length=nrow(spdf))
            for (i in 1:length(tbl)) tbl[[i]] <- var1
            theObject <- SpatialCrossTabs(locations=spdf,xtabs=tbl)
            return(theObject)
          }
)

#' @rdname new_spxt
#' @export
setMethod(f="new_spxt",
          signature=c("SpatialPointsDataFrame","character","character"),
          definition=function(spdf,var1,var2)
          {
            pts <- coordinates(spdf)
            tbl <- vector(mode='list',length=nrow(pts))
            inject <- table(spdf@data[,var1],spdf@data[,var2])
            names(dimnames(inject)) <- c(var1,var2)
            for (i in 1:length(tbl)) {
              inject <- inject * 0
              inject[as.character(spdf@data[i,var1]),as.character(spdf@data[i,var2])] <- 1
              tbl[[i]] <- inject
            }
            theObject <- SpatialCrossTabs(locations=pts,xtabs=tbl)
            proj4string(theObject) <- CRS(proj4string(spdf))
            return(theObject)
          }
)

#' Print out basic (short form) information for a SpatialCrossTabs object
#' @param object A SpatialCrossTabs object
#' @note  Not usually called explicitly
#' @export
setMethod(f='show',
          signature = 'SpatialCrossTabs',
          definition = function(object) {
            cat("Spatial Crosstabulation object\n")
            cat(sprintf("Number of locations: %d\n",length(object@xtabs)))
            cat(sprintf("Dimension: %d x %d\n",nrow(object@xtabs[[1]]),ncol(object@xtabs[[1]])))
            cat(sprintf("Proj4 String: %s\n",object@proj4@projargs))
          }
)

#' Access the nth crosstabulation in a SpatialCrossTabs object
#' @param x A SpatialCrossTabs object
#' @param i the index of the crosstabulation to extract
#' @return The \code{i}th crosstabulation in \code{x}.
setMethod(f='[',
          signature = c('SpatialCrossTabs','numeric'),
          definition = function(x,i) x@xtabs[[i]]
)

#' Assign the nth crosstabulation in a SpatialCrossTabs object
#' @param x A SpatialCrossTabs object
#' @param i the index of the crosstabulation to extract
#' @param A contingency table
setMethod(f='[<-',
          signature = c('SpatialCrossTabs','numeric','missing','table'),
          definition = function(x,i,value) {
            x@xtabs[[i]] <- value
            x
          }
)

#' Sum of contingency tables
#'
#' @param spxt A SpatialCrossTabs object
#' @return sum of all contingency tables in \code{spxt}
#' @name tabsum
NULL


#' @rdname tabsum
#' @export
setGeneric(name="tabsum",
           def=function(spxt)
           {
             standardGeneric("tabsum")
           }
)

#' @rdname tabsum
#' @export
setMethod(f='tabsum',
          signature = c('SpatialCrossTabs'),
          definition = function(spxt) {
            res <- spxt@xtabs[[1]]*0
            for (i in 1:length(spxt@xtabs)) res <- res + spxt@xtabs[[i]]
            res
          }
)

#' Return coordinates of a \code{SpatialCrossTabs} object
#' @export
#' @param obj A \code{SpatialCrossTabs} object
#' @return The coordinates as a \code{c(n,2)} matrix
setMethod(f='coordinates',
          signature = c('SpatialCrossTabs'),
          definition = function(obj,...) {
            obj@locations
          }
)


#' Weighted sum of contingency tables
#'
#' @param spxt A SpatialCrossTabs object
#' @param w A set of weights (optional - default is \code{1/n})
#' @return Weighted sum of all contingency tables in \code{spxt}
#' @name wtabsum
NULL

#' @rdname wtabsum
#' @export
setGeneric(name="wtabsum",
           def=function(spxt,w)
           {
             standardGeneric("wtabsum")
           }
)

#' @rdname wtabsum
#' @export
setMethod(f='wtabsum',
          signature = c('SpatialCrossTabs','numeric'),
          definition = function(spxt,w) {
            res <- spxt@xtabs[[1]]*0
            if (length(w) != length(spxt@xtabs)) stop('Weights must be same length as number of xtabs.')
            for (i in 1:length(spxt@xtabs)) res <- res + spxt@xtabs[[i]]*w[i]
            res
          }
)

#' @rdname wtabsum
#' @export
setMethod(f='wtabsum',
          signature = c('SpatialCrossTabs','missing'),
          definition = function(spxt) {
            res <- spxt@xtabs[[1]]*0
            n <- length(spxt@xtabs)
            for (i in 1:length(spxt@xtabs)) res <- res + spxt@xtabs[[i]]/n
            res
          }
)

#' Geographically Weighted contingency tables
#'
#' @param spxt A SpatialCrossTabs object
#' @param xtab_pts A \code{SpatialPoints} object (optional - if omitted points are as in \code{spxt})
#' @param bw Bandwidth
#' @return Weighted sum of all contingency tables in \code{spxt}
#' @name gwxtab
NULL

#' @rdname gwxtab
#' @export
setGeneric(name="gwxtab",
           def=function(spxt,xtab_pts,bw)
           {
             standardGeneric("gwxtab")
           }
)

#' @rdname gwxtab
#' @export
setMethod(f='gwxtab',
          signature = c('SpatialCrossTabs','missing','numeric'),
          definition = function(spxt,bw) {
            mat <- bisq2(spxt@locations,spxt@locations,bw)
            n <- length(spxt@xtabs)
            gwtabs <- vector(mode='list',length=n)
            for (i in 1:n) gwtabs[[i]] <- wtabsum(spxt,mat[,i])
            SpatialCrossTabs(locations=spxt@locations,xtabs=gwtabs)
          }
)

#' @rdname gwxtab
#' @export
setMethod(f='gwxtab',
          signature = c('SpatialCrossTabs','SpatialPoints','numeric'),
          definition = function(spxt,xtab_pts,bw) {
            mat <- bisq2(spxt@locations,coordinates(xtab_pts),bw)
            n <- length(xtab_pts)
            gwtabs <- vector(mode='list',length=n)
            for (i in 1:n) gwtabs[[i]] <- wtabsum(spxt,mat[,i])
            SpatialCrossTabs(locations=coordinates(xtab_pts),xtabs=gwtabs)
          }
)

#' Geographically Weighted contingency table probing functions
#'
#' @param spxt A SpatialCrossTabs object
#' @param bw Bandwidth
#' @param mode Kind of parameters to give to probe: \code{scalar} for \code{(x,y)} pair, \code{vector} for \code{X} length 2 vector, \code{matrix} for \code{c(n,2)} matrix, or \code{Spatial} for a a spatial object
#' @param melt 'melting' operator
#' @return Function to return geographically weighted sum of all contingency tables in \code{spxt} centred on probe argument.
#'
#' @name gwxtab_probe
NULL

#' @rdname gwxtab_probe
#' @export
setGeneric(name="gwxtab_probe",
           def=function(spxt,bw,mode,melt)
           {
             standardGeneric("gwxtab_probe")
           }
)

#' @rdname gwxtab_probe
#' @export
setMethod(f='gwxtab_probe',
          signature = c('SpatialCrossTabs','numeric','missing','missing'),
          definition = function(spxt,bw) {
            function (x,y) {
              wt <- bisq(spxt@locations,c(x,y),bw)
              wtabsum(spxt,wt)
            }
          }
)

#' @rdname gwxtab_probe
#' @export
setMethod(f='gwxtab_probe',
          signature = c('SpatialCrossTabs','bandwidth','missing','missing'),
          definition = function(spxt,bw) {
            if (!is.adaptive(bw)) {
              bw2 <- bw@value * cvt[units_used(bw),units_used(spxt)]
              function (x,y) {
                wt <- bisq(spxt@locations,c(x,y),bw2)
                wtabsum(spxt,wt)
              }
            } else {
              nn <- bw@value
              function (x,y) {
                xx <- matrix(c(x,y),1,2)
                bw2 <- nnd(spxt@locations,xx,nn)
                wt <- bisq(spxt@locations,c(x,y),bw2)
                wtabsum(spxt,wt)
              }
            }
          }
)

#' @rdname gwxtab_probe
#' @export
setMethod(f='gwxtab_probe',
          signature = c('SpatialCrossTabs','bandwidth','missing','function'),
          definition = function(spxt,bw,melt) {
            if (!is.adaptive(bw)) {
              bw2 <- bw@value * cvt[units_used(bw),units_used(spxt)]
              function (x,y) {
                wt <- bisq(spxt@locations,c(x,y),bw2)
                melt(wtabsum(spxt,wt))
              }
            } else {
              nn <- bw@value
              function (x,y) {
                xx <- matrix(c(x,y),1,2)
                bw2 <- nnd(spxt@locations,xx,nn)
                wt <- bisq(spxt@locations,c(x,y),bw2)
                melt(wtabsum(spxt,wt))
              }
            }
          }
)


#' @rdname gwxtab_probe
#' @export
setMethod(f='gwxtab_probe',
          signature = c('SpatialCrossTabs','numeric','character','missing'),
          definition = function(spxt,bw,mode) {
            switch(mode,
              scalar=function (x,y) {
                wt <- bisq(spxt@locations,c(x,y),bw)
                wtabsum(spxt,wt)
              },
              vector = function (x) {
                wt <- bisq(spxt@locations,x,bw)
                wtabsum(spxt,wt)
              },
              matrix = function (x) {
                n <- nrow(x)
                res <- array(0,c(dim(spxt@xtabs[[1]]),n))
                for (i in 1:n) {
                  wt <- bisq(spxt@locations,x[i,],bw)
                  res[,,i] <- wtabsum(spxt,wt)
                }
                res
              },
              Spatial = function(x) {
                xx <- coordinates(x)
                n <- nrow(xx)
                res <- vector(mode='list',length=n)
                if (! is.projected(x)) stop('The Spatial objected must be projected')
                for (i in 1:n) {
                  wt <- bisq(spxt@locations,x[i,],bw)
                  res[[i]] <- melt(wtabsum(spxt,wt))
                }
                cstr <- guess_constructor(x)
                dfres <- data.frame(do.call(rbind,res))
                cstr(x,dfres,match.ID=FALSE)
              }
            )
          }
)


#' @rdname gwxtab_probe
#' @export
setMethod(f='gwxtab_probe',
          signature = c('SpatialCrossTabs','numeric','character','function'),
          definition = function(spxt,bw,mode,melt) {
            switch(mode,
                   scalar=function (x,y) {
                     wt <- bisq(spxt@locations,c(x,y),bw)
                     melt(wtabsum(spxt,wt))
                   },
                   vector = function (x) {
                     wt <- bisq(spxt@locations,x,bw)
                     melt(wtabsum(spxt,wt))
                   },
                   matrix = function (x) {
                     n <- nrow(x)
                     res <- vector(mode='list',length=n)
                     for (i in 1:n) {
                       wt <- bisq(spxt@locations,x[i,],bw)
                       res[[i]] <- melt(wtabsum(spxt,wt))
                     }
                     res
                   },
                   Spatial = function(x) {
                     xx <- coordinates(x)
                     n <- nrow(xx)
                     res <- vector(mode='list',length=n)
                     if (! is.projected(x)) stop('The Spatial objected must be projected')
                     for (i in 1:n) {
                       wt <- bisq(spxt@locations,xx[i,],bw)
                       res[[i]] <- melt(wtabsum(spxt,wt))
                     }
                     cstr <- guess_constructor(x)
                     dfres <- data.frame(do.call(rbind,res))
                     cstr(x,dfres,match.ID=FALSE)
                   }
            )
          }
)

#' @rdname gwxtab_probe
#' @export
setMethod(f='gwxtab_probe',
          signature = c('SpatialCrossTabs','bandwidth','character','function'),
          definition = function(spxt,bw,mode,melt) {
            if (!is.adaptive(bw)) {
              bw2 <- bw@value * cvt[units_used(bw),units_used(spxt)]
              switch(mode,
                   scalar=function (x,y) {
                     wt <- bisq(spxt@locations,c(x,y),bw2)
                     melt(wtabsum(spxt,wt))
                   },
                   vector = function (x) {
                     wt <- bisq(spxt@locations,x,bw2)
                     melt(wtabsum(spxt,wt))
                   },
                   matrix = function (x) {
                     n <- nrow(x)
                     res <- vector(mode='list',length=n)
                     for (i in 1:n) {
                       wt <- bisq(spxt@locations,x[i,],bw2)
                       res[[i]] <- melt(wtabsum(spxt,wt))
                     }
                     res
                   },
                   Spatial = function(x) {
                     xx <- coordinates(x)
                     n <- nrow(xx)
                     res <- vector(mode='list',length=n)
                     if (! is.projected(x)) stop('The Spatial objected must be projected')
                     for (i in 1:n) {
                       wt <- bisq(spxt@locations,xx[i,],bw2)
                       res[[i]] <- melt(wtabsum(spxt,wt))
                     }
                     cstr <- guess_constructor(x)
                     dfres <- data.frame(do.call(rbind,res))
                     cstr(x,dfres,match.ID=FALSE)
                   }
              )
            } else {
              nn <- bw@value
              switch(mode,
                     scalar=function (x,y) {
                       xx <- matrix(c(x,y),1,2)
                       bw2 <- nnd(spxt@locations,xx,nn)
                       wt <- bisq(spxt@locations,c(x,y),bw2)
                       melt(wtabsum(spxt,wt))
                     },
                     vector = function (x) {
                       xx <- matrix(xx,1,2)
                       bw2 <- nnd(spxt@locations,xx,nn)
                       wt <- bisq(spxt@locations,x,bw2)
                       melt(wtabsum(spxt,wt))
                     },
                     matrix = function (x) {
                       n <- nrow(x)
                       res <- vector(mode='list',length=n)
                       bw2 <- nnd(spxt@locations,x,nn)
                       for (i in 1:n) {
                         wt <- bisq(spxt@locations,x[i,],bw2[i])
                         res[[i]] <- melt(wtabsum(spxt,wt))
                       }
                       res
                     },
                     Spatial = function(x) {
                       xx <- coordinates(x)
                       n <- nrow(xx)
                       res <- vector(mode='list',length=n)
                       bw2 <- nnd(spxt@locations,xx,nn)
                       if (! is.projected(x)) stop('The Spatial objected must be projected')
                       for (i in 1:n) {
                         wt <- bisq(spxt@locations,xx[i,],bw2[i])
                         res[[i]] <- melt(wtabsum(spxt,wt))
                       }
                       cstr <- guess_constructor(x)
                       dfres <- data.frame(do.call(rbind,res))
                       cstr(x,dfres,match.ID=FALSE)
                     }
              )
            }
          }
)



#' Return \code{proj4} string for map of projection
#' @export
setMethod(f=proj4string,
          signature='SpatialCrossTabs',
          definition = function(obj) obj@proj4@projargs
            )


#' Set \code{proj4} string for map of projection
#' @name set_proj4string
NULL

#' @export
#' @rdname set_proj4string
setMethod(f='proj4string<-',
          signature= c('SpatialCrossTabs','CRS'),
          definition = function(obj,value) {
            obj@proj4 <- value
            obj
            }
)


#' @export
#' @rdname set_proj4string
setMethod(f='proj4string<-',
          signature= c('SpatialCrossTabs','character'),
          definition = function(obj,value) {
            obj@proj4 <- CRS(value)
            obj
          }
)

#' @export
#' @rdname set_proj4string




#' Check whether SpatialCrossTabs object uses projected  coordinates
#'
#' This is exactly the same as the \code{is.projected} method for \code{Spatial} objects
#' @export
setMethod(f='is.projected',
          signature= c('SpatialCrossTabs'),
          definition = function (obj)
          {
            p4str <- proj4string(obj)
            if (is.na(p4str) || !nzchar(p4str))
              return(as.logical(NA))
            else {
              res <- grep("longlat", p4str, fixed = TRUE)
              if (length(res) == 0)
                return(TRUE)
              else return(FALSE)
            }
          }
)


#' Check whether Raster object uses projected  coordinates
#'
#' This is exactly the same as the \code{is.projected} method for \code{Spatial} objects
#' @export
setMethod(f='is.projected',
          signature= c('Raster'),
          definition = function (obj)
          {
            p4str <- proj4string(obj)
            if (is.na(p4str) || !nzchar(p4str))
              return(as.logical(NA))
            else {
              res <- grep("longlat", p4str, fixed = TRUE)
              if (length(res) == 0)
                return(TRUE)
              else return(FALSE)
            }
          }
)



#' spTransform method for \code{SpatialCrossTabs}
#'
#' Changes the map projection for a \code{SpatialCrossTabs} object.
#' @param x The \code{SpatialCrossTabs} object.
#' @param CRSobj Either a \code{CRS} object, a \code{PROJ4} string or a \code{Spatial} object to extract the \code{PROJ4} string from.
#' @return A \code{SpatialCrossTabs} object with transformed coordinates.
#'
#' @name projection_transform
NULL


#' spTransform method for \code{SpatialCrossTabs}
#' @export
setMethod(f='spTransform',
          signature= c('SpatialCrossTabs','CRS'),
          definition = function (x,CRSobj)
          {
            p4str <- proj4string(x)
            spdf <- SpatialPoints(coordinates(x),proj4string = CRS(p4str))
            spdf <- sp::spTransform(spdf,CRSobj)
            x2 <- x
            x2@locations <- coordinates(spdf)
            x2@proj4 <- CRSobj
            return(x2)
          }
)


#' @rdname projection_transform
#' @export
setMethod(f='spTransform',
          signature= c('SpatialCrossTabs','CRS'),
          definition = function (x,CRSobj)
          {
            p4str <- proj4string(x)
            spdf <- SpatialPoints(coordinates(x),proj4string = CRS(p4str))
            spdf <- sp::spTransform(spdf,CRSobj)
            x2 <- x
            x2@locations <- coordinates(spdf)
            x2@proj4 <- CRSobj
            return(x2)
          }
)

#' @rdname projection_transform
#' @export
setMethod(f='spTransform',
          signature= c('SpatialCrossTabs','character'),
          definition = function (x,CRSobj)
          {
            p4str <- proj4string(x)
            spdf <- SpatialPoints(coordinates(x),proj4string = CRS(p4str))
            spdf <- sp::spTransform(spdf,CRS(CRSobj))
            x2 <- x
            x2@locations <- coordinates(spdf)
            x2@proj4 <- CRS(CRSobj)
            return(x2)
          }
)

#' @rdname projection_transform
#' @export
setMethod(f='spTransform',
          signature= c('SpatialCrossTabs','Spatial'),
          definition = function (x,CRSobj)
          {
            p4str <- proj4string(x)
            spdf <- SpatialPoints(coordinates(x),proj4string = CRS(p4str))
            spdf <- sp::spTransform(spdf,CRS(proj4string(CRSobj)))
            x2 <- x
            x2@locations <- coordinates(spdf)
            x2@proj4 <- CRS(proj4string(CRSobj))
            return(x2)
          }
)


#' Geographically Weighted contingency table sampling functions
#'
#' @param spxt A SpatialCrossTabs object
#' @param bw Bandwidth
#' @param melt 'melting' operator
#' @return Function to return geographically weighted sum of all contingency tables in \code{spxt} centred on probe argument.
#'
#' @name gwxtab_sample
NULL

#' @rdname gwxtab_sample
#' @export
setGeneric(name="gwxtab_sample",
           def=function(spdf,spxt,bw,melt)
           {
             standardGeneric("gwxtab_sample")
           }
)

#' @rdname gwxtab_sample
#' @export
setMethod(f='gwxtab_sample',
          signature = c('Spatial','SpatialCrossTabs','numeric','function'),
          definition = function(spdf,spxt,bw,melt) gwxtab_probe(spxt,bw,'Spatial',melt)(spdf)
)

#' @rdname gwxtab_sample
#' @export
setMethod(f='gwxtab_sample',
          signature = c('Spatial','SpatialCrossTabs','bandwidth','function'),
          definition = function(spdf,spxt,bw,melt) gwxtab_probe(spxt,bw,'Spatial',melt)(spdf)
)


#' @rdname gwxtab_sample
#' @export
setMethod(f='gwxtab_sample',
          signature = c('Raster','SpatialCrossTabs','bandwidth','function'),
          definition = function(spdf,spxt,bw,melt) {
            X <- coordinates(spdf)
            vals <- do.call(c,gwxtab_probe(spxt,bw,'matrix',melt)(X))
            new_spdf <- spdf
            values(new_spdf) <- vals
            return(new_spdf)
          }
)

#' @export
setMethod("units_used",
          "SpatialCrossTabs",
          function(x) get_units(x)
)

#' @export
setMethod("units_used",
          "SpatialPointsDataFrame",
          function(x) get_units(x)
)

#' @export
setMethod("units_used",
          "SpatialPoints",
          function(x) get_units(x)
)

#' @export
setMethod("units_used",
          "SpatialPolygonsDataFrame",
          function(x) get_units(x)
)

#' @export
setMethod("units_used",
          "SpatialPolygons",
          function(x) get_units(x)
)

#' @export
setMethod("units_used",
          "Raster",
          function(x) get_units(x)
)

