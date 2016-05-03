#' An S4 class to represent a bandwidth for gw-xtabs eyc.
#'
#' @slot value The value of the bandwidth
#' @slot units The units of length used to specify a non-adaptive bandwidth
#' @slot adaptive A logical value indicating whether bandwidth is adaptive
#' @export
bandwidth <- setClass(
  # Set the name for the class
  "bandwidth",

  # Define the slots
  slots = c(
    value = "numeric",
    units = "character",
    adaptive = "logical"
  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    value = 0,
    units = "km",
    adaptive = FALSE
  ),

  # Make a function that can test to see if the data is consistent.
  # This is not called if you have an initialize function defined!
  validity=function(object)
  { if (!is.numeric(object@value)) return('Non-numeric bw value')
    if (!is.character(object@units)) return('Non-character units.')
    if (! objects@units %in% c('km','m','mi','nm','deg','rad')) return('Unrecognised units - must be "km","mi","nm","deg","rad" or NA')
    if (!is.logical(object@adaptive)) return('Whether adaptive must be logical')
    return(TRUE)
  }
)

setMethod("initialize",
          "bandwidth",
          function(.Object, value = 0, units=as.character(NA), adaptive=FALSE) {
            .Object@value = value
            if (is.na(units) || units %in%  c('km','m','mi','nm','deg','rad')) {
              .Object@units = units
            } else {
              stop('Unrecognised units: must be "km","mi","m",nm","deg","rad" or NA')
            }
            .Object@adaptive = adaptive
            .Object
          })

#' Adaptive bandwidth
#'
#' Create an adaptive bandwidth object
#'
#' @param n - Number of nearest neighbours
#' @export
adapt <- function(n) new('bandwidth',value=n,adaptive=TRUE)

# Internal helper - conversions
units_list <- c('km','m','mi','nm','deg','rad')
cvt <- matrix(NA,length(units_list),length(units_list))
colnames(cvt) <- units_list
rownames(cvt) <- units_list
diag(cvt) <- 1.0
cvt['km','m'] <- 1000 ; cvt['m','km'] <- 1.0/cvt['km','m']
cvt['km','mi'] <- 0.621371; cvt['mi','km'] <- 1.0/cvt['km','mi']
cvt['km','nm'] <- 0.539957; cvt['nm','km'] <- 1.0/cvt['km','nm']
cvt['rad','km'] <- 6371.008; cvt['km','rad'] <- 1.0/cvt['rad','km']
cvt['deg','rad'] <- pi/180.0; cvt['rad','deg'] <- 1.0/cvt['deg','rad']
cvt['km','deg'] <- cvt['km','rad'] * cvt['rad','deg']; cvt['deg','km'] <- 1.0/cvt['km','deg']
cvt['m','mi'] <- cvt['km','mi']/1000; cvt['mi','m'] <- 1.0/cvt['m','mi']
cvt['m','nm'] <- cvt['km','nm']/1000; cvt['nm','m'] <- 1.0/cvt['m','nm']
cvt['m','deg'] <- cvt['km','deg']/1000; cvt['deg','m'] <- 1.0/cvt['m','deg']
cvt['m','rad'] <- cvt['km','rad']/1000; cvt['rad','m'] <- 1.0/cvt['m','rad']
cvt['deg','nm'] <- 60; cvt['nm','deg'] <- 1.0/cvt['deg','nm']
cvt['mi',c('nm','deg','rad')] <- cvt['km',c('nm','deg','rad')] * cvt['mi','km']
cvt[c('nm','deg','rad'),'mi']  <- 1.0/cvt['mi',c('nm','deg','rad')]
cvt['rad','nm'] <- 10800/pi; cvt['nm','rad'] <- pi/10800;

# Internal helper - find units from a proj4string
get_units <- function(x) {
  p4str <- proj4string(x)
  if (is.na(p4str)) return(as.character(NA))
  ustr <- regmatches(p4str,regexpr("units=\\w*",p4str))
  if (length(ustr) == 0) return('deg')
  return(sub("units=","",ustr))
}

#' Fixed bandwidth
#'
#' Create a fixed bandwidth object
#'
#' @param n - value (distance) of bandwidth
#' @param units - Distance units used ('km','m','mi','deg','rad','nm' or NA)
#' @export
fixed <- function(n,units='km') new('bandwidth',value=n,units=units)

#' @export
setGeneric('is.adaptive',
           function(bwo) standardGeneric('is.adaptive')
)

#' @export
setMethod("is.adaptive",
          "bandwidth",
          function(bwo) bwo@adaptive
          )

#' @export
setMethod("is.adaptive",
          "numeric",
          function(bwo) FALSE
)

#' @export
setMethod('show',
          'bandwidth',
          function(object) {
            if (is.adaptive(object)) {
              cat(paste0('Adaptive bandwidth - ',object@value,' nearest neighbours.'))
            } else {
              cat(paste0('Fixed bandwidth - ',object@value,object@units,'.'))
            }
          })


#' @export
setGeneric('units_used',
           function(x) standardGeneric('units_used')
)


#' @export
setMethod("units_used",
          "bandwidth",
          function(x) x@units
)


nnd <- function(X1,X2,n) {
  n_int <- n %/% 1
  n_fra <- n %% 1
  d1 <- RANN::nn2(X1,X2,n_int)$nn.dists[,n_int]
  d2 <- RANN::nn2(X1,X2,n_int+1)$nn.dists[,n_int+1]
  return(d1*(1-n_fra) + d2*n_fra)
}
