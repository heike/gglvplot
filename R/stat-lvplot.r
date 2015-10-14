#' Compute letter value summary table.
#'
#' @param x numeric vector
#' @param qu quantiles to compute
#' @param out binary vector of outliers (\code{TRUE} for outlier,
#'   \code{FALSE} otherwise)
#' @inheritParams determineDepth
#' @keywords internal
#' @return
#'  \item{letter.val}{letter value statistic, distinguishes between upper and
#'    lower LV statistic for all statistics but the median}
#'  \item{conf.int}{confidence interval of corresponding letter value
#'    statistic}
#'  \item{out}{list of defined outliers}
outputLVplot <- function(x,qu,k,out,alpha) {
  n <- length(x)

  depth <- getDepth(k,n)

  LV <- cbind(depth,lower=qu[k:1],upper=qu[k-1+1:k])
  conf <- confintLV(x, k, alpha=alpha)

  dimnames <- nameLV(k)
  row.names(LV) <- dimnames[[1]]
  row.names(conf) <- dimnames[[2]]

  result <- list(letter.val = LV, conf.int= conf,outliers = which(out))
  return(result)
}

# Get a vector of data depths for k letter values
getDepth <- function(k, n) {
  # compute letter values based on depth
  depth <- rep(0,k)
  depth[1] <- (1+n)/2

  if (k > 1) {
    for (j in 2:k) depth[j] <- (1 + floor(depth[j - 1])) / 2
  }
  depth
}

# Calculate first k letter values for vector x
calcLV <- function(x, k) {
  n <- length(x)
  depth <- getDepth(k, n)

  y <- sort(x)
  d <- c(rev(depth),n-depth+1)
  qu <- (y[floor(d)] + y[ceiling(d)])/2
  # floor and ceiling is the same for .0 values
  # .5 values yield average of two neighbours

  # k, k+1 is the median
  # report only once
  qu[-k]
}

# Determine names of first k letter values
# output is (1) list of names starting with 'M'edian, and
# (2) a vector of letter values ordered by rank from lower kth letter value to upper k letter value
nameLV <- function(k) {
  # list of
  #  k letter values (starting with median)
  #  lower/upper letter values ordered from lowest to highest
  lvs <- NULL
  conf <- "M"
  if (k > 1) {
    lvs <- c(LETTERS[6:1], LETTERS[c(26:14, 12:7)])[1:(k - 1)]
    conf <- c(paste(rev(lvs), "l", sep = ""), "M", paste(lvs, "u", sep = ""))
  }
  list(LV = c("M", lvs), conf = conf)
}


#' Determine depth of letter values needed for n observations.
#'
#' @details Supply one of \code{k}, \code{alpha} or \code{perc}.
#'
#' @param n number of observation to be shown in the LV boxplot
#' @param k number of letter value statistics used
#' @param alpha if supplied, depth k is calculated such that (1-\code{alpha})100% confidence
#'   intervals of an LV statistic do not extend into
#'   neighboring LV statistics.
#' @param perc if supplied, depth k is adjusted such that \code{perc} percent
#'   outliers are shown
determineDepth <- function(n, k = NULL, alpha = NULL, perc = NULL) {
  if (!is.null(k)) {
    stopifnot(is.numeric(k) && length(k) == 1)
    k <- as.integer(k)
  } else if (!is.null(perc)) {
    # we're aiming for perc percent of outlying points
    stopifnot(is.numeric(perc) && length(perc) == 1)

    k <- ceiling((log2(n))+1) - ceiling((log2(n*perc*0.01))+1) + 1
  } else if (!is.null(alpha)) {
    # confidence intervals around an LV statistic
    # should not extend into surrounding LV statistics
    stopifnot(is.numeric(alpha) && length(alpha) == 1)
    stopifnot(alpha > 0 && alpha < 1)

    #    cat(sprintf("two rules: %d (ceiling) %d (floor)", ceiling((log2(n))-log2(2*qnorm(alpha+(1-alpha)/2)^2)),
    #            floor(log2(n)) - floor(log2(2*qnorm(1-(1-alpha)/2)^2))))
    k <- floor(log2(n)) - floor(log2(2*qnorm(1-(1-alpha)/2)^2))
  } else {
    stop("Must specify one of k, alpha, perc", call. = FALSE)
  }

  max(k, 1L)
}

#' Compute table of k letter values for vector x
#'
#' Calculate k letter values from numeric vector x together with an 100(1-alpha)\% confidence interval.
#' @param x input numeric vector
#' @param k number of letter values to compute
#' @param alpha alpha-threshold for confidence level
#' @return data frame of letter value (upper and lower value, if not the median), depth, and confidence interval.
lvtable <- function(x, k, alpha=0.95) {
  n <- length(x)
  if (2^k > n) k <- ceiling(log2(n)) + 1

  # depths for letter values
  depth <- getDepth(k, n)

  # letter value
  qu <- calcLV(x,k)

  tab <- matrix(c(c(rev(depth), depth[-1]), qu), ncol = 2,
                dimnames = list(nameLV(k)[[2]], c("depth","LV")))

  # confidence limits
  conf <- confintLV(x, k, alpha = alpha)

  cbind(tab, conf)
}


# confidence interval for k letter values
confintLV <- function(x, k, alpha=0.95) {
  n <- length(x)
  y <- sort(x)

  depth <- getDepth(k,n)
  extend <- ceiling(0.5 *sqrt(2*depth-1) * qnorm(alpha+(1-alpha)/2))
  low <- depth - extend
  high <- depth + extend
  clow <- pmax(1,ceiling(low))
  flow <- pmax(1,floor(low))
  chigh <- pmin(n, ceiling(high))
  fhigh <- pmin(n, floor(high))

  lvllow <- rev(rowMeans(cbind(y[clow],y[flow]), na.rm=T))
  if (length(lvllow) == 0) lvllow <- NA
  lvlhigh <- rev(rowMeans(cbind(y[chigh],y[fhigh]), na.rm=T))
  if (length(lvlhigh) == 0) lvlhigh <- NA
  # no 1 is the median - that's the last element in lvl
  lvulow <- rowMeans(cbind(y[n-chigh],y[n-fhigh]), na.rm=T)[-1]
  lvuhigh <- rowMeans(cbind(y[n-clow],y[n-flow]), na.rm=T)[-1]

  conf <- cbind(c(lvllow, lvulow), c(lvlhigh, lvuhigh))

  colnames(conf) <- c(paste((1 - alpha) / 2 * 100, "%" ,sep = "") ,
                      paste((alpha + (1 - alpha) / 2) * 100, "%", sep = ""))
  conf
}

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}


#' @rdname geom_lvplot
#' @param conf confidence level
#' @param percent numeric value: percent of data in outliers
#' @param k number of letter values shown
#' @param na.rm If \code{FALSE} (the default), removes missing values with
#'    a warning.  If \code{TRUE} silently removes missing values.
#' @section Computed/reported variables:
#' \describe{
#'   \item{k}{Number of Letter Values used for the display}
#'   \item{LV}{Name of the Letter Value}
#'   \item{width}{width of the interquartile box}
#' }
#' @export
stat_lvplot <- function(mapping = NULL, data = NULL, geom = "lvplot",
  position = "dodge", na.rm = TRUE, conf = 0.95, percent = NULL, k = NULL, show.legend = NA,
  inherit.aes = TRUE, ...)
{
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = StatLvplot,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      conf = conf,
      k = k,
      percent = percent,
      na.rm = na.rm,
      ...
    )
  )
}

#' @export
StatLvplot <- ggplot2::ggproto("StatLvplot", ggplot2::Stat,
  required_aes = c("x", "y"),
  non_missing_aes = "weight",

  setup_params = function(data, params) {
    params$width <- params$width %||% resolution(data$x) * 0.75
    params
  },

  compute_group = function(data, scales, width = NULL, na.rm = FALSE, k = NULL, conf=0.95, percent=NULL) {
    n <- nrow(data)

    k <- determineDepth(n, k, alpha=conf, percent)
    # compute letter values and outliers
    stats <- calcLV(data$y, k)
    outliers <- data$y < min(stats) | data$y > max(stats)
    res <- outputLVplot(data$y, stats, k, outliers, alpha=conf)

    df <- data.frame(res$letter.val)
    df$k <- k
    df$LV <- factor(row.names(df), levels=c("M", LETTERS[6:1], LETTERS[c(26:14, 12:7)]))
    df$ci <- list(res$conf.int)
    df$outliers <- list(data$y[res$outliers])

    df$x <- if (is.factor(data$x)) data$x[1] else mean(range(data$x))
    df$width <- width
    df$relvarwidth <- sqrt(n)
    yrange <- range(data$y)
    if (any(outliers)) {
      yrange <- range(c(res$letter.val[k,-1], data$y[!outliers]), na.rm = TRUE)
    }
    df$ymin <- yrange[1]
    df$ymax <- yrange[2]

    df
  }
)
