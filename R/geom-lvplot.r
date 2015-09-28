#' Side-by-side LV boxplots.
#'
#' An extension of standard boxplots which draws k letter statistics.
#' Conventional boxplots (Tukey 1977) are useful displays for conveying rough
#' information about the central 50\% of the data and the extent of the data.
#' For moderate-sized data sets (\eqn{n < 1000}), detailed estimates of tail
#' behavior beyond the quartiles may not be trustworthy, so the information
#' provided by boxplots is appropriately somewhat vague beyond the quartiles,
#' and the expected number of ``outliers'' and ``far-out'' values for a
#' Gaussian sample of size \eqn{n} is often less than 10 (Hoaglin, Iglewicz,
#' and Tukey 1986). Large data sets (\eqn{n \approx 10,000-100,000}) afford
#' more precise estimates of quantiles in the tails beyond the quartiles and
#' also can be expected to present a large number of ``outliers'' (about
#' \eqn{0.4 + 0.007 n}).
#' The letter-value box plot addresses both these shortcomings: it conveys
#' more detailed information in the tails using letter values, only out to the
#' depths where the letter values are reliable estimates of their
#' corresponding quantiles (corresponding to tail areas of roughly
#' \eqn{2^{-i}}); ``outliers'' are defined as a function of the most extreme
#' letter value shown. All aspects shown on the letter-value boxplot are
#' actual observations, thus remaining faithful to the principles that
#' governed Tukey's original boxplot.
#'
#' @seealso \code{\link{stat_quantile}} to view quantiles conditioned on a
#'   continuous variable.
#' @inheritParams ggplot2::geom_point
#' @param geom,stat Use to override the default connection between
#'   \code{geom_lvplot} and \code{stat_lvplot}.
#' @param outlier.colour Override aesthetics used for the outliers. Defaults
#'   come from \code{geom_point()}.
#' @param outlier.shape Override aesthetics used for the outliers. Defaults
#'   come from \code{geom_point()}.
#' @param outlier.size Override aesthetics used for the outliers. Defaults
#'   come from \code{geom_point()}.
#' @param outlier.stroke Override aesthetics used for the outliers. Defaults
#'   come from \code{geom_point()}.
#' @param varwidth if \code{FALSE} (default) draw boxes that are the same size for each group. If
#'   \code{TRUE}, boxes are drawn with widths proportional to the
#'   square-roots of the number of observations in the groups (possibly
#'   weighted, using the \code{weight} aesthetic).
#' @export
#' @references McGill, R., Tukey, J. W. and Larsen, W. A. (1978) Variations of
#'     box plots. The American Statistician 32, 12-16.
#' @examples
#' library(ggplot2)
#' p <- ggplot(mpg, aes(class, hwy))
#' p + geom_lvplot(aes(fill=..LV..)) + scale_fill_brewer()
#' p + geom_lvplot() + geom_jitter(width = 0.2)
#' p + geom_lvplot(alpha=1, aes(fill=..LV..)) + scale_fill_brewer()
#'
#' # Outliers
#' p + geom_lvplot(varwidth = TRUE, aes(fill=..LV..)) + scale_fill_brewer()
#' p + geom_lvplot(fill = "white", colour = "#3366FF")
#' p + geom_lvplot(outlier.colour = "red", outlier.shape = 1)
#'
#' # Plots are automatically dodged when any aesthetic is a factor
#' p + geom_lvplot(aes(fill = drv))
#'
#' \dontrun{
#' # just for now: read On_Time data into object ot from lvplot paper
#' library(RColorBrewer)
#' cols <- c("white",brewer.pal(9, "Greys"), "Red", "Pink")
#' reds <- brewer.pal(5, "Reds")[-1]
#' cols <- c("grey50", cols[2:3], reds[1], cols[4:6], reds[2], cols[7:9], reds[3], cols[10:11])
#'
#' ggplot(data=ot) + geom_lvplot(aes(x=UniqueCarrier, y=sqrt(TaxiOut+TaxiIn),
#'   fill=..LV..), alpha=1) +
#'   scale_fill_manual(values=cols)
#' ggplot(data=ot) + geom_lvplot(aes(x=factor(DayOfWeek), y=sqrt(TaxiOut+TaxiIn),
#'   fill=..LV..), alpha=1) +
#'   scale_fill_manual(values=cols) + facet_wrap(~UniqueCarrier)
#' }
geom_lvplot <- function(mapping = NULL, data = NULL, stat = "lvplot",
  position = "dodge", outlier.colour = "black", outlier.shape = 19,
  outlier.size = 1.5, outlier.stroke = 0.5,
  varwidth = FALSE, show.legend = NA, inherit.aes = TRUE, ...)
{
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomLvplot,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      outlier.colour = outlier.colour,
      outlier.shape = outlier.shape,
      outlier.size = outlier.size,
      outlier.stroke = outlier.stroke,
      varwidth = varwidth,
      ...
    )
  )
}

#' @importFrom grid grobTree
GeomLvplot <- ggplot2::ggproto("GeomLvplot", ggplot2::Geom,
  setup_data = function(data, params) {
  #  browser()
    data$width <- data$width %||%
      params$width %||% (resolution(data$x, FALSE) * 0.9)

    if (!is.null(data$outliers)) {
      suppressWarnings({
        out_min <- vapply(data$outliers, min, numeric(1))
        out_max <- vapply(data$outliers, max, numeric(1))
      })

      data$ymin_final <- pmin(out_min, data$ymin)
      data$ymax_final <- pmax(out_max, data$ymax)
    }

    # if `varwidth` not requested or not available, don't use it
    if (is.null(params) || is.null(params$varwidth) || !params$varwidth || is.null(data$relvarwidth)) {
      data$xmin <- data$x - data$width / 2
      data$xmax <- data$x + data$width / 2
    } else {
      # make `relvarwidth` relative to the size of the largest group
      data$relvarwidth <- data$relvarwidth / max(data$relvarwidth)
      data$xmin <- data$x - data$relvarwidth * data$width / 2
      data$xmax <- data$x + data$relvarwidth * data$width / 2
    }
    data$width <- NULL
    if (!is.null(data$relvarwidth)) data$relvarwidth <- NULL

    data
  },

  draw_group = function(data, panel_scales, coord,
                        outlier.colour = "black", outlier.shape = 19,
                        outlier.size = 1.5, outlier.stroke = 0.5,
                        varwidth = FALSE) {
    common <- data.frame(
      colour = data$colour,
      size = data$size,
      linetype = data$linetype,
      fill = alpha(data$fill, data$alpha),
      group = data$group,
      stringsAsFactors = FALSE
    )

    i <- seq_len(data$k[1]-1)-1
    width <- data$xmax - data$xmin
    offset <- c(0, (i / (2 * data$k[1])))  * width
    box <- data.frame(
      xmin = data$xmin + offset,
      xmax = data$xmax - offset,
      ymin = data$lower,
      ymax = data$upper,
      alpha = data$alpha,
      LV = data$LV,
      common,
      stringsAsFactors = FALSE
    )
    box <- box[nrow(box):1,]

    medians <- subset(box, LV=="M")
    medians <-   transform(medians,
        colour = fill,
        x = xmin,
        xend = xmax,
        y = ymin,
        yend = ymax
      )

    if (!is.null(data$outliers) && length(data$outliers[[1]] >= 1)) {
      outliers <- data.frame(
        y = data$outliers[[1]],
        x = data$x[1],
        colour = outlier.colour %||% data$colour[1],
        shape = outlier.shape %||% data$shape[1],
        size = outlier.size %||% data$size[1],
        stroke = outlier.stroke %||% data$stroke[1],
        fill = NA,
        alpha = NA,
        stringsAsFactors = FALSE
      )
      outliers_grob <- GeomPoint$draw_panel(outliers, panel_scales, coord)
    } else {
      outliers_grob <- NULL
    }

    ggplot2:::ggname("geom_lvplot", grobTree(
      outliers_grob,
      GeomRect$draw_panel(box, panel_scales, coord),
      GeomSegment$draw_panel(medians, panel_scales, coord)
    ))
  },

  draw_key = ggplot2::draw_key_rect,

  default_aes = ggplot2::aes(weight = 1, colour = "grey70", fill = "grey70", size = 0.5,
    alpha = 1, shape = 19, linetype = "solid", outlier.colour = "black",
    outlier.shape = 19, outlier.size = 1.5, outlier.stroke = 0.5),

  required_aes = c("x", "k", "LV")
)
