# Original VSUP code used with kind permission of Claus Wilke.
# See https://github.com/clauswilke/multiscales for more details

# -------------------------------------------------------------------------


zip <- function(...) purrr::transpose(list(...))

"%||%" <- function(a, b) if (!is.null(a)) a else b

# -------------------------------------------------------------------------




#' Train range for bivariate scale
#'
#' @importFrom scales train_continuous
#' @importFrom purrr transpose
#' @param new New data on which to train.
#' @param existing Existing range
#' @return A tibble containing two columns, `range1` and `range2`, each representing the trained
#' continuous range based on the new and existing data. This function is used to update or define
#' the scales of a bivariate analysis by considering both new input data and any existing range
#' specifications.
#' @export

train_bivariate <- function(new, existing = NULL) {

  if (is.null(new)) return(existing)
  range1 <- train_continuous(unlist(purrr::transpose(new)[[1]]), existing$range1)
  range2 <- train_continuous(unlist(purrr::transpose(new)[[2]]), existing$range2)
  tibble(range1, range2)
}




Range <- ggproto("Range", NULL,
                 range = NULL,
                 reset = function(self) {
                   self$range <- NULL
                 }
)

#' @rdname bivariate_range
#' @usage NULL
#' @export
RangeBivariate <- ggproto("RangeBivariate", Range,
                          train = function(self, x) {
                            self$range <- train_bivariate(x, self$range)
                          }
)


#' Constructor for bivariate range object
#' @export
bivariate_range <- function() {
  ggproto(NULL, RangeBivariate)
}

# -------------------------------------------------------------------------


# -------------------------------------------------------------------------


#' @rdname bivariate_scale
#' @usage NULL
#'
#' @importFrom scales rescale
#' @importFrom scales censor
#' @importFrom scales identity_trans
#' @importFrom scales zero_range
#' @importFrom dplyr tibble
#' @importFrom purrr transpose
#' @export
ScaleBivariate <- ggproto("ScaleBivariate",
                          Scale,
                          range = bivariate_range(),
                          rescaler = list(scales::rescale, scales::rescale),
                          oob = scales::censor,
                          #trans = list(identity_trans, dentity_trans),

                          is_discrete = function() FALSE,
                          is_bivariate = function() TRUE,

                          train = function(self, x) {
                            if (length(x) == 0) return()
                            self$range$train(x)
                          },

                          transform = function(self, x) {
                            ## fix for data frames
                            if (!is.list(x)) {
                              stop("For bivariate scale, aesthetic needs to be a list of two data columns. Did you forget `zip()`?", call. = FALSE)
                            }
                            x1 <- unlist(purrr::transpose(x)[[1]])
                            x2 <- unlist(purrr::transpose(x)[[2]])

                            x1 <- self$trans[[1]]$transform(x1)
                            x2 <- self$trans[[2]]$transform(x2)

                            ## fix for data frames
                            zip(x1, x2)
                          },

                          map = function(self, x, limits = self$get_limits()) {
                            ## fix for data frames
                            x1 <- unlist(purrr::transpose(x)[[1]])
                            x2 <- unlist(purrr::transpose(x)[[2]])

                            x1 <- self$rescaler[[1]](self$oob(x1, range = limits[[1]]), from = limits[[1]])
                            x2 <- self$rescaler[[2]](self$oob(x2, range = limits[[2]]), from = limits[[2]])

                            scaled <- self$palette(x1, x2)

                            ifelse(!is.na(scaled), scaled, self$na.value)
                          },

                          #  if scale contains a NULL, use the default scale range
                          #  if scale contains a NA, use the default range for that axis, otherwise
                          #  use the user defined limit for that axis
                          get_limits = function(self) {
                            if (self$is_empty()) return(dplyr::tibble(limits1 = c(0, 1), limits2 = c(0, 1)))

                            if (is.null(self$limits)) {
                              return(dplyr::tibble(limits1 = self$range$range[[1]], limits2 = self$range$range[[2]]))
                            } else {
                              limits1 <- ifelse(!is.na(self$limits[[1]]), self$limits[[1]], self$range$range[[1]])
                              limits2 <- ifelse(!is.na(self$limits[[2]]), self$limits[[2]], self$range$range[[2]])
                              return(dplyr::tibble(limits1, limits2))
                            }
                          },

                          get_breaks = function(self, limits = self$get_limits()) {
                            breaks1 <- self$get_breaks_1d(1, limits[[1]])
                            breaks2 <- self$get_breaks_1d(2, limits[[2]])

                            list(breaks1 = breaks1, breaks2 = breaks2)
                          },

                          # breaks for one data dimension
                          get_breaks_1d = function(self, i = 1, limits = self$get_limits()[[i]]) {
                            if (self$is_empty()) return(numeric(0))

                            # Limits in transformed space need to be converted back to data space
                            limits <- self$trans[[i]]$inverse(limits)

                            if (is.null(self$breaks)) {
                              return(NULL)
                            } else if (identical(self$breaks[[i]], NA)) {
                              stop("Invalid breaks specification. Use NULL, not NA")
                            } else if (scales::zero_range(as.numeric(limits))) {
                              breaks <- limits[[i]][1]
                            } else if (is.waive(self$breaks[[i]])) {
                              breaks <- self$trans[[i]]$breaks(limits)
                            } else if (is.function(self$breaks[[i]])) {
                              breaks <- self$breaks[[i]](limits)
                            } else {
                              breaks <- self$breaks[[i]]
                            }

                            breaks <- scales::censor(self$trans[[i]]$transform(breaks), self$trans[[i]]$transform(limits),
                                             only.finite = FALSE)
                            breaks
                          },

                          get_labels = function(self, breaks = self$get_breaks()) {
                            labels1 <- self$get_labels_1d(1, breaks[[1]])
                            labels2 <- self$get_labels_1d(2, breaks[[2]])

                            list(labels1 = labels1, labels2 = labels2)
                          },

                          # labels for one data dimension
                          get_labels_1d = function(self, i = 1, breaks = self$get_breaks()[[i]]) {
                            if (is.null(breaks)) return(NULL)

                            breaks <- self$trans[[i]]$inverse(breaks)

                            if (is.null(self$labels[[i]])) {
                              return(NULL)
                            } else if (identical(self$labels[[i]], NA)) {
                              stop("Invalid labels specification. Use NULL, not NA", call. = FALSE)
                            } else if (is.waive(self$labels[[i]])) {
                              labels <- self$trans[[i]]$format(breaks)
                            } else if (is.function(self$labels[[i]])) {
                              labels <- self$labels[[i]](breaks)
                            } else {
                              labels <- self$labels[[i]]
                            }
                            if (length(labels) != length(breaks)) {
                              stop("Breaks and labels are different lengths")
                            }
                            labels
                          },


                          clone = function(self) {
                            new <- ggproto(NULL, self)
                            new$range <- bivariate_range()
                            new
                          }
)


#' Constructor for bivariate scale object
#'
#' @inheritParams ggplot2::continuous_scale
#' @param limits Data frame with two columns of length two each defining the limits for the two data dimensions.
#' @param trans Either one transformation applied to both data dimensions or list of two transformations, one
#'   for each data dimension. Transformations can be given as either the name of a transformation object
#'   or the object itself. See [`ggplot2::continuous_scale()`] for details.
#' @param rescaler Either one rescaling function applied to both data dimensions or list of two rescaling functions,
#'   one for each data dimension.
#'
#' @importFrom scales as.trans
#' @export
bivariate_scale <- function(aesthetics, palette, name = waiver(),
                            breaks = waiver(), labels = waiver(), limits = NULL,
                            rescaler = scales::rescale, oob = scales::censor, expand = waiver(), na.value = NA_real_,
                            trans = "identity", guide = "none", super = ScaleBivariate,
                            scale_name = "bivariate_scale") {

  breaks <- bivariatize_arg(breaks, "breaks")
  labels <- bivariatize_arg(labels, "labels")

  trans <- bivariatize_arg(trans, "trans")
  trans[[1]] <- scales::as.trans(trans[[1]])
  trans[[2]] <- scales::as.trans(trans[[2]])

  rescaler <- bivariatize_arg(rescaler, "rescaler")

  if (!is.null(limits)) {
    # Check that limits are data frame or list with two columns of two values
    if (!is.list(limits)) {
      stop("Limits argument has to be a data frame or list of vectors", call. = FALSE)
    } else if (length(limits) != 2 || length(limits[[1]]) != 2 || length(limits[[2]]) != 2) {
      stop("Limits need to be two values each for both data dimensions", call. = FALSE)
    }

    # limits are given and valid, need to transform
    limits <- tibble(
      limits1 = trans[[1]]$transform(limits[[1]]),
      limits2 = trans[[2]]$transform(limits[[2]])
    )
  }

  ggproto(
    NULL, super,
    call = match.call(),

    aesthetics = aesthetics,
    scale_name = scale_name,
    palette = palette,

    range = bivariate_range(),
    limits = limits,
    trans = trans,
    na.value = na.value,
    expand = expand,
    rescaler = rescaler,
    oob = oob,

    name = name,
    breaks = breaks,

    labels = labels,
    guide = guide
  )
}

bivariatize_arg <- function(arg, name = "argument") {
  if (!is.null(oldClass(arg)) || is.function(arg) || is.atomic(arg)) {
    return(list(arg, arg))
  }

  if (!is.list(arg) || length(arg) != 2) {
    stop(paste0("In `bivariate_scale()`, argument `", name, "` needs to be given either as one argument applied to both data dimensions or as a list of exactly two arguments."), call. = FALSE)
  }

  arg
}



# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

#' Variance suppressing uncertainty palette
#'
#' Returns a palette function that turns `v` (value) and `u` (uncertainty) (both between 0 and 1) into
#' colors.
#' @param values Color values to be used at minimum uncertainty. Needs to be a vector of
#'   length `2^unc_levels`.
#' @param unc_levels Number of discrete uncertainty levels. The number of discrete colors
#'   at each level doubles.
#' @param max_light Maximum amount of lightening
#' @param max_desat Maximum amount of desaturation
#' @param pow_light Power exponent of lightening
#' @param pow_desat Power exponent of desaturation
#'
#' @importFrom scales colour_ramp
#' @importFrom colorspace desaturate
#' @importFrom colorspace lighten
#' @return A function that takes two parameters, `v` (value) and `u` (uncertainty), both expected to be
#' in the range of 0 to 1, and returns a color. This color is determined by the specified `values` colors
#' at minimum uncertainty, and modified according to the given `v` and `u` parameters to represent
#' uncertainty by adjusting lightness and saturation. The resulting function is useful for creating
#' color palettes that can encode both value and uncertainty in visualizations.
#'
#' @export
pal_vsup <- function(values, unc_levels = 4, max_light = 0.9, max_desat = 0, pow_light = 0.8, pow_desat = 1) {
  n <- 2^(unc_levels - 1)
  if (length(values) != n) {
    stop(length(values), " colors are provided but ", n, " colors are needed for ", unc_levels, " uncertainty levels.", call. = FALSE)
  }

  ramp <- scales::colour_ramp(values)

  # v = value, 0: small, 1: large
  # u = uncertainty, 0: completely certain, 1: completely uncertain
  map_to_discrete <- function(v, u) {
    j <- 1 + floor((1 - u) * unc_levels)
    j <- ifelse(j >= unc_levels, unc_levels, j)

    val_levels <- 2^(j-1) # total number of value levels at that uncertainty
    i <- 1 + floor(v * val_levels)
    i <- ifelse( i >= val_levels, val_levels, i)

    list(i = i,
         j = j,
         v = ((i - 0.5)/val_levels - 0.5/n)*n/(n - 1),
         u = 1 - (j - 1)/(unc_levels - 1))
  }

  function(v, u){
    x <- map_to_discrete(v, u)
    v <- x$v
    u <- x$u # need maximum lightening for 0 certainty

    range01 <- function(x,a,b){
       ((b-a)*((x-min(x)))/(max(x)-min(x)))+a
     }

    # limit maximal desaturation and lightening
    des_amt <- max_desat*u^pow_desat
    light_amt <- max_light*u^pow_light
    cols_des <- colorspace::desaturate(ramp(v), des_amt)
    #cols_des <- colorspace::lighten(ramp(v), range01(des_amt, 0, 0.9))
    nas <- is.na(light_amt)
    light_amt[nas] <- 0
    ifelse(nas, NA, colorspace::lighten(cols_des, light_amt, space = "HLS"))
  }
}


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

#' Colourfan guide
#'
#' @param title Title
#' @param title.x.position Title x position
#' @param title.y.position Title y position
#' @param title.theme Title theme
#' @param title.hjust Title hjust
#' @param title.vjust Title vjust
#' @param label Label
#' @param label.theme Label theme
#' @param barwidth Barwidth
#' @param barheight Barheight
#' @param nbin Number of bins
#' @param reverse Reverse
#' @param order order
#' @param available_aes Available aesthetics
#' @param ... Extra paramters
#'
#' @export
guide_colourfan <- function(

  # title
  title = waiver(),
  title.x.position = "top",
  title.y.position = "right",
  title.theme = NULL,
  title.hjust = 0.5,
  title.vjust = NULL, ## can be deleted?

  # label
  label = TRUE,
  label.theme = NULL,

  # bar
  barwidth = NULL,
  barheight = NULL,
  nbin = 32,

  # general
  reverse = FALSE,
  order = 0,
  available_aes = c("colour", "color", "fill"),

  ...) {

  if (!is.null(barwidth) && !grid::is.unit(barwidth)) barwidth <- unit(barwidth, default.unit)
  if (!is.null(barheight) && !grid::is.unit(barheight)) barheight <- unit(barheight, default.unit)

  structure(list(
    # title
    title = title,
    title.x.position = title.x.position,
    title.y.position = title.y.position,
    title.theme = title.theme,
    title.hjust = title.hjust,
    title.vjust = title.vjust,

    # label
    label = label,
    label.theme = label.theme,

    # bar
    barwidth = barwidth,
    barheight = barheight,
    nbin = nbin,

    # general
    reverse = reverse,
    order = order,

    # parameter
    available_aes = available_aes,
    ...,
    name = "colourfan"),
    class = c("guide", "colourfan")
  )
}


#' @export
guide_train.colourfan <- function(guide, scale, aesthetic = NULL) {

  # do nothing if scale are inappropriate
  if (length(intersect(scale$aesthetics, guide$available_aes)) == 0) {
    warning("colorfan guide needs appropriate scales: ",
            paste(guide$available_aes, collapse = ", "))
    return(NULL)
  }
  if (!scale$is_bivariate()) {
    warning("colorfan guide needs bivariate scales.")
    return(NULL)
  }

  # create tick positions and labels
  breaks <- scale$get_breaks()
  if (length(breaks[[1]]) == 0 && length(breaks[[2]]) == 0 ||
      all(is.na(breaks[[1]])) && all(is.na(breaks[[2]])))
    return()
  labels <- scale$get_labels(breaks)

  guide$ticks1 <- dplyr::tibble(value = breaks[[1]], label = labels[[1]])
  guide$ticks2 <- dplyr::tibble(value = breaks[[2]], label = labels[[2]])

  # needed to make guide show, even if this is not how we keep track of labels and ticks
  key <- as.data.frame(
    setNames(list(NA), aesthetic %||% scale$aesthetics[1]),
    stringsAsFactors = FALSE
  )
  guide$key <- key

  # fan specification
  limits <- scale$get_limits()
  v1 <- seq(limits[[1]][1], limits[[1]][2], length = guide$nbin)
  if (length(v1) == 0) {
    v1 = unique(limits[[1]])
  }
  v2 <- seq(limits[[2]][1], limits[[2]][2], length = guide$nbin)
  if (length(v2) == 0) {
    v2 = unique(limits[[2]])
  }
  # fan data matrix
  guide$fan <- expand.grid(x = v1, y = v2)
  guide$fan$colour <- scale$map(zip(guide$fan$x, guide$fan$y))

  # keep track of individual values along x and y also
  guide$fan.x <- v1
  guide$fan.y <- v2

  #guide$hash <- with(guide, digest::digest(list(title, ticks1, ticks2, name)))
  guide$hash <- with(guide,  rlang::hash(list(title, ticks1, ticks2, name)))
  guide
}





# simply discards the new guide
#' @export
guide_merge.colourfan <- function(guide, new_guide) {
  guide
}

# this guide is not geom-based.

#' @importFrom gtable gtable
#' @importFrom gtable gtable_add_grob
#'
#'
#' @export
guide_geom.colourfan <- function(guide, layers, default_mapping) {
  # Layers that use this guide
  # guide_layers <- plyr::llply(layers, function(layer) {
  #   matched <- matched_aes(layer, guide, default_mapping)
  #
  #   if (length(matched) && ((is.na(layer$show.legend) || layer$show.legend))) {
  #     layer
  #   } else {
  #     # This layer does not use this guide
  #     NULL
  #   }
  # })

  guide_layers <- lapply(layers, function(layer) {
    matched <- matched_aes(layer, guide, default_mapping)

    if (length(matched) && ((is.na(layer$show.legend) || layer$show.legend))) {
      layer
    } else {
      # This layer does not use this guide
      NULL
    }
  })

  # Remove this guide if no layer uses it
  #if (length(plyr::compact(guide_layers)) == 0) guide <- NULL
  if (length(purrr::compact(guide_layers)) == 0) guide <- NULL

  guide
}

#' @export
guide_gengrob.colourfan <- function(guide, theme) {
  title.x.position <- guide$title.x.position %||% "top"
  title.y.position <- guide$title.y.position %||% "right"

  fanwidth <- width_cm(theme$legend.key.width * 5)
  fanheight <- height_cm(theme$legend.key.height * 5)
  nbreak <- nrow(guide$key)

  # make the fan grob (`grob.fan`)
  grob.fan <- colourfan_grob(guide$fan$colour, nrow = guide$nbin, ncol = guide$nbin)

  # make ticks and labels
  # tick.x.pos <- rescale(
  #   guide$ticks1$value,
  #   c(0.5, guide$nbin - 0.5),
  #   guide$fan.x[c(1, length(guide$fan.x))]
  # ) / guide$nbin

  # tick.y.pos <- rescale(
  #   guide$ticks2$value,
  #   c(guide$nbin - 0.5, 0.5),
  #   guide$fan.y[c(1, length(guide$fan.y))]
  #   #guide$fan.y[c(15, length(guide$fan.y))]
  #   #guide$fan.y[c(26, length(guide$fan.y))]
  # ) / (guide$nbin)

  # this is where to change the legend tick positions
  tick.x.pos <- seq(0,1, length.out = 5)
  a <- c(0, .25, .5, .75, 1)
  #a<- rev(a)
  #a <- a + 0.125
  tick.y.pos <- a


  label.x.pos <- transform_radial(dplyr::tibble(x = tick.x.pos, y = 1), yoff = 0.04)
  label.y.pos <- transform_radial(dplyr::tibble(x = 1, y = tick.y.pos),
                                 # yoff = 0,
                                  xoff = 0.04)

  # get the label theme
  label.theme <- guide$label.theme %||% calc_element("legend.text", theme)

  # We break inheritance for hjust and vjust, because that's more intuitive here; it still allows manual
  # setting of hjust and vjust if desired. The alternative is to ignore hjust and vjust altogether, which
  # seems worse
  if (is.null(guide$label.theme$hjust) && is.null(theme$legend.text$hjust)) label.theme$hjust <- NULL
  if (is.null(guide$label.theme$vjust) && is.null(theme$legend.text$vjust)) label.theme$vjust <- NULL

  # label.theme in param of guide_legend() > theme$legend.text.align > default
  hjust <- label.theme$hjust %||% 0.5
  vjust <- label.theme$vjust %||% 0.5

  if (!guide$label) # are we drawing labels?
    grob.label.x <- NULL
  else {
    x <- unit(fanwidth*label.x.pos$x, "cm")
    y <- unit(fanheight*label.x.pos$y, "cm")
    margin_x <- FALSE
    margin_y <- FALSE

    label <- guide$ticks1$label

    # If any of the labels are quoted language objects, convert them
    # to expressions. Labels from formatter functions can return these
    if (any(vapply(label, is.call, logical(1)))) {
      label <- lapply(
        label,
        function(l) {
          if (is.call(l)) substitute(expression(x), list(x = l))
          else l
        }
      )
      label <- do.call(c, label)
    }
    grob.label.x <- element_grob(
      element = label.theme,
      label = label,
      x = x,
      y = y,
      hjust = 0.5,
      vjust = 0,
      margin_x = margin_x,
      margin_y = margin_y
    )
    grob.label.x <- ggname("guide.label.x", grob.label.x)
  }

  label.x.width <- width_cm(grob.label.x)
  label.x.height <- height_cm(grob.label.x)

  if (!guide$label) # are we drawing labels?
    grob.label.y <- NULL
  else {
    x <- unit(fanwidth*label.y.pos$x, "cm")
    y <- unit(fanheight*label.y.pos$y, "cm")
    margin_x <- FALSE
    margin_y <- FALSE

    label <- guide$ticks2$label

    # If any of the labels are quoted language objects, convert them
    # to expressions. Labels from formatter functions can return these
    if (any(vapply(label, is.call, logical(1)))) {
      label <- lapply(
        label,
        function(l) {
          if (is.call(l)) substitute(expression(x), list(x = l))
          else l
        }
      )
      label <- do.call(c, label)
    }
    grob.label.y <- element_grob(
      element = label.theme,
      label = label,
      x = x,
      y = y,
      hjust = 0,
      vjust = 0.5,
      margin_x = margin_x,
      margin_y = margin_y
    )
    grob.label.y <- ggname("guide.label.y", grob.label.y)
  }

  label.y.width <- width_cm(grob.label.y)
  label.y.height <- height_cm(grob.label.y)

  # make titles

  # obtain the theme for the legend title. We need this both for the title grob
  # and to obtain the title fontsize.
  title.theme <- guide$title.theme %||% calc_element("legend.title", theme)

  title.hjust <- guide$title.hjust %||% theme$legend.title.align %||% title.theme$hjust %||% 0
  title.vjust <- guide$title.vjust %||% title.theme$vjust %||% 0.5

  # make title grobs if needed
  title.x.label <- guide$title[1]
  if (is.null(title.x.label) || is.na(title.x.label)) {
    title.x.position <- "none"
  } else {
    grob.title.x <- ggname(
      "guide.title.x",
      element_grob(
        title.theme,
        label = title.x.label,
        hjust = title.hjust,
        vjust = title.vjust,
        margin_x = TRUE,
        margin_y = TRUE
      )
    )
    title.x.width <- width_cm(grob.title.x)
    title.x.height <- height_cm(grob.title.x)
  }

  title.y.label <- guide$title[2]
  if (is.null(title.y.label) || is.na(title.y.label)) {
    title.y.position <- "none"
  } else {
    title.y.pos <- transform_radial(
      tibble(x = 1, y = 0.5), xoff = 0.55
    )

    grob.title.y <- element_grob(
      element = title.theme,
      label = title.y.label,
      x = unit(fanwidth * title.y.pos$x, "cm"),
      y = unit(fanheight * title.y.pos$y, "cm"),
      hjust = 0.4,
      vjust = 0,
      angle = 60,
      margin_x = FALSE,
      margin_y = FALSE
    )
    title.y.width <- width_cm(grob.title.y)
    title.y.height <- height_cm(grob.title.y)
  }

  # gap between keys etc
  # the default horizontal and vertical gap need to be the same to avoid strange
  # effects for certain guide layouts
  title_fontsize <- title.theme$size %||% calc_element("legend.title", theme)$size %||% 0
  hgap <- width_cm(theme$legend.spacing.x  %||% (0.25 * unit(title_fontsize, "pt")))
  vgap <- height_cm(theme$legend.spacing.y %||% (0.25 * unit(title_fontsize, "pt")))

  # legend padding
  padding <- grid::convertUnit(theme$legend.margin %||% margin(), "cm")

  # we set up the entire legend as an 11x11 table which contains:
  # margin, title, gap, labels, ticks, fan, ticks, labels, gap, title, margin
  # depending on where titles and labels are added, some cells remain empty

  widths <- c(padding[4], 0, 0, 0, 0, fanwidth, 0, 0, 0, 0, padding[2])
  heights <- c(padding[1], 0, 0, 0, 0, fanheight, 0, 0, 0, 0, padding[3])


  heights[4] <- label.x.height - fanheight*(1 - min(label.x.pos$y))
  widths[8] <- label.y.width - fanwidth*(1 - min(label.y.pos$x))

  # titles
  grob.title.x.top <- NULL
  grob.title.x.bottom <- NULL
  if (title.x.position %in% c("top", "both")) {
    heights[2] <- title.x.height
    heights[3] <- vgap
    grob.title.x.top <- justify_grobs(
      grob.title.x,
      hjust = title.hjust,
      vjust = title.vjust,
      int_angle = title.theme$angle,
      debug = title.theme$debug
    )
  }
  if (title.x.position %in% c("bottom", "both")) {
    heights[10] <- title.x.height
    heights[9] <- vgap
    grob.title.x.bottom <- justify_grobs(
      grob.title.x,
      hjust = title.hjust,
      vjust = title.vjust,
      int_angle = title.theme$angle,
      debug = title.theme$debug
    )
  }

  grob.title.y.left <- NULL
  grob.title.y.right <- NULL
  if (title.y.position %in% c("right", "both")) {
    grob.title.y.right <- grob.title.y
  }

  # background
  grob.background <- element_render(theme, "legend.background")

  gt <- gtable(widths = unit(widths, "cm"), heights = unit(heights, "cm"))
  gt <- gtable_add_grob(
    gt, grob.background, name = "background", clip = "off",
    t = 1, r = -1, b = -1, l = 1
  )
  gt <- gtable_add_grob(
    gt, grob.fan, name = "fan", clip = "off",
    t = 6, r = 6, b = 6, l = 6
  )
  if (!is.null(grob.title.x.top)) {
    gt <- gtable_add_grob(
      gt, grob.title.x.top, name = "title.x.top", clip = "off",
      t = 2, r = 6, b = 2, l = 6
    )
  }
  if (!is.null(grob.label.x)) {
    gt <- gtable_add_grob(
      gt, grob.label.x, name = "label.x.top", clip = "off",
      t = 6, r = 6, b = 6, l = 6
    )
  }
  if (!is.null(grob.title.x.bottom)) {
    gt <- gtable_add_grob(
      gt, grob.title.x.bottom, name = "title.x.bottom", clip = "off",
      t = 10, r = 6, b = 10, l = 6
    )
  }
  if (!is.null(grob.title.y.left)) {
    gt <- gtable_add_grob(
      gt, grob.title.y.left, name = "title.y.left", clip = "off",
      t = 6, r = 2, b = 6, l = 2
    )
  }
  if (!is.null(grob.title.y.right)) {
    gt <- gtable_add_grob(
      gt, grob.title.y.right, name = "title.y.right", clip = "off",
      t = 6, r = 6, b = 6, l = 6
    )
  }
  if (!is.null(grob.label.y)) {
    gt <- gtable_add_grob(
      gt, grob.label.y, name = "label.y.top", clip = "off",
      t = 6, r = 6, b = 6, l = 6
    )
  }


  gt
}

#' @importFrom grid polygonGrob
#' @importFrom grid gpar
#' @importFrom grid is.grob
#' @importFrom grid convertWidth
#' @importFrom grid grobWidth
#' @importFrom grid is.unit
#' @importFrom grid convertHeight
#' @importFrom grid grobHeight
#'
#' @return  A `grob` object representing a color fan. This `grob` can be added to a grid-based plot
#' or a ggplot2 object to visualize a range of colors in a fan-like structure. Each segment of the fan
#' corresponds to a color specified in the `colours` parameter, allowing for an intuitive representation
#' of color gradients or palettes.

#' @export
#' @rdname guide_colourfan
guide_colorfan <- guide_colourfan


colourfan_grob <- function(colours, nrow, ncol, nmunch = 10) {
  # the trick is that we first make square polygons and then transform coordinates
  dx <- 1 / ncol
  dy <- 1 / nrow

  # grid of base points
  x <- rep((0:(ncol-1))/ncol, nrow)
  y <- rep(((nrow-1):0)/nrow, each = ncol)

  # turn into polygon boundaries
  x <- unlist(lapply(x, function(x) c(x+dx*(0:nmunch)/nmunch, x+dx*(nmunch:0)/nmunch)))
  y <- unlist(lapply(y, function(y) c(rep(y, nmunch + 1), rep(y+dy, nmunch + 1))))
  id <- rep(1:(nrow*ncol), each = 22)

  # now transform coordinates and make polygon
  data <- transform_radial(tibble(x, y))
  grid::polygonGrob(data$x, data$y, id, gp = grid::gpar(fill = colours, col = colours, lwd = 0.5, lty = 1))
}


# map square into fan
# assumes x and y run from 0 to 1
# x runs left to right
# y runs top to bottom
transform_radial <- function(data, xoff = 0, yoff = 0) {
  phi <- (data$x * 60 - 30)*(pi/180)
  Y <- (data$y + yoff) * cos(phi) - xoff * sin(60*pi/360)
  X <- (data$y + yoff) * sin(phi) + 0.5 + xoff * cos(60*pi/360)
  tibble(x = X, y = Y)
}


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------


width_cm <- function(x) {
  if (grid::is.grob(x)) {
    grid::convertWidth(grid::grobWidth(x), "cm", TRUE)
  } else if (grid::is.unit(x)) {
    grid::convertWidth(x, "cm", TRUE)
  } else if (is.list(x)) {
    vapply(x, width_cm, numeric(1))
  } else {
    stop("Unknown input")
  }
}

height_cm <- function(x) {
  if (grid::is.grob(x)) {
    grid::convertHeight(grid::grobHeight(x), "cm", TRUE)
  } else if (grid::is.unit(x)) {
    grid::convertHeight(x, "cm", TRUE)
  } else if (is.list(x)) {
    vapply(x, height_cm, numeric(1))
  } else {
    stop("Unknown input")
  }
}

matched_aes <- function(layer, guide, defaults) {
  all <- names(c(layer$mapping, if (layer$inherit.aes) defaults, layer$stat$default_aes))
  geom <- c(layer$geom$required_aes, names(layer$geom$default_aes))
  matched <- intersect(intersect(all, geom), names(guide$key))
  matched <- setdiff(matched, names(layer$geom_params))
  setdiff(matched, names(layer$aes_params))
}

# not copied for now
# element_render <- ggplot2:::element_render
# ggname <- ggplot2:::ggname
# justify_grobs <- ggplot2:::justify_grobs
# is.waive <- ggplot2:::is.waive


# ---- ggname --------------------------------------------------------------
# Vendored and renamed from ggplot2 (original internal function)
#
# Description:
# This function adds a name to a grob object using a prefix, useful for debugging and layout management.

ggname <- function (prefix, grob)
{
  grob$name <- grid::grobName(grob, prefix)
  grob
}


# ---- element_render ------------------------------------------------------
# Vendored and renamed from ggplot2 (original internal function)
#
# Description:
# Renders a theme element (like axis text or legend title) into a grob. If the element is missing,
# returns a zeroGrob and informs the user.

element_render <- function (theme, element, ..., name = NULL)
{
  el <- ggplot2::calc_element(element, theme)
  if (is.null(el)) {
    cli::cli_inform("Theme element {.var {element}} is missing")
    return(ggplot2::zeroGrob())
  }
  grob <- ggplot2::element_grob(el, ...)
  ggname(paste(element, name, sep = "."), grob)
}

# ---- justify_grobs -------------------------------------------------------
# Vendored and modified from ggplot2 (original internal function)
#
# Description:
# Justifies one or more grobs to a specific x/y position using horizontal/vertical justification.
# Includes internal helpers for angle-based adjustment and input validation.

stop_input_type_alter <- function (x, what, ..., allow_na = FALSE, allow_null = FALSE,
                                   show_value = TRUE, arg = caller_arg(x), call = caller_env())
{
  cli <- rlang::env_get_list(nms = c("format_arg", "format_code"),
                             last = topenv(), default = function(x) sprintf("`%s`",
                                                                            x), inherit = TRUE)
  if (allow_na) {
    what <- c(what, cli$format_code("NA"))
  }
  if (allow_null) {
    what <- c(what, cli$format_code("NULL"))
  }
  if (length(what)) {
    what <- oxford_comma(what)
  }
  message <- sprintf("%s must be %s, not %s.", cli$format_arg(arg),
                     what, obj_type_friendly(x, value = show_value))
  abort(message, ..., call = call, arg = arg)
}



rotate_just_alter <-
  function (angle, hjust, vjust)
  {
    angle <- (angle %||% 0)%%360
    if (is.character(hjust)) {
      hjust <- match(hjust, c("left", "right")) - 1
      hjust[is.na(hjust)] <- 0.5
    }
    if (is.character(vjust)) {
      vjust <- match(vjust, c("bottom", "top")) - 1
      vjust[is.na(vjust)] <- 0.5
    }
    size <- vctrs::vec_size_common(angle, hjust, vjust)
    angle <-  vctrs::vec_recycle(angle, size)
    hjust <-  vctrs::vec_recycle(hjust, size)
    vjust <-  vctrs::vec_recycle(vjust, size)
    case <- findInterval(angle, c(0, 90, 180, 270, 360))
    hnew <- hjust
    vnew <- vjust
    is_case <- which(case == 2)
    hnew[is_case] <- 1 - vjust[is_case]
    vnew[is_case] <- hjust[is_case]
    is_case <- which(case == 3)
    hnew[is_case] <- 1 - hjust[is_case]
    vnew[is_case] <- 1 - vjust[is_case]
    is_case <- which(case == 4)
    hnew[is_case] <- vjust[is_case]
    vnew[is_case] <- 1 - hjust[is_case]
    list(hjust = hnew, vjust = vnew)
  }



justify_grobs <- function (grobs, x = NULL, y = NULL, hjust = 0.5, vjust = 0.5,
                                 int_angle = 0, debug = FALSE)
{
  if (!inherits(grobs, "grob")) {
    if (is.list(grobs)) {
      return(lapply(grobs, justify_grobs, x, y, hjust,
                    vjust, int_angle, debug))
    }
    else {
      stop_input_type_alter(grobs,
                            as_cli("an individual {.cls grob} or list of {.cls grob} objects"))
    }
  }
  if (inherits(grobs, "zeroGrob")) {
    return(grobs)
  }
  just <- rotate_just_alter(int_angle, hjust, vjust)
  x <- x %||% unit(just$hjust, "npc")
  y <- y %||% unit(just$vjust, "npc")
  if (isTRUE(debug)) {
    children <- grid::gList(grid::rectGrob(gp = grid::gpar(fill = "lightcyan",
                                                           col = NA)), grobs)
  }
  else {
    children = grid::gList(grobs)
  }
  result_grob <-  grid::gTree(children = children, vp = grid::viewport(x = x,
                                                                 y = y,
                                                                 width = grid::grobWidth(grobs),
                                                                 height = grid::grobHeight(grobs),
                                                                 just = unlist(just)))
  if (isTRUE(debug)) {
    grid::grobTree(result_grob, grid::pointsGrob(x, y, pch = 20, gp = grid::gpar(col = "mediumturquoise")))
  }
  else {
    result_grob
  }
}


# ---- is.waive -----------------------------------------------------------
# Vendored from ggplot2 (trivial internal function)
#
# Description:
# Checks whether an object is a ggplot2 'waiver' object, used to represent default behaviour.

is.waive <- function(x) {
  inherits(x, "waiver")
}


