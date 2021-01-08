# Setup environment ------------------------------------------------------------

## Import library
library(reshape2)
library(lme4)
library(ggplot2)
library(car)
library(readxl)
library(reshape2)

## Helper functions
`%=>%` = function(lhs, rhs) {
    rhs = substitute(rhs)
    if (is.symbol(rhs)) rhs = as.call(c(rhs, quote(..)))
    if (length(rhs) == 1) rhs = as.call(c(rhs[[1L]], quote(..)))
    eval(rhs, envir = list(.. = lhs), enclos = parent.frame())
}

`%>=>%` = function(lhs, rhs) {
    rhs = substitute(rhs)
    if (is.symbol(rhs)) rhs = as.call(c(rhs, quote(..)))
    if (length(rhs) == 1) rhs = as.call(c(rhs[[1L]], quote(..)))
    eval(rhs, envir = list(.. = lhs), enclos = parent.frame())
    lhs
}

`%->%` = function(lhs, rhs) {
    invisible(eval.parent(substitute((rhs = lhs))))
}

mutate = function(.data, ...) { # . is used to deal with partial match in ...
    .cond = vapply(substitute(...()), .x ->> deparse(.x, 500), NA_character_)
    names(.cond) = ifelse(names(.cond) == "", .cond, names(.cond))
    for (i in seq_along(.cond)) .data[, names(.cond)[i]] = eval(
        str2lang(.cond[i]), envir = .data, enclos = parent.frame()
    ) #Don't change to lappy for realtime update!
    return(.data)
}

## Import scripts
	"/crossvalidation.R" %=>% paste0(getwd(), ..) %=>% source
	"/goftest.R" %=>% paste0(getwd(), ..) %=>% source

# Import data ------------------------------------------------------------------

## Data - weather
"Data/together.xlsx" %=>%
	readxl::read_excel(.., sheet = "weather", na = "NA") %=>%
	setNames(.., c("trial", "wpi", "temp", "rh", "dew")) %->%
	datw1

## Data - field
readxl::read_excel("Data/together.xlsx", sheet = "inc", na = "NA") %=>%
	setNames(.., c("loc", "trial", "block", "dis", "ini", 1:4)) %=>%
	melt(..,
		id = c("loc", "trial", "block", "dis", "ini"),
		var = "wpi",
		value.name = "presence"
	) %=>%

# Model fitting ----------------------------------------------------------------

## Format data
	mutate(..,
		dpi = as.numeric(as.character(wpi)) * 7 - ifelse(loc == "C", 2, 0),
		loc = as.factor(loc),
		ini = factor(ini, labels = c("low", "high")),
		trial = as.factor(trial),
		block = as.factor(block)
	) %=>%

## Merge field and weather data
	merge(.., datw1, all.x = T) %=>%
	na.omit %->% df5 %=>%

## Fit a GLMM model
	glmer(
		presence ~ dis * dpi + ini + temp + dis:temp +(1 | loc / trial / block),
		data = ..,
		family = binomial(link = "logit")
	) %->% fit %>=>%

## Model summary
	print(summary(..)) %>=>%

## Significance of fixed effects
	print(Anova(..)) %>=>%

## Goodness-0f-fit
	gofTest(..) %=>%

## Model accuracy and cross validation
	xVal(.., leave = 10, fold = 10) %>=>%

## Likelihood ratio / Model selection
	print(anova(
		update(.., formula = presence ~ (1 | block)), # Null model
		# Main effect model
		update(.., formula = presence ~ dis + dpi + ini + temp + (1 | block)),
		..
	)) %=>%

## Save model
	save(.., fit.rda)

# Model analysis ---------------------------------------------------------------

## Fidelity of prediction pot
widths = function(v) {
  w = sort(unique(v))
  x = w[-length(w)] + diff(w) / 2
  x = c(0, x, 2 * x[length(x)] - x[length(x)-1])
  y = sapply(seq_along(w), function(a) 2 * abs(max(w[a] - x[a], x[a+1] - w[a])))
  z = setNames(y, w)
  return(unname(z[as.character(v)]))
}

data.frame(fit = fitted(fit), model.frame(fit)) %=>%
	with(.., aggregate(list(presence, fit), by = list(dis, dpi, loc), mean)) %=>%
	setNames(.., c("dis", "dpi", "loc", "presence", "fitted")) %=>%
	reshape2::melt(.., id = c("dis", "dpi", "loc"), var = "source") %=>%
	mutate(.., week = (dpi + 2) %/% 7) %=>% #because commercial is oberved early
	split(.., ..$loc) %=>%
	lapply(.., .z ->> mutate(.z, height = widths(dis))) %=>%
	do.call(rbind, ..) %=>%
	ggplot(.., aes(week, dis, height = height)) +
	geom_tile(aes(fill = value)) +
	facet_grid(
		loc~source, space = "free", scale = "free", labeller = labeller(
			source = c(presence = "Observed", fitted = "Predicted"),
			loc = c(C = "Commercial greenhouse", G = "GCREC greenhouse")
		)
	) +
	scale_x_continuous("Week post-inoculation",
		breaks = 1:4,
		expand = c(0, 0)
	) +
	scale_y_continuous("Distance from point of inoculation (cm)",
		breaks = x ->> if (x[2] > 200) c(3, 33 * 1:8) else
			c(8, 24, 41, 63, 84, 104, 131, 155, 179), # Distance collected
		expand = c(0, 0)
	) +
	scale_fill_gradient2("Presence of\n bacteria",
		low = "forestgreen", mid = "yellow", high = "red2", midpoint = 0.5
	) +
	coord_cartesian(ylim = c(0, NA)) +
	theme_classic() +
	theme(
		axis.line = element_blank(),
		panel.spacing.x = unit(4, "mm"),
		panel.spacing.y = unit(10, "mm"),
		strip.background = element_blank(),
		plot.title = element_text(hjust = 0.5),
		axis.title = element_text("OpenSans", size = 16),
		axis.text = element_text("OpenSans", size = 14),
		legend.title = element_text("OpenSans", face = "bold", size = 16),
		legend.text = element_text("OpenSans", size = 14),
		strip.text = element_text("OpenSans", face = "bold", size = 16)
	) %=>%
	ggsave("fidelity.png", .., type = "cairo", width = 8, height = 8)

## Prediction plot (over distance and time at various conditions)

expand.grid(
	dpi = 0:12,
	dis = seq(0, 200, 1),
	temp = c(22, 26),
	ini = c("low", "high")
) %=>%
	data.frame(fit = predict(fit, .., type = "resp", re.form = NA), ..) %=>%
	mutate(..,
		fitb = ifelse(fit < 0.1, 0, 1),
		symp = ifelse(ini == "low", 8, 5),
		temp = ifelse(temp < 24, "cold temperature", "warm temperature"),
		ini = ifelse(ini == "low", "low inoculum", "high inoculum")
	) %->% dfxf %=>%
	..[..$fit >= 0.05 & ..$dpi == ..$symp, ] %=>%
	split(.., list(..$temp, ..$ini)) %=>%
	lapply(.., .z ->> .z[which.max(.z$dis), ]) %=>%
	do.call(rbind, ..) %=>%
	mutate(..,
		lab = "Presence: 0.05", # Annotation text
		ang = c(36, 32, 20, 18), # Rotate annotation
		x = 10,
		y = c(85, 162, 40, 73) # Height of annotation
	) %->% annot

ggplot(dfxf, aes(x = dpi, y = dis, fill = fit)) +
	geom_raster(interpolate = T) +
	geom_contour(
		aes(z = fit, col = factor(..level.. == 0.05)),
		breaks = 0.05 * 0:20
	) +
	facet_grid( ini ~ temp) +
	scale_x_continuous("Days postinoculation",
		breaks = 0:14,
		expand = c(0, 0)
	) +
	scale_y_continuous("Distance (cm)",
		breaks = seq(0, 200, 50),
		expand = c(0, 0),
		sec.axis = dup_axis()
	) +
	scale_colour_manual(values = c("transparent", "red"), guide = F) +
	scale_fill_gradientn("Pathogen\nincidence",
		colors = c("white", "grey70", "grey50", "grey30", "grey10", "black"),
		limits = c(0, 1)
	) +
	geom_vline(aes(xintercept = symp), dummy, linetype = 2, col = "black") +
	geom_hline(aes(yintercept = dis), dummy, linetype = 2, col = "black") +
	geom_text(
		aes(x = x, y = y, label = lab, fill = NULL, angle = ang),
		data = dummy, col = "red", size = 6) +
	theme(
		strip.background = element_blank(),
		strip.text = element_text(family = "OpenSans", face = "bold", size = 24),
		strip.placement = "outside",
		strip.switch.pad.grid = unit(1, "lines"),
		panel.border = element_rect(fill = NA, color = "grey10"),
		panel.spacing = unit(2, "lines"),
		axis.title = element_text(family = "OpenSans", face = "bold", size = 22),
		axis.text = element_text(family = "OpenSans", size = 20),
		legend.title = element_text(family = "OpenSans", face = "bold", size = 24),
		legend.title.align = 0.5,
		legend.text = element_text(family = "OpenSans", size = 22),
		legend.key.height = unit(1, "in"),
	) %=>%
	ggsave(
		filename = "pred2.pdf", .., device = "cairo_pdf",
		width = 15, height = 12, dpi = 300
	)

## Effect of temperature plot

expand.grid(
	ini = c("low", "high"),
	dis = seq(0, 1000, 0.5),
	temp = seq(20, 27, 0.1)
) %=>%
	mutate(.., dpi = ifelse(ini == "low", 8, 5)) %=>%
	data.frame(fit = predict(fit, .., type = "resp", re.form = NA), ..) %=>%
	split(.., list(..$temp, ..$ini)) %=>%
	lapply(.., .z ->> .z[which.min(abs(0.05 - .z$fit)), ]) %=>%
	do.call(rbind, ..) %=>%
	ggplot(.., aes(temp, dis, col = ini)) +
		geom_smooth(size = 1.2, se = F) +
		ggtitle(expression(paste(
			italic("X. perforans"), " dispersal by latent period"
		))) +
		scale_x_continuous("Temperature (°C)",
			breaks = 20:27,
			expand = c(0, 0.2)
		) +
		scale_y_continuous("Distance from point of inoculation (cm)",
			breaks = seq(0, 160, 20),
			expand = c(0.01, 0),
			limits = c(0, 160)) +
		scale_color_manual("Initial inoculum load:",
			values = c("low" = "#0a9574", "high" = "#2772d7"),
			labels = c("low" = "Low", "high" = "High")
		) +
		theme_classic() +
		theme(
			axis.title = element_text(family = "OpenSans", size = 20),
			axis.text = element_text(family = "OpenSans", size = 18),
			legend.title = element_text(family = "OpenSans", size = 20),
			legend.title.align = 0.5,
			legend.text = element_text(family = "OpenSans", size = 18),
			legend.position = "bottom",
			plot.title = element_text(
				size = 26,
				hjust = 0.5,
				margin = margin(0, 0, 10, 0)
			)
		) %=>%
	ggsave(
		filename = "temperature.pdf", .., device = "cairo_pdf",
		width = 8, height = 8, dpi = 300
	)

# Effect of temperature --------------------------------------------------------

## Percentage increase

tl = 22; th = 23; dl = 5; dh = 10

dltl = paste(dl, tl, sep = ".")
dlth = paste(dl, th, sep = ".")
dhtl = paste(dh, tl, sep = ".")
dhth = paste(dh, th, sep = ".")

expand.grid(
    dpi = c(dl, dh),
    dis = seq(0, 200, 0.01),
    temp = c(tl, th),
    ini = c("high")
) %=>%
    data.frame(fit = predict(fit, .., type = "resp", re.form = NA), ..) %=>%
    mutate(.., dev = abs(fit - 0.05)) %=>%
    split(.., list(..$dpi, ..$temp)) %=>%
    sapply(.., .z ->> .z$dis[which.min(.z$dev)]) %=>%
    scales::label_percent()(
        ((..[dhth] - ..[dlth]) - (..[dhtl] - ..[dltl])) / (..[dhtl] - ..[dltl])
    )