# get-coef.R
# Helper function for making coefficient table ####
# Kevin Chen
# May 28, 2020

get.coef <- function(
	model.path = to_drive_D(here::here(
		"left truncation/resources",
		ifelse(length(probs) - 1 == 20, "mod.gam.rds", "small.mod.gam.rds"))),
	output.dir = here::here("reports/left truncation/survival to 1985/resources"),
	output.name = NULL) {

	if (is.null(output.name)) {
		output.name <- ifelse(length(probs) - 1 == 20, "gam.tab.rds", "small.gam.tab.rds")
		if (grepl("glm\\.rds", model.path)) {
			output.name <- gsub("gam\\.tab", "glm.tab", output.name)
		}
	}

	dir.create(output.dir, F, T)

	if (!grepl("/$", output.dir)) {
		output.dir <- paste0(output.dir, "/")
	}
	if (!grepl("^/", output.name)) {
		output.name <- gsub("^/", "", output.name)
	}

	# For case counts
	dat <- as.data.table(as.data.frame(cohort_analytic[year <= 1985]))

	dat <- dat[y == 1,.(
		`Years_since_hire`,
		`Calendar_year`,
		`Age`,
		Sex, Race, Plant,
		`Time_spent_in_assembly`,
		`Time_spent_machining`,
		`Time_spent_off`,
		`Cumulative_time_off`,
		`Year_of_hire`,
		`Cumulative_soluble_exposure`,
		`Cumulative_straight_exposure`,
		`Cumulative_synthetic_exposure`,
		`Employment_status`
	), by = .(studyno)]

	# Load model results
	tmp.mod <- readRDS(model.path)

	# Get coefficients
	if (grepl("glm\\.rds", model.path)) {
		tmp.tab <- summary(tmp.mod)$coefficients
	}
	if (grepl("gam\\.rds", model.path)) {
		library(mgcv)
		tmp.tab <- rbind(summary(tmp.mod)$p.table,
										 summary(tmp.mod)$s.table)
	}

	tmp.tab <- as.data.frame(cbind(
		rownames(tmp.tab), as.data.frame(tmp.tab)))

	setDT(tmp.tab)
	names(tmp.tab) <- c("Covariate", "Estimate", "SE", "z", "p")

	# splined covariates
	if (grepl("gam\\.rds", model.path)) {
		splines.tab <- rbindlist(list(
			# Years since hire (Spline)
			tmp.tab[grepl("sincehire\\.years", Covariate),],
			# # Calendar year (Spline)
			# tmp.tab[grepl("^s\\(year\\)$", Covariate),],
			# Age (Spline)
			tmp.tab[grepl("age.year2", Covariate),]

		))

		splines.tab$Covariate <- c("Splined time since hire",
																 # "Splined calendar year",
																 "Splined age")

		names(splines.tab) <- c("Covariate", "Estimated $df$", "Reference $df$", "$\\Chi^2$", "p")
		}

	# Cleaned up table
	tmp.tab <- tmp.tab[, .(
		Covariate = gsub("`|\\(|\\[", "", substr(Covariate, 1, sapply(gregexpr("`|\\(|\\[", Covariate), max))),
		level = gsub(",", ", ", gsub("`", "", substring(Covariate, sapply(gregexpr("`|\\(|\\[", Covariate), max)))),
		OR = paste0(formatC(round(exp(Estimate), digits = 2),
							format = 'f', digits = 2)),
		SE = paste0(formatC(round(exp(SE), digits = 2),
							format = 'f', digits = 2)),
		`95\\% CI` = paste0("(",
			formatC(round(exp(Estimate - qnorm(0.975) * SE), digits = 2),
							format = 'f', digits = 2),
			", ",
			formatC(round(exp(Estimate + qnorm(0.975) * SE), digits = 2),
							format = 'f', digits = 2), ")"))]
	tmp.tab[1, `:=`(Covariate = "Intercept", level = "")]
	tmp.tab[grepl("Race", level), `:=`(Covariate= "Race", level = "Black")]
	tmp.tab[grepl("Sex", level), `:=`(Covariate= "Sex", level = "Female")]
	tmp.tab[grepl("Plant", level), `:=`(Covariate = "Plant", level = 2:3)]
	tmp.tab[grepl("Employment_status", level), `:=`(Covariate = "Employment_status", level = "Left work")]
	tmp.tab[grepl("Employment_status", level), `:=`(Covariate = "Employment_status", level = "Left work")]

	tmp.tab <- rbindlist(list(
		# Intercept
		tmp.tab[1,],
		# Years since hire
		if (grepl("glm\\.rds", model.path)) {
			data.table(
			Covariate = c('Years_since_hire'),
			level = levels(dat$`Years_since_hire`)[1])},
		if (grepl("glm\\.rds", model.path)) {
		tmp.tab[grepl("Years_since_hire", Covariate),]},
		# # Calendar year
		# if (grepl("glm\\.rds", model.path)) {
		# data.table(
		# 	Covariate = c('Calendar_year'),
		# 	level = levels(dat$`Calendar_year`)[1]
		# )},
		# if (grepl("glm\\.rds", model.path)) {
		# tmp.tab[grepl("Calendar_year", Covariate),]},
		# Age
		if (grepl("glm\\.rds", model.path)) {
		data.table(
			Covariate = c('Age'),
			level = levels(dat$`Age`)[1]
		)},
		if (grepl("glm\\.rds", model.path)) {
		tmp.tab[grepl("Age", Covariate),]},
		# Plant
		data.table(
			Covariate = c('Plant'),
			level = levels(dat$Plant)[1]
		),
		tmp.tab[grepl("Plant", Covariate),],
		# Race
		data.table(
			Covariate = c('Race'),
			level = levels(dat$Race)[1]
		),
		tmp.tab[grepl("Race", Covariate),],
		# Sex
		data.table(
			Covariate = c('Sex'),
			level = levels(dat$`Sex`)[1]
		),
		tmp.tab[grepl("Sex", Covariate),],
		# Time spent machining
		data.table(
			Covariate = c('Time_spent_in_assembly'),
			level = levels(dat$`Time_spent_in_assembly`)[1]
		),
		tmp.tab[grepl("Time_spent_in_assembly", Covariate),],
		# Time spent machining
		data.table(
			Covariate = c('Time_spent_machining'),
			level = levels(dat$`Time_spent_machining`)[1]
		),
		tmp.tab[grepl("Time_spent_machining", Covariate),],
		# Time spent off
		data.table(
			Covariate = c('Time_spent_off'),
			level = levels(dat$`Time_spent_off`)[1]
		),
		tmp.tab[grepl("Time_spent_off", Covariate),],
		# Cumulative time off
		data.table(
			Covariate = c('Cumulative_time_off'),
			level = levels(dat$`Cumulative_time_off`)[1]
		),
		tmp.tab[grepl("Cumulative_time_off", Covariate),],
		# Year of hire
		data.table(
			Covariate = c('Year_of_hire'),
			level = levels(dat$`Year_of_hire`)[1]
		),
		tmp.tab[grepl("Year_of_hire", Covariate),],
		# Cumulative soluble exposure
		data.table(
			Covariate = c('Cumulative_soluble_exposure'),
			level = levels(dat$`Cumulative_soluble_exposure`)[1]
		),
		tmp.tab[grepl("Cumulative_soluble_exposure", Covariate),],
		# Cumulative straight exposure
		data.table(
			Covariate = c('Cumulative_straight_exposure'),
			level = levels(dat$`Cumulative_straight_exposure`)[1]
		),
		tmp.tab[grepl("Cumulative_straight_exposure", Covariate),],
		# Cumulative synthetic exposure
		data.table(
			Covariate = c('Cumulative_synthetic_exposure'),
			level = levels(dat$`Cumulative_synthetic_exposure`)[1]
		),
		tmp.tab[grepl("Cumulative_synthetic_exposure", Covariate),],
		# Employment status
		data.table(
			Covariate = c('Employment_status'),
			level = levels(dat$`Employment_status`)[1]
		),
		tmp.tab[grepl("Employment_status", Covariate),]
	),
	fill = T)

	tmp.tab[, level := gsub("\\[-Inf,0\\]", "0", level)]
	tmp.tab <- as.data.table(as.data.frame(tmp.tab))

	# jobloss age event counts
	n <- unlist(sapply(unique(tmp.tab$Covariate)[-1], function(x) {
		table(dat[, x, with = F])}))
	length(n); nrow(tmp.tab) - 1

	tmp.tab$n <- c("", n)

	tmp.tab <- tmp.tab[,.(
		Covariate = gsub("_", " ", Covariate),
		level,
		`$n$` = n,
		OR,
		`95\\% CI`,
		SE
	)]

	saveRDS(tmp.tab,
					file = paste0(output.dir, output.name))

	if (grepl("gam\\.rds", model.path)) {
		splines.tab[,`:=`(level = gsub("_", " ", Covariate))]
		saveRDS(splines.tab,
					file = paste0(output.dir, gsub("\\.tab", ".spline_tab", output.name)))
	}

	print(tmp.tab)
}