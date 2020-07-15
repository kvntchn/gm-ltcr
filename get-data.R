# get-data.R
# Get data for estimating survival in each year up to 1985
# Kevin Chen
# May 22, 2020

library(here)

if (!('cohort' %in% ls(envir = .GlobalEnv))) {
source(here::here('../gm-wrangling/wrangling', '00-Hello.R'))
}

# Get person-year dataset ####
# rm(cohort_analytic)
if (!"cohort_analytic" %in% ls()) {
	cohort_analytic <- get.cohort_analytic(
		outcome_type = "mortality",
		exposure.lag = 0,
		deathage.max = NULL,
		end.year = 1984,
		hire.year.min = 1938
	)
	setorder(cohort_analytic, studyno, year)

	# Drop unnecessary data ####
	cohort_analytic <- cohort_analytic[
		year <= 1984 &
			wh == 1 & nohist == 0 & immortal == 0 &
			right.censored == 0 &
			yin >= as.Date("1938-01-01") &
			yin < as.Date("1983-01-01"),]

	# Indicator for at work or not
	cohort_analytic[year <= year(yout), `:=`(`At work` = 1)]
	cohort_analytic[year > year(yout), `:=`(`At work` = 0)]
	cohort_analytic[year > year(yout), `:=`(
		machining.gan = 0,
		machining.han = 0,
		machining.san = 0,
		assembly.gan = 0,
		assembly.han = 0,
		assembly.san = 0,
		off.gan = 0,
		off.han = 0,
		off.san = 0)]
	cohort_analytic[is.na(machining.gan) & year == year(yout), `:=`(
		machining.gan = 0,
		machining.han = 0,
		machining.san = 0,
		assembly.gan = 0,
		assembly.han = 0,
		assembly.san = 0,
		off.gan = 0,
		off.han = 0,
		off.san = 0)]

	# years since hire variable
	cohort_analytic[,`:=`(sincehire.years = 1:.N + 3), by = .(studyno)]

	# pretty yin
	cohort_analytic[,`:=`(yin16 = yin16 + 1900)]

	# Clean up work history variables
	cohort_analytic[, `:=`(
		off = apply(data.frame(off.gan + off.han + off.san, 1), 1, min),
		machining = apply(data.frame(machining.gan + machining.han + machining.san, 1), 1, min),
		assembly = apply(data.frame(assembly.gan + assembly.han + assembly.san, 1), 1, min)
	)]

	# Fill in for years 1993 and 1994 and get cumulative off
	cohort_analytic[, `:=`(
		off = zoo::na.locf(off),
		machining = zoo::na.locf(machining),
		assembly = zoo::na.locf(assembly),
		cumulative_off = cumsum(zoo::na.locf(off))
	), by = .(studyno)]

	# outcome variable
	cohort_analytic[, y := `All natural causes`]
}

length(table(cohort_analytic$studyno)); sum(cohort_analytic$y)

# Make categorical variables for ease ####
# probs <- seq(0, 1, 0.25)
probs <- seq(0, 1, 0.1)

if ("Age" %in% names(cohort_analytic)) {
	cohort_analytic <- cohort_analytic[,-"Age", with = F]}
if ("Time_spent_machining" %in% names(cohort_analytic)) {
	cohort_analytic <- cohort_analytic[, -c(
		"Years_since_hire",
		"Calendar_year",
		"Time_spent_in_assembly",
		"Time_spent_machining",
		"Time_spent_off",
		"Cumulative_time_off",
		"Year_of_hire",
		"Cumulative_soluble_exposure",
		"Cumulative_straight_exposure",
		"Cumulative_synthetic_exposure",
		"Employment_status"
	), with = F]
}

cohort_analytic[,`:=`(
	age = age.year2/365,
	Employment_status = factor(ifelse(year <= year(yout), "At work", "Left work"),
														 levels = c("At work", "Left work")),
	Plant = factor(plant, levels = as.character(1:3)),
	Race = factor(race, levels = c("White", "Not white"), labels = c("White", "Black")),
	Sex = factor(sex, levels = c("M", "F"), labels = c("Male", "Female"))
)]

covariate.quantile <- apply(
	cohort_analytic[, .(
		sincehire.years, year, age,
		assembly, machining, off, cumulative_off,
		yin16, cum_soluble, cum_straight, cum_synthetic
	)], 2, function(x) {
		quantile.tmp <- quantile(x[cohort_analytic[, y == 1]], probs)
		if (quantile.tmp[1] == quantile.tmp[2] &
				quantile.tmp[1] == 0) {quantile.tmp[1] <- -Inf} else {
					quantile.tmp[1] <- min(x)
				}
		quantile.tmp[length(quantile.tmp)] <- max(x)
		return(quantile.tmp)
	})

cohort_analytic[,(c(
	"Years_since_hire",
	"Calendar_year",
	"Age",
	"Time_spent_in_assembly",
	"Time_spent_machining",
	"Time_spent_off",
	"Cumulative_time_off",
	"Year_of_hire",
	"Cumulative_soluble_exposure",
	"Cumulative_straight_exposure",
	"Cumulative_synthetic_exposure"
)):=lapply(
	colnames(covariate.quantile), function(x) {
		x <- cut(get(x), unique(covariate.quantile[,x]), include.lowest = T, dig.lab = 4)
		levels(x) <- gsub("\\[-Inf,0\\]", "0", levels(x))
		x
	})]

setorder(cohort_analytic, studyno, year)