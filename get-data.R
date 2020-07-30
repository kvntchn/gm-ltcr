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

# Refactoring and other transformations
cohort_analytic[,`:=`(
	age = age.year2/365,
	Age.exp = exp(age.year2/365),
	Employment_status = factor(ifelse(year <= year(yout), "At work", "Left work"),
														 levels = c("At work", "Left work")),
	Plant = factor(plant, 1:3, paste("Plant", 1:3)),
	Race = factor(race, levels = c("White", "Black")),
	Sex = factor(sex, levels = c("M", "F"), labels = c("Male", "Female")),
	employment.years = time_length(difftime(
		apply(data.frame(
			as.Date(paste0(year , "-12-31")),
			yoc,
			yod,
			yout
		), 1, min, na.rm = T),
		yin,
	), 'year')
)]
# cohort_analytic[year > year(yout), .(year, Duration_of_employment, yin, yout), by = .(studyno)]

# Make categorical variables ####
categorize.cols <- c(
	"employment.years", "sincehire.years", "year", "age",
	"assembly", "machining", "off", "cumulative_off",
	"yin16", "cum_soluble", "cum_straight", "cum_synthetic")
new_cat.cols <- c(
	"Duration_of_employment", "Years_since_hire", "Calendar_year", "Age",
	"Time_spent_in_assembly", "Time_spent_machining", "Time_spent_off", "Cumulative_time_off",
	"Year_of_hire", "Cumulative_soluble_exposure", "Cumulative_straight_exposure", "Cumulative_synthetic_exposure"
)

probs <- list(
	# Fewer categories
	sapply(categorize.cols, function(x) {
		if (x %in% paste0("cum_", c("soluble", "straight", "synthetic"))) {
			c(0, 2/3, 3/4, 1)
		} else {seq(0, 1, 1/4)}
	}, simplify = F),
	# More categories (deciles)
	"_deciles" = sapply(categorize.cols, function(x) {
		seq(0, 1, 1/10)
	}, simplify = F))

# Make new categorical variables for each cutpoint set of interest
for (i in 1:length(probs)) {
	
	# Remove "new" categorical columns if they were made earlier
	# This would happen if this script were already sourced previously
	sapply(paste0(new_cat.cols, names(probs)[i]), function(remove.cols) {
		if (remove.cols %in% names(cohort_analytic)) {
			cohort_analytic <<- cohort_analytic[, -remove.cols, with = F]
		}}) 
	
	# Get covariate cutpoints (or quantiles)
	covariate.cutpoints <- sapply(categorize.cols, function(cat.col = categorize.cols[1]) {
		x <- unlist(cohort_analytic[, cat.col, with = F])
		quantile.tmp <- quantile(x[cohort_analytic$y == 1], probs[[i]][[cat.col]])
		if (quantile.tmp[1] == quantile.tmp[2] & quantile.tmp[1] == 0) {
			quantile.tmp[1] <- -Inf} else {
				quantile.tmp[1] <- min(x)
			}
		if (cat.col %in% c("cumulative_off", "cum_soluble", "cum_straight", "cum_synthetic")) {
			quantile.tmp <- c(-Inf, 0, quantile.tmp)
		}
		quantile.tmp[length(quantile.tmp)] <- max(x)
		return(sort(unique(quantile.tmp)))
	}, simplify = F)
	
	# User-specified cutpoints
	# Keeping 1st and last, but changing middle cutpoints
	change.middle <- function(var.name, ...) {
		sort(c(..., covariate.cutpoints[[var.name]][c(1, length(covariate.cutpoints[[var.name]]))]))
	}
	
	if (length(table(sapply(probs[[i]], length))) != 1) {
		covariate.cutpoints$year <- change.middle("year", 1970, 1980)
		covariate.cutpoints$age <- change.middle('age', 55, 70)
		covariate.cutpoints$yin16  <- change.middle("yin16", 1945, 1960)
		covariate.cutpoints$employment.years <- change.middle("employment.years", 10, 20)
		covariate.cutpoints$cumulative_off <- change.middle("cumulative_off", 0, 1)
	}
	
	cohort_analytic[, (paste0(new_cat.cols, names(probs)[i])):=lapply(
		names(covariate.cutpoints), function(x = names(covariate.cutpoints)[1]) {
			x <- cut(get(x), unique(covariate.cutpoints[[x]]), include.lowest = T, dig.lab = 4)
			levels(x) <- gsub("\\[-Inf,0\\]", "0", levels(x))
			x
		})]
}

setorder(cohort_analytic, studyno, year)