# cumulative_incidence.R
# Parametric estimation of cumulative incidence
# June 15, 2020
# Kevin Chen

# Plan ####
# 1. Use pooled logistic regression to crudely estimate the J = 2 subdistribution hazard functions (j = 1 for cancer incidence and j = 2 for natural cause mortality), given covariates W and exposure X
# 	h_j (t | W, X) = lim_{δ -> 0} 1/δ * Pr{ t < T < t + δ, J = j | T > T or (T < t and J ≠ j), W, X }
# 		≈ 1/1 * Pr{ t < T < t + 1, J = j | T > T or (T < t and J ≠ j), W, X }
# 2. For each person-year, get p_j by evaluating the jth model at (X = x, W = w) where x is the counter-factual exposure level of interest, and w is the observed covariate combination
# 3. For each person-year, take the the product p_1 * p_2
# 	(Would this to be something like a model-based estimate of S(t)h(t) where h(t) is the hazard of cancer incidence weighted by that for natural cause mortality, marginal over covariates for a particular exposure level x? Since folks do not contribute person time after experiencing the event, can we take S(t) = 1 for at-risk time and S(t) = 0 for not at-risk time?)
# 4. Take the sum of products
# 	(Would this be a rough estimate of some weighted cumulative incidence ∫S(t)h(t)dt ?)
# 5. Take the ratio of (1) the cumulative incidence estimate under an exposure level of interest to (2) that for the reference level of exposure.


if (!grepl("gm", getwd(), ignore.case = T)) {
	if ("here" %in% .packages()) {
		detach("package:here", unload = T)
	}
	setwd('eisen/gm')
	library(here)
}

# Table 1 helper function
source(here::here("reports/table1.R"))

# Get analytic data ####
if (!("exposure.lag" %in% ls())) {
	exposure.lag <- 20
	end.year <- 2015
	# rm(cohort_analytic)
}

# Get cohort analytic ####
if (
	!('cohort_analytic' %in% ls())
) {
	outcome.type <- 'incidence'
	source(here::here('../gm-wrangling/', '05-Get-Exposure-Outcome.R'))
	cohort_analytic <- get.cohort_analytic(
		outcome_type = outcome.type,
		exposure.lag = exposure.lag,
		deathage.max = NULL,
		end.year = end.year,
		hire.year.min = 1938,
		use_seer = T
	)
	setorder(cohort_analytic, studyno, year)
	cohort_analytic[, `:=`(yin.gm = date.to.gm(yin))]

	# # Keep only people who appear in the exposure data
	# cohort_analytic <-
	# 	cohort_analytic[studyno %in% unique(exposure$studyno)]

	# PICK YOUT ####
	cohort_analytic[, jobloss.date := employment_end.date]

	# Exposure after leaving work is 0
	cohort_analytic[year > (year(jobloss.date) + exposure.lag), `:=`(
		straight = 0,
		soluble = 0,
		synthetic = 0)]
	# NA fill
	cohort_analytic[year <= (year(jobloss.date) + exposure.lag), `:=`(
		straight = zoo::na.locf(straight),
		soluble = zoo::na.locf(soluble),
		synthetic = zoo::na.locf(synthetic)
	), by = .(studyno)]
	cohort_analytic[, `:=`(
		cum_straight = cumsum(straight),
		cum_soluble = cumsum(soluble),
		cum_synthetic = cumsum(synthetic)
	), by = .(studyno)]

	# Which columns ####
	col.names <- c(
		"studyno",
		"age.year1",
		"age.year2",
		"year1",
		"year2",
		grep("canc\\_", names(cohort_analytic), value = T),
		"straight",
		"cum_straight",
		"soluble",
		"cum_soluble",
		"synthetic",
		"cum_synthetic",
		"year",
		"yin.gm",
		"yin",
		"yrin",
		"yrin16",
		"race",
		"finrace",
		"plant",
		grep("ddiag", names(cohort_analytic), value = T),
		"yod",
		"yoc",
		"yob",
		"sex",
		"dateout.date",
		"employment_end.date",
		"employment_end.date.legacy",
		"yout",
		"yout_recode",
		"jobloss.date",
		"All causes",
		"Chronic obstructive pulmonary disease",
		"All external causes",
		"nohist", "wh", "immortal", "right.censored",
		"possdiscr_new", "flag77", "oddend", "status15")

	# Drop unnecessary data ####
	cohort_analytic <- cohort_analytic[
		wh == 1 & nohist == 0 & immortal == 0 & right.censored == 0,
		col.names, with = F]
}

# Get no seer data ####
if (
	!('cohort_noseer' %in% ls())
) {
	cohort_noseer <- get.cohort_analytic(
		outcome_type = outcome.type,
		exposure.lag = exposure.lag,
		deathage.max = NULL,
		end.year = end.year,
		hire.year.min = 1938,
		use_seer = F
	)
	setorder(cohort_noseer, studyno, year)
	cohort_noseer[, `:=`(yin.gm = date.to.gm(yin))]

	# # Keep only people who appear in the exposure data
	# cohort_noseer <-
	# 	cohort_noseer[studyno %in% unique(exposure$studyno)]

	# PICK YOUT ####
	cohort_noseer[, jobloss.date := employment_end.date]

	# Exposure after leaving work is 0
	cohort_noseer[year > (year(jobloss.date) + exposure.lag), `:=`(
		straight = 0,
		soluble = 0,
		synthetic = 0)]
	# NA fill
	cohort_noseer[year <= (year(jobloss.date) + exposure.lag), `:=`(
		straight = zoo::na.locf(straight),
		soluble = zoo::na.locf(soluble),
		synthetic = zoo::na.locf(synthetic)
	), by = .(studyno)]
	cohort_noseer[, `:=`(
		cum_straight = cumsum(straight),
		cum_soluble = cumsum(soluble),
		cum_synthetic = cumsum(synthetic)
	), by = .(studyno)]

	# Drop unnecessary data ####
	cohort_noseer <- cohort_noseer[
		wh == 1 & nohist == 0 & immortal == 0 & right.censored == 0,
		col.names, with = F]
}

# Restrict cohort to those still alive at the start of 1973
cohort_cancer <- cohort_analytic[
	plant != 3 &
	sex == "M" & wh == 1 & nohist == 0 &
		possdiscr_new != 3 &
		flag77 == 0 & oddend == 0 &
		(is.na(yod) | yod >= as.Date("1973-01-01")) &
		(is.na(ddiag_first) | ddiag_first >= as.Date("1973-01-01")) &
		yin < as.Date("1986-01-01") &
		immortal == 0,]


# Get incidence key and other helpers
cohort2 <- as.data.frame(matrix(ncol = length(col.names)))
source(here::here("cancer incidence", "incidence.R"))

# Sort
setorder(cohort_cancer, studyno, year)
cohort_cancer[plant != 3,`:=`(
	I = 1:.N,
	N = .N
), by = .(studyno)]

# Check Ns
length(table(cohort_cancer$studyno))
nrow(cohort_cancer[year(ddiag_first) >= 1973 & canc_first == 1, ])
nrow(cohort_cancer[canc_first == 1])
nrow(cohort_cancer[canc_first == 1 & (ddiag_first <= ddiag_first | is.na(ddiag_first))])

# Get reshaped analytic dataset for cancer ####
if (!"first" %in% incidence.key$code) {
	incidence.key <- rbindlist(list(
		incidence.key,
		data.table(
			code = "first",
			description = "All cancers",
			var.name = "canc_first",
			date.name = "ddiag_first")))
	}

get.coxph(
	cohort_name = "cohort_cancer",
	outcomes = which(incidence.key$code == "first"),
	start.year = 1973,
	time_scale = "year", run_model = F)

cancer.dat <- as.data.frame(first.dat)
setDT(cancer.dat)
rm(first.dat)


# Make yoi with class numeric
cancer.dat[!is.na(yoi), `:=`(
	yoi.gm = date.to.gm(yoi))]

# Fewer age categories
cancer.dat <- cancer.dat[Sex == "M"]
cancer.dat[,`:=`(Age = cut(age.year2/365.25,
											 c(-Inf, 60, 70, 80, Inf) - 1))]

# Check levels
cancer.dat$Race %>% levels

# Pretty Age cagetories
cancer.dat[, `:=`(Age = factor(
	Age, levels = levels(Age), labels = {
		# Get upper and lower bounds (delimited by comma)
		lower <- as.numeric(substr(levels(Age), 2, unlist(gregexpr(",", levels(Age))) - 1))
		upper <- as.numeric(substr(levels(Age), unlist(gregexpr(",", levels(Age))) + 1,
															 nchar(levels(Age)) - 1))

		# Make new levels by pasting bounds
		new.levels <- paste(lower + 1, "to", upper)

		# Correct endpoints
		new.levels[!is.finite(upper)] <- paste(lower[!is.finite(upper)] + 1, "and older")
		new.levels[!is.finite(lower)] <- paste(upper[!is.finite(lower)], "and younger")

		new.levels
	}))]
cancer.dat$Age %>% table

# Overall incidence
overall.tab <- cancer.dat[Plant != 3, .(
	N = length(table(studyno)),
	py = sum(year2 - year1),
	cases = sum(status),
	rate = sum(status)/(sum(year2 - year1)/365.25),
	risk = sum(status)/length(table(studyno))
), by = .(year)]
# Get CIs
overall.tab[,`:=`(
	rate.upper = exp(log(rate) + qnorm(0.975) * sqrt(rate + (1 - rate)/py)),
	rate.lower = exp(log(rate) - qnorm(0.975) * sqrt(rate + (1 - rate)/py)),
	risk.upper = risk + qnorm(0.975) * sqrt(risk + (1 - risk)/N),
	risk.lower = risk - qnorm(0.975) * sqrt(risk + (1 - risk)/N)
)]

# Age-specific stuff
age_specific.tab <- cancer.dat[Plant != 3, .(
	N = length(table(studyno)),
	py = sum(year2 - year1),
	cases = sum(status),
	rate = sum(status)/(sum(year2 - year1)/365.25),
	risk = sum(status)/length(table(studyno))
), by = .(year, Age)]
# Get CIs
age_specific.tab[,`:=`(
	rate.upper = exp(log(rate) + qnorm(0.975) * sqrt(rate + (1 - rate)/py)),
	rate.lower = exp(log(rate) - qnorm(0.975) * sqrt(rate + (1 - rate)/py)),
	risk.upper = risk + qnorm(0.975) * sqrt(risk + (1 - risk)/N),
	risk.lower = risk - qnorm(0.975) * sqrt(risk + (1 - risk)/N)
)]

# Make plots ####
overall.ggplot <- ggplot() +
	geom_line(data = overall.tab, aes(
		x = year,
		y = rate * 100000), size = 0.5) +
	labs(x = "Calendar year",
			 y = "Incidence (per 100,000 p$\\cdot$years)") +
	geom_rug(data = cancer.dat[status == 1], aes(
		x = yoi.gm
	), size = 0.01) +
	mytheme

by_age.ggplot <- ggplot() +
	geom_line(data = age_specific.tab, aes(
		x = year,
		y = rate * 100000,
		color = Age), size = 0.5) +
	# geom_smooth(se = T, span = 0.6, size = 0.2, alpha = 0.2) +
	# facet_wrap(. ~ Age, ncol = 2) +
	labs(x = "Calendar year",
			 y = "Incidence (per 100,000 p$\\cdot$years)") +
	geom_rug(data = cancer.dat[status == 1], aes(
		x = yoi.gm
	), size = 0.01) +
	mytheme + theme(
		legend.position = "bottom",
		legend.margin = margin()
	) + guides(color = guide_legend(nrow = 2))

# # Render
# tikz(here::here(paste0("/Users/kevinchen/eisen/GM/left truncation/resources", "FU from ", min(cancer.dat$year)), "all-cancers_all-ages.tex"),
# 		 standAlone = T, height = 2, width = 3)
print(overall.ggplot)
# dev.off()
#
# tikz(here::here(paste0("/Users/kevinchen/eisen/GM/left truncation/resources", "FU from ", min(cancer.dat$year)), "all-cancers_by-age.tex"),
# 		 standAlone = T, height = 2.5, width = 3)
print(by_age.ggplot)
# dev.off()
#
# # Compile
# lualatex(directory = here::here(paste0("/Users/kevinchen/eisen/GM/left truncation/resources", "FU from ", min(cancer.dat$year))))

# Add rows after death ####
# Make index
setorder(cancer.dat, studyno, year)
cancer.dat$I <- NA; cancer.dat$N <- NA
cancer.dat[is.na(yoi), `:=`(
	I = 1:.N,
	N = .N
	), by = .(studyno)]

# Make rows to add, carrying forward last covariate combination
to_add <- merge(
	# Rows to add
	cancer.dat[I == N, .(
		year,
		studyno)][year != end.year, .(
			year = (year + 1):2015
		), by = .(studyno)],
	# Covariates in last row
	cancer.dat[I == N, -"year"],
	on = "studyno",
	all.x = T)

# Add rows
cancer.dat <- rbindlist(list(
	cancer.dat,
	to_add
), use.names = T, idcol = "subdistribution")
setorder(cancer.dat, studyno, year)
cancer.dat[, subdistribution := (0:1)[as.numeric(factor(subdistribution, levels = 1:2))]]

# Check levels
cancer.dat[,`:=`(
			Race = relevel(factor(Race), "White"),
			Plant = relevel(factor(Plant), "1"),
			Age = cut(age.year2/365,
								unique(c(min(age.year2/365),
									quantile(cancer.dat[status == 1, age.year2/365],
													 c(0.25, 0.5, 0.75)),
									max(age.year2/365))),
								include.lowest = T)
			)]

# # Descriptive table
# source(here::here("reports", "table1.R"))
# popchar1.tab1 <- get.tab1(
# 	df = cancer.dat[subdistribution == 0], exposure_lag = 20, incidence = T)
# saveRDS(popchar1.tab1,
# 				file = here::here(
# 					'left truncation/resources/cif', "popchar1.tab1.rds"))
#

# Run PLR ####
library(mgcv)
names(cancer.dat) <- gsub(" ", "_", names(cancer.dat))
# tmp.plr <- mgcv::gam(status ~
# 			Cumulative_straight +
# 			Cumulative_soluble_5 +
# 			Cumulative_synthetic +
# 			 + Plant + Race + Year_of_hire +
# 			Age +
# 			s(year, bs = "ps") +
# 			s(yin.gm, bs = "ps"),
# 		data = cancer.dat)
#
# # Save PLR
# saveRDS(tmp.plr,
# 				file = paste0(
# 					to_drive_D(here::here('left truncation/resources/')),
# 					"cancer.plr.rds"
# 				))

readRDS(paste0(
					to_drive_D(here::here('left truncation/resources/')),
					"cancer.plr.rds"))

# Get coefficient table
get.coef(outcomes = which(incidence.key$code == "first"),
				 cohort_name = "cancer.dat",
				 analytic.name = "cancer.dat",
				 new_dat = F,
				 mod.directory = to_drive_D(here::here('left truncation/resources/')),
				 mod.name = "cancer.plr.rds",
				 directory = here::here('reports/left truncationweighted cumulative incidence/resources/'),
				 file.prefix = "cancer.plr")

# Get survival estimates ####
survival_by_year <- box_read(679752455177)
setDT(survival_by_year)
