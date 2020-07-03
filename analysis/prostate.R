# prostate.R
# Exploring Prostate Cancer incidence by calendar year
# May 5, 2020
# Kevin Chen

# Aims ####
# Graphs calendar time (by year) on the x-axis and prostate cancer incidence on the y-axis
#  - Overall in one plot
#  - Stratified by age category at diagnosis (try 5-year windows, or maybe 10)
# 1. Run analyses in this dataset, starting in 1973.
# 2. Run analyses in this dataset, starting in 1985 and excluding the workers who had prostate cancer in SEER before 1985.
# 3. Run analyses in the original cancer dataset starting in 1985, but excluding workers who had prostate cancer in SEER before 1985 (so the plant 3 people are still misclassified)
# 4. Repeat all analyses (including the original one), using inverse weights for probability of surviving to start of follow-up (1973 or 1985)
# 5. Compare all results.

library(here)

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
	source(here::here('../gm-wrangling/wrangling', '05-Get-Exposure-Outcome.R'))
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
	cohort_analytic[, jobloss.date := yout]
	
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
	cohort_noseer[, jobloss.date := yout]
	
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
cohort_prostate <- cohort_analytic[
	sex == "M" & wh == 1 & nohist == 0 &
		possdiscr_new != 3 &
		flag77 == 0 & oddend == 0 &
		(is.na(yod) | yod >= as.Date("1973-01-01")) &
		(is.na(ddiag_pr) | ddiag_pr >= as.Date("1973-01-01")) &
		yin < as.Date("1986-01-01") &
		immortal == 0,]


# Get incidence key and other helpers
cohort2 <- as.data.frame(matrix(ncol = length(col.names)))
source(here::here("cancer incidence", "incidence.R"))

# Sort
setorder(cohort_prostate, studyno, year)
cohort_prostate[plant != 3,`:=`(
	I = 1:.N,
	N = .N
), by = .(studyno)]
# Check Ns
length(table(cohort_prostate[plant != 3]$studyno))
nrow(cohort_prostate[plant != 3 & year(ddiag_pr) >= 1973 & canc_pr == 1, ])
nrow(cohort_prostate[canc_first == 1 & plant != 3])
nrow(cohort_prostate[canc_first == 1 & (ddiag_first <= ddiag_pr | is.na(ddiag_pr)) & plant != 3])
# Any cancer other than prosate
cohort_prostate[canc_first == 1 & (ddiag_first < ddiag_pr | is.na(ddiag_pr)) & plant != 3, .(
	studyno,
	year,
	ddiag_first,
	"Not prostate" = apply(
		rbindlist(lapply(incidence.key[!code %in% c("pr", "copd", "external"), code], function(code) {
			if (code == "ma") {
				as.data.frame(matrix(canc_ma == 1 & ddiag_ma == ddiag_first & (ddiag_pr > ddiag_first | is.na(ddiag_pr)), ncol = .N))
			} else {as.data.frame(matrix(get(paste0("canc_", code)) == 1 & get(paste0("ddiag_", code)) == ddiag_first, ncol = .N))}
		})), 2, function(x) {as.numeric(sum(x) > 0)}),
	canc_pr, canc_ma, canc_which_first, ddiag_pr, ddiag_ma
)][`Not prostate` == 1]#[sapply(canc_which_first, function(x) {"ma" %in% x & !"pr" %in% x})]
# All cancers (some individuals represented more than once)
cohort_prostate[I == N & plant != 3, .(Site = {
	site <- unlist(canc_which_first)
	if (is.null(site)) {as.character(NA)} else {incidence.key$description[match(site, incidence.key$code)]}
}),
by = .(studyno)]$Site %>%
	table %>% as.data.frame %>% arrange(-Freq)

# No plant 3
cohort_prostate <- cohort_prostate[plant != 3]

# Get reshaped analytic dataset for prostate cancer ####
get.coxph(
	cohort_name = "cohort_prostate",
	outcomes = which(incidence.key$code == "pr"),
	start.year = 1973,
	time_scale = "year", run_model = F)

# Make yoi with class numeric
pr.dat[!is.na(yoi), `:=`(
	yoi.gm = date.to.gm(yoi))]

# Fewer age categories
pr.dat <- pr.dat[Sex == "M"]
pr.dat[,`:=`(Age = cut(age.year2/365.25,
											 c(-Inf, 60, 70, 80, Inf) - 1))]

# Check levels
pr.dat$Race %>% levels

# Pretty Age cagetories
pr.dat$Age %>% table
pr.dat[, `:=`(Age = factor(
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

# Overall incidence
overall.tab <- pr.dat[Plant != 3, .(
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
age_specific.tab <- pr.dat[Plant != 3, .(
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
	geom_rug(data = pr.dat[status == 1], aes(
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
	geom_rug(data = pr.dat[status == 1], aes(
		x = yoi.gm
	), size = 0.01) +
	mytheme + theme(
		legend.position = "bottom",
		legend.margin = margin()
	) + guides(color = guide_legend(nrow = 2))

# # Render
# tikz(here::here(paste0("/Users/kevinchen/eisen/GM/left truncation/resources", "FU from ", min(pr.dat$year)), "prostate-cancer_all-ages.tex"),
# 		 standAlone = T, height = 2, width = 3)
# print(overall.ggplot)
# dev.off()
#
# tikz(here::here(paste0("/Users/kevinchen/eisen/GM/left truncation/resources", "FU from ", min(pr.dat$year)), "prostate-cancer_by-age.tex"),
# 		 standAlone = T, height = 2.5, width = 3)
# print(by_age.ggplot)
# dev.off()
#
# # Compile
# lualatex(directory = here::here(paste0("/Users/kevinchen/eisen/GM/left truncation/resources", "FU from ", min(pr.dat$year))))

# Model 1: Plants 2 and 3 ####
# Run analyses in this dataset, just like I did for the original cancer cohort but starting in 1973.
pr1.dat <- as.data.table(as.data.frame(pr.dat))
pr1.dat[,`:=`(
	Race = relevel(factor(Race), "White"),
	Plant = relevel(factor(Plant), "1")
)]

# tmp.coxph <- coxph(
# 	Surv(age.year1, age.year2, status) ~
# 		`Cumulative straight` +
# 		`Cumulative soluble 5` +
# 		`Cumulative synthetic` +
# 		pspline(year, df = 0)  +
# 		pspline(yin.gm, df = 0) +
# 		Race + Plant,
# 	data = pr1.dat,
# 	method = "efron")
#
# # Save model  ####
# saveRDS(tmp.coxph,
# 				file = paste0(
# 					to_drive_D(here::here('left truncation/resources/')),
# 					"pr_mod1.coxph.rds"
# 				))
#
# get.coef(outcomes = which(incidence.key$code == "pr"),
# 				 cohort_name = "cohort_prostate",
# 				 analytic.name = "pr1.dat",
# 				 new_dat = F,
# 				 messy_sol = NULL,
# 				 mod.directory = to_drive_D(here::here('left truncation/resources/')),
# 				 mod.name = "pr_mod1.coxph.rds",
# 				 directory = here::here('reports/left truncation/prostate cancer/resources/'),
# 				 file.prefix = "pr_mod1")
#
# # Descriptive table
# source(here::here("reports", "table1.R"))
# popchar1.tab1 <- get.tab1(df = pr1.dat, exposure_lag = 20, incidence = T)
# saveRDS(popchar1.tab1,
# 				file = here::here(
# 					'left truncation/resources/', "popchar1.tab1.rds"))

# Model 2: Start in 1985, exclude those with prostate cancer ####
# Run analyses in this dataset starting in 1985, excluding the workers who had prostate cancer in SEER before 1985.
length(table(pr.dat[,][yoi < as.Date("1985-01-01"), studyno]))
pr2.dat <- as.data.table(as.data.frame(pr.dat[(is.na(yoi) | yoi >= as.Date("1985-01-01")) &
																								year >= 1985,]))

pr2.dat[,`:=`(
	Race = relevel(factor(Race), "White"),
	Plant = relevel(factor(Plant), "1")
)]

# tmp.coxph <- coxph(
# 	Surv(age.year1, age.year2, status) ~
# 		`Cumulative straight` +
# 		`Cumulative soluble 5` +
# 		`Cumulative synthetic` +
# 		pspline(year, df = 0)  +
# 		pspline(yin.gm, df = 0) +
# 		Race + Plant,
# 	data = pr2.dat,
# 	method = "efron")
#
# # Save model  ####
# saveRDS(tmp.coxph,
# 				file = paste0(
# 					to_drive_D(here::here('left truncation/resources/')),
# 					"pr_mod2.coxph.rds"
# 				))
#
# get.coef(outcomes = which(incidence.key$code == "pr"),
# 				 cohort_name = "cohort_prostate",
# 				 analytic.name = "pr2.dat",
# 				 new_dat = F,
# 				 messy_sol = NULL,
# 				 mod.directory = to_drive_D(here::here('left truncation/resources/')),
# 				 mod.name = "pr_mod2.coxph.rds",
# 				 directory = here::here('reports/left truncation/prostate cancer/resources/'),
# 				 file.prefix = "pr_mod2")
#
# # Descriptive table
# popchar2.tab1 <- get.tab1(df = pr2.dat, exposure_lag = 20, incidence = T)
# saveRDS(popchar2.tab1,
# 				file = here::here('left truncation/resources/',
# 													"popchar2.tab1.rds"))


# Model 3 ####
# Run analyses starting in 1985, but excluding workers who had prostate cancer in SEER before 1985 WITH WEIGHTS
length(table(pr2.dat[,][yoi < as.Date("1985-01-01"), studyno]))
pr3.dat <- as.data.table(as.data.frame(pr2.dat))

# Get weights
weights.year.min <- 1941
if (weights.year.min == 1941) {
	# Survival: 1941 to 1985
	survival_to_85 <- box_read(674612305081)}
if (weights.year.min == 1970) {
	# Survival: 1970 to 1985
	survival_to_85 <- box_read(671315790911)
}
setDT(survival_to_85)
pr3.dat <- merge(pr3.dat,
								 survival_to_85[dead_by_85 == 0, .(
								 	studyno,
								 	p = big.gam.prob)],
								 on = "studyno",
								 all.x = T)

# If hired in 1985, go ahead and give them full weight
pr3.dat[is.na(p) & year(yin) + 3 > 1985, p := 1]

# tmp.coxph <- coxph(
# 	Surv(age.year1, age.year2, status) ~
# 		`Cumulative straight` +
# 		`Cumulative soluble 5` +
# 		`Cumulative synthetic` +
# 		pspline(year, df = 0)  +
# 		pspline(yin.gm, df = 0) +
# 		Race + Plant,
# 	weights = p,
# 	cluster = studyno,
# 	data = pr3.dat,
# 	method = "efron")
#
# # Save model  ####
# saveRDS(tmp.coxph,
# 				file = paste0(
# 					to_drive_D(here::here('left truncation/resources/')),
# 					"pr_mod3.coxph.rds"
# 				))
#
# get.coef(outcomes = which(incidence.key$code == "pr"),
# 				 cohort_name = "cohort_prostate",
# 				 analytic.name = "pr3.dat",
# 				 new_dat = F,
# 				 messy_sol = NULL,
# 				 mod.directory = to_drive_D(here::here('left truncation/resources/')),
# 				 mod.name = "pr_mod3.coxph.rds",
# 				 directory = here::here('reports/left truncation/prostate cancer/resources/'),
# 				 file.prefix = "pr_mod3")

# Descriptive table
popchar3.tab1 <- get.tab1(df = pr3.dat, exposure_lag = 20, incidence = T)
saveRDS(popchar3.tab1,
				file = here::here('left truncation/resources/',
													"popchar3.tab1.rds"))

# Model 4 ####
# Run analyses starting in 1985, but excluding workers who had prostate cancer in SEER before 1985 WITH STABILIZED WEIGHTS
overall.p <- 1 - mean(pr3.dat[,.(p = p[1]), by = .(studyno)]$p)
pr3.dat[,`:=`(
	sw = overall.p/(1 - p)
), by = .(studyno)]

# Truncate at 90th percentile
sw.max <- quantile(pr3.dat[, .(sw = sw[1]), by = studyno]$sw,
									 .99)
pr3.dat[sw > sw.max,`:=`(
	sw = sw.max
), by = .(studyno)]

pr3.dat[, .(
	sw = sw[1],
	age = age.year2[year == 1985]/365
), by = studyno]$sw %>% summary

(pr3.dat[year == 1985]$age.year2/365) %>% summary

pr3.dat[year == 1985,`:=`(
	age1985.cat = cut(age.year2/365,
										c(25, 48, 104),
										include.lowest = T)
)]

pr3.dat[, .(
	sw = sw[1],
	age = age.year2[year == 1985]/365,
	`year of hire` = yin.gm[1],
	`cumulative soluble` = cum_soluble[year == 1985],
	age.cat = age1985.cat[year == 1985]
), by = studyno][`cumulative soluble` > 0] %>% ggplot(
	aes(x = `cumulative soluble`, y = sw)) + geom_point(
		alpha = 0.6) +
	# coord_cartesian(xlim = c(0, 60)) +
	mytheme.web

# tmp.coxph <- coxph(
# 	Surv(age.year1, age.year2, status) ~
# 		`Cumulative straight` +
# 		`Cumulative soluble 5` +
# 		`Cumulative synthetic` +
# 		pspline(year, df = 0)  +
# 		pspline(yin.gm, df = 0) +
# 		Race + Plant,
# 	weights = sw,
# 	cluster = studyno,
# 	data = pr3.dat,
# 	method = "efron")
#
# summary(tmp.coxph)
#
# # Save model  ####
# saveRDS(tmp.coxph,
# 				file = paste0(
# 					to_drive_D(here::here('left truncation/resources/')),
# 					"pr_mod4.coxph.rds"
# 				))
#
# get.coef(outcomes = which(incidence.key$code == "pr"),
# 				 cohort_name = "cohort_prostate",
# 				 analytic.name = "pr3.dat",
# 				 new_dat = F,
# 				 messy_sol = NULL,
# 				 mod.directory = to_drive_D(here::here('left truncation/resources/')),
# 				 mod.name = "pr_mod4.coxph.rds",
# 				 directory = here::here('reports/left truncation/prostate cancer/resources/'),
# 				 file.prefix = "pr_mod4")

# Model 5 ####
# Run analyses starting in 1985, but excluding workers who had prostate cancer in SEER before 1985 with inverse prob of being ALIVE
# overall.p <- mean(pr3.dat$p)
pr3.dat[,`:=`(
	alive.sw = overall.p/(p)
), by = .(studyno)]

ggplot(pr3.dat[, .(sw = alive.sw[1]), by = studyno], aes(x = sw)) + geom_histogram() + mytheme.web

pr3.dat[,`:=`(
	sw = c(sw[1], rep(1, .N - 1))
), by = .(studyno)]

# tmp.coxph <- coxph(
# 	Surv(age.year1, age.year2, status) ~
# 		`Cumulative straight` +
# 		`Cumulative soluble 5` +
# 		`Cumulative synthetic` +
# 		pspline(year, df = 0)  +
# 		pspline(yin.gm, df = 0) +
# 		Race + Plant,
# 	weights = sw,
# 	cluster = studyno,
# 	data = pr3.dat,
# 	method = "efron")
#
# # Save model  ####
# saveRDS(tmp.coxph,
# 				file = paste0(
# 					to_drive_D(here::here('left truncation/resources/')),
# 					"pr_mod5.coxph.rds"
# 				))
#
# get.coef(outcomes = which(incidence.key$code == "pr"),
# 				 cohort_name = "cohort_prostate",
# 				 analytic.name = "pr3.dat",
# 				 new_dat = F,
# 				 messy_sol = NULL,
# 				 mod.directory = to_drive_D(here::here('left truncation/resources/')),
# 				 mod.name = "pr_mod5.coxph.rds",
# 				 directory = here::here('reports/left truncation/prostate cancer/resources/'),
# 				 file.prefix = "pr_mod5")