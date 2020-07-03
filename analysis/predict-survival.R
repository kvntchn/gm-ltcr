# predict-survival.R
# Pooled logistic regression for estimating survival in each year up to 1985
# Kevin Chen
# May 22, 2020

# Get person-year dataset ####
source(here::here('../gm-wrangling/wrangling', '05-Get-Exposure-Outcome.R'))

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
probs <- seq(0, 1, 0.25)
# probs <- seq(0, 1, 0.05)

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

# # Fit models ####
# mod.gam <- mgcv::gam(
# 	y ~ s(sincehire.years, bs='ps') +
# 		# s(year, bs='ps') +
# 		s(age.year2, bs='ps') +
# 		# `Years_since_hire` +
# 		# `Calendar_year` +
# 		# `Age` +
# 		Race + Plant + Sex +
# 		`Time_spent_in_assembly` +
# 		`Time_spent_machining` +
# 		`Time_spent_off` +
# 		`Cumulative_time_off` +
# 		`Year_of_hire` +
# 		`Cumulative_soluble_exposure` +
# 		`Cumulative_straight_exposure` +
# 		`Cumulative_synthetic_exposure` +
# 		`Employment_status`,
# 	data = cohort_analytic[,],
# 	family = binomial)
#
# summary(mod.gam)
#
# mod.glm <- glm(
# 	y ~ `Years_since_hire` +
# 		# `Calendar_year` +
# 		`Age` +
# 		Race + Plant + Sex +
# 		`Time_spent_in_assembly` +
# 		`Time_spent_machining` +
# 		`Time_spent_off` +
# 		`Cumulative_time_off` +
# 		`Year_of_hire` +
# 		`Cumulative_soluble_exposure` +
# 		`Cumulative_straight_exposure` +
# 		`Cumulative_synthetic_exposure` +
# 		`Employment_status`,
# 	data = cohort_analytic[,],
# 	family = binomial)
# summary(mod.glm)
#
# # Save models
# saveRDS(mod.gam, to_drive_D(here::here(
# 	paste0("left truncation/resources", "/FU from ", min(cohort_analytic$year),
# 				 " to ", max(cohort_analytic$year)),
# 	ifelse(length(probs) - 1 > 4, "mod.gam.rds", "small.mod.gam.rds"))))
# saveRDS(mod.glm, to_drive_D(here::here(
# 	paste0("left truncation/resources", "/FU from ", min(cohort_analytic$year),
# 				 " to ", max(cohort_analytic$year)),
# 	ifelse(length(probs) - 1 > 4, "mod.glm.rds", "small.mod.glm.rds"))))

# Get models
mod.gam <- readRDS(to_drive_D(here::here(
	paste0("left truncation/resources", "/FU from ", min(cohort_analytic$year),
				 " to ", max(cohort_analytic$year)),
	ifelse(length(probs) - 1 > 4, "mod.gam.rds", "small.mod.gam.rds"))))
mod.glm <- readRDS(to_drive_D(here::here(
	paste0("left truncation/resources", "/FU from ", min(cohort_analytic$year),
				 " to ", max(cohort_analytic$year)),
	ifelse(length(probs) - 1 > 4, "mod.glm.rds", "small.mod.glm.rds"))))

# Extract and save model output ####
source(here::here("left truncation", "get-coef.R"))
get.coef(
	model.path = to_drive_D(here::here(
		paste0("left truncation/resources", "/FU from ", min(cohort_analytic$year),
				 " to ", max(cohort_analytic$year)),
		ifelse(length(probs) - 1 > 4, "mod.gam.rds", "small.mod.gam.rds"))))
get.coef(
	model.path = to_drive_D(here::here(
		paste0("left truncation/resources", "/FU from ", min(cohort_analytic$year),
				 " to ", max(cohort_analytic$year)),
		ifelse(length(probs) - 1 > 4, "mod.glm.rds", "small.mod.glm.rds"))))

# Probability of being alive ####
# (1 - Pr(subject died at t = 1)) * (1 - Pr(subject died at t = 2)) * ... * (1 - Pr(subject died at t = t))
length(mod.gam$fitted.values); length(mod.glm$fitted.values); nrow(cohort_analytic[])
setorder(cohort_analytic, studyno, year)
cohort_analytic[, `:=`(
	gam.fitted = mod.gam$fitted.values,
	glm.fitted = mod.glm$fitted.values
)]
prob.tab <- cohort_analytic[, .(
	gam.prob = cumprod(1 - gam.fitted)[.N],
	glm.prob = cumprod(1 - glm.fitted)[.N],
	mortality = as.numeric(max(y) > 0),
	Age = Age[.N],
	`Year of hire` = Year_of_hire[1],
	Race = factor(race[1], levels = c("White", "Not white"), labels = c("White", "Black")),
	Sex = factor(sex[1], levels = c("M", "F"), labels = c("Men", "Women")),
	`Cumulative time off` = Cumulative_time_off[.N],
	`Cumulative straight` = Cumulative_straight_exposure[.N],
	`Cumulative soluble` = Cumulative_soluble_exposure[.N],
	`Cumulative synthetic` = Cumulative_synthetic_exposure[.N],
	`Employment status` = Employment_status[.N]
), by = .(studyno)]

# Probabilities by stratum ####
# By Age and employment status
saveRDS(
	prob.tab[mortality == 0, .(
		`Model 1` = mean(glm.prob),
		`Model 2` = mean(gam.prob)
	), by = .(Age)][order(Age),],
	here::here(paste0("reports/left truncation/survival to 1985/resources", "/FU from ", min(cohort_analytic$year),
				 " to ", max(cohort_analytic$year)),
						 paste0(ifelse(length(probs) > 5, "", "small."), "survival_by_age.rds"))
)

saveRDS(
	prob.tab[mortality == 0, .(
		`Model 1` = mean(glm.prob),
		`Model 2` = mean(gam.prob)
	), by = .(`Employment status`)][order(`Employment status`),],
	here::here(paste0("reports/left truncation/survival to 1985/resources", "/FU from ", min(cohort_analytic$year),
				 " to ", max(cohort_analytic$year)),
						 paste0(ifelse(length(probs) > 5, "", "small."), "survival_by_employment-status.rds"))
)
# By Race and Sex and Year of hire
saveRDS(
	prob.tab[mortality == 0, .(
		`Model 1` = mean(glm.prob),
		`Model 2` = mean(gam.prob)
	), by = .(Race)][order(Race)],
	here::here(paste0("reports/left truncation/survival to 1985/resources", "/FU from ", min(cohort_analytic$year),
				 " to ", max(cohort_analytic$year)),
						 paste0(ifelse(length(probs) > 5, "", "small."), "survival_by_race.rds"))
)

saveRDS(
	prob.tab[mortality == 0, .(
		`Model 1` = mean(glm.prob),
		`Model 2` = mean(gam.prob)
	), by = .(Sex)][order(Sex)],
	here::here(paste0("reports/left truncation/survival to 1985/resources", "/FU from ", min(cohort_analytic$year),
				 " to ", max(cohort_analytic$year)),
						 paste0(ifelse(length(probs) > 5, "", "small."), "survival_by_sex.rds"))
)

saveRDS(
	prob.tab[mortality == 0, .(
		`Model 1` = mean(glm.prob),
		`Model 2` = mean(gam.prob)
	), by = .(`Year of hire`)][order(`Year of hire`)],
	here::here(paste0("reports/left truncation/survival to 1985/resources", "/FU from ", min(cohort_analytic$year),
				 " to ", max(cohort_analytic$year)),
						 paste0(ifelse(length(probs) > 5, "", "small."), "survival_by_yin.rds"))
)

# By MWF exposure
saveRDS(
	prob.tab[mortality == 0, .(
		`Model 1` = mean(glm.prob),
		`Model 2` = mean(gam.prob)
	), by = .(`Cumulative straight`)][order(`Cumulative straight`)],
	here::here(paste0("reports/left truncation/survival to 1985/resources", "/FU from ", min(cohort_analytic$year),
				 " to ", max(cohort_analytic$year)),
						 paste0(ifelse(length(probs) > 5, "", "small."), "survival_by_straight.rds"))
)

saveRDS(
	prob.tab[mortality == 0, .(
		`Model 1` = mean(glm.prob),
		`Model 2` = mean(gam.prob)
	), by = .(`Cumulative soluble`)][order(`Cumulative soluble`)],
	here::here(paste0("reports/left truncation/survival to 1985/resources", "/FU from ", min(cohort_analytic$year),
				 " to ", max(cohort_analytic$year)),
						 paste0(ifelse(length(probs) > 5, "", "small."), "survival_by_soluble.rds"))
)

saveRDS(
	prob.tab[mortality == 0, .(
		`Model 1` = mean(glm.prob),
		`Model 2` = mean(gam.prob)
	), by = .(`Cumulative synthetic`)][order(`Cumulative synthetic`)],
	here::here(paste0("reports/left truncation/survival to 1985/resources", "/FU from ", min(cohort_analytic$year),
				 " to ", max(cohort_analytic$year)),
						 paste0(ifelse(length(probs) > 5, "", "small."), "survival_by_synthetic.rds"))
)


# ROC ####
library(pROC)
gam.roc <- roc(prob.tab[,.(mortality, p = 1 - gam.prob)],
							 mortality, p, ci = T)
glm.roc <- roc(prob.tab[,.(mortality, p = 1 - glm.prob)],
							 mortality, p, ci = T)

# # CIs for curve
# gam.sens <- ci.se(gam.roc, seq(0, 1, 0.025))
# glm.sens <- ci.se(glm.roc, seq(0, 1, 0.025))
#
# # Save CIs
# saveRDS(gam.sens,
# 				here::here(paste0("left truncation/resources", "/FU from ", min(cohort_analytic$year)),
# 									 paste0(ifelse(length(probs) > 5, "", "small."), "roc.gam.rds")))
# saveRDS(glm.sens,
# 				here::here(paste0("left truncation/resources", "/FU from ", min(cohort_analytic$year)),
# 									 paste0(ifelse(length(probs) > 5, "", "small."), "roc.gam.rds")))
#
# # Get CIs
# gam.sens <- readRDS(here::here(paste0("left truncation/resources", "/FU from ", min(cohort_analytic$year)),
# 															 paste0(ifelse(length(probs) > 5, "", "small."), "roc.gam.rds")))
# glm.sens <- readRDS(here::here(paste0("left truncation/resources", "/FU from ", min(cohort_analytic$year)),
# 															 paste0(ifelse(length(probs) > 5, "", "small."), "roc.gam.rds")))

# Plot
rbindlist(list(
	data.frame(
		Sensitivity = glm.roc$sensitivities,
		Specificity = glm.roc$specificities,
		AUC = as.character(round(as.numeric(glm.roc$auc), 3))),
	data.frame(
		Sensitivity = gam.roc$sensitivities,
		Specificity = gam.roc$specificities,
		AUC = as.character(round(as.numeric(gam.roc$auc), 3)))),
	idcol = "Model") %>% mutate(
		Model = factor(Model, levels = 1:2, labels = c("Model 1 (no splines)", "Model 2 (with splines)")),
		AUC = relevel(factor(AUC), as.character(round(as.numeric(glm.roc$auc), 3)))
	) -> roc.ggtab

setDT(roc.ggtab); setorder(roc.ggtab, Model, Specificity)
n_i <- 50
roc.ggtab[,`:=`(I = c(rep(1:n_i, length(Specificity) %/% n_i),
											seq(1, length.out = length(Specificity) - length(Specificity) %/% n_i * n_i))
), by = .(Model)]

roc.ggtab[I == 1] %>% ggplot(aes(
	x = Specificity, y = Sensitivity
)) + geom_step(aes(color = Model), size = 0.5) + geom_segment(
	x = -1, y = 0, xend = 0, yend = 1, color = "black", linetype = 2, size = 0.5) +
	geom_rect(aes(xmin = 0.5, xmax = 0.5, ymin = 0.5, ymax = 0.5, fill = AUC), alpha = 0) +
	guides(fill = guide_legend(
		override.aes = list(alpha = 1))) +
	scale_x_reverse() + mytheme -> roc.ggplot

roc.ggplot

# Compare two methods ####
gam.roc$auc; glm.roc$auc

# # Compile plot
# library(tikzDevice)
# tikz(file = here::here(paste0("reports/left truncation/survival to 1985/resources", "/FU from ", min(cohort_analytic$year),
# 				 " to ", max(cohort_analytic$year)),
# 											 ifelse(length(probs) - 1 > 4, "roc.tex", "small.roc.tex")),
# 		 height = 3, width = 4.75, standAlone = T)
# roc.ggplot
# dev.off()
#
# lualatex(pattern = ifelse(length(probs) - 1 > 4, "^roc\\.tex", "^small\\.roc\\.tex"),
# 				 directory = here::here(paste0("reports/left truncation/survival to 1985/resources", "/FU from ", min(cohort_analytic$year)),
# 				 " to ", max(cohort_analytic$year)),
# 				 break_after = 60)

# Make table with results from all models ####
prob.tab.og <- prob.tab
if (length(probs) > 5) {
	probs <- seq(0, 1, 0.25)
} else {probs <- seq(0, 1, 0.05)}


if ("Age" %in% names(cohort_analytic)) {
	cohort_analytic <- cohort_analytic[,-"Age", with = F]}
if ("Time_spent_machining" %in% names(cohort_analytic)) {
	cohort_analytic <- cohort_analytic[, -c(
		"Years_since_hire",
		"Calendar_year",
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
														 levels = c("At work", "Left work"))
)]

covariate.quantile <- apply(
	cohort_analytic[, .(
		sincehire.years, year, age,
		machining, off, cumulative_off,
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
	"Time_spent_machining",
	"Time_spent_off",
	"Cumulative_time_off",
	"Year_of_hire",
	"Cumulative_soluble_exposure",
	"Cumulative_straight_exposure",
	"Cumulative_synthetic_exposure"
)):=lapply(
	colnames(covariate.quantile), function(x) {
		cut(get(x), unique(covariate.quantile[,x]), include.lowest = T, dig.lab = 4)
	})]

# Get models
mod.gam <- readRDS(to_drive_D(here::here(
	paste0("left truncation/resources/FU from ", min(cohort_analytic$year),
				 " to ", max(cohort_analytic$year)),
	ifelse(length(probs) - 1 > 4, "mod.gam.rds", "small.mod.gam.rds"))))
mod.glm <- readRDS(to_drive_D(here::here(
	paste0("left truncation/resources/FU from ", min(cohort_analytic$year),
				 " to ", max(cohort_analytic$year)),
	ifelse(length(probs) - 1 > 4, "mod.glm.rds", "small.mod.glm.rds"))))

# Probability of being alive ####
# (1 - Pr(subject died at t = 1)) * (1 - Pr(subject died at t = 2)) * ... * (1 - Pr(subject died at t = t))
length(mod.gam$fitted.values); length(mod.glm$fitted.values); nrow(cohort_analytic[])
setorder(cohort_analytic, studyno, year)
cohort_analytic[, `:=`(
	tmp.gam.fitted = mod.gam$fitted.values,
	tmp.glm.fitted = mod.glm$fitted.values
)]
names(cohort_analytic)[grep("^tmp\\.g", names(cohort_analytic))] <-
	gsub("tmp", ifelse(length(probs) > 5, "big", "small"), grep("^tmp\\.g", names(cohort_analytic), value = T))
names(cohort_analytic)[grep("^g.*fitted", names(cohort_analytic))] <-
	paste0(ifelse(length(probs) > 5, "small.", "big."), grep("^g.*fitted", names(cohort_analytic), value = T))

prob.tab <- cohort_analytic[, .(
	gam.prob = cumprod(1 - gam.fitted)[.N],
	glm.prob = cumprod(1 - glm.fitted)[.N],
	mortality = as.numeric(max(y) > 0),
	Age = Age[.N],
	`Year of hire` = Year_of_hire[1],
	Race = factor(race[1], levels = c("White", "Not white"), labels = c("White", "Black")),
	Sex = factor(sex[1], levels = c("M", "F"), labels = c("Men", "Women")),
	`Cumulative time off` = Cumulative_time_off[.N],
	`Cumulative straight` = Cumulative_straight_exposure[.N],
	`Cumulative soluble` = Cumulative_soluble_exposure[.N],
	`Cumulative synthetic` = Cumulative_synthetic_exposure[.N],
	`Employment status` = Employment_status[.N]
), by = .(studyno)]

names(prob.tab)[grepl("gam.prob", names(prob.tab))] <- paste0(ifelse(length(probs) > 5, "big.", "small."), "gam.prob")
names(prob.tab)[grepl("glm.prob", names(prob.tab))] <- paste0(ifelse(length(probs) > 5, "big.", "small."), "glm.prob")

names(prob.tab.og)[grepl("gam.prob", names(prob.tab))] <- paste0(ifelse(length(probs) <= 5, "big.", "small."), "gam.prob")
names(prob.tab.og)[grepl("glm.prob", names(prob.tab))] <- paste0(ifelse(length(probs) <= 5, "big.", "small."), "glm.prob")

survival.tab <- merge(prob.tab[,grep("\\.prob|studyno", names(prob.tab)), with = F],
				prob.tab.og[,grep("\\.prob|dead|studyno", names(prob.tab)), with = F],
				on = "studyno")
box_write(
	merge(prob.tab[,grep("\\.prob|studyno", names(prob.tab)), with = F],
				prob.tab.og[,grep("\\.prob|dead|studyno", names(prob.tab)), with = F],
				on = "studyno"),
	file_name = paste0("survival_", min(cohort_analytic$year), "_to_", max(cohort_analytic$year), ".csv"),
	dir_id = 113431246688,
	description = paste0("Columns named `big.***.prob` were computed using the model with covariates of up to 20 levels. Columns named `small.***.prob` were computed using the model with covariates of no more than 4 levels. Probabilities computed with the model with splines are named `***.gam.prob`. Those with the model with categorically-coded variables only are named `***.glm.prob`. The column `mortality` indicates death due to natural causes in any year before ", max(cohort_analytic$year),"."))

# Save cumulative survival by year
by_year.tab <- cohort_analytic[, .(
	studyno,
	year,
	small.gam.prob = cumprod(1 - small.gam.fitted),
	small.glm.prob = cumprod(1 - small.glm.fitted),
	big.gam.prob = cumprod(1 - big.gam.fitted),
	big.glm.prob = cumprod(1 - big.glm.fitted),
	natural_cause_mortality = y,
	Age = Age
)]

# file id: 679752455177
box_write(by_year.tab,
	file_name = paste0("survival_", min(cohort_analytic$year), "_to_", max(cohort_analytic$year), "_by-year.csv"),
	dir_id = 113431246688,
	description = paste0("Columns named `big.***.prob` were computed using the model with covariates of up to 20 levels. Columns named `small.***.prob` were computed using the model with covariates of no more than 4 levels. Probabilities computed with the model with splines are named `***.gam.prob`. Those with the model with categorically-coded variables only are named `***.glm.prob`. The column `mortality` indicates death due to natural causes in any year before ", max(cohort_analytic$year),"."))
