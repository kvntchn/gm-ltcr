# survival-sl.R
# Estimate survival in each year up to 1985 using SuperLearner
# Kevin Chen
# July 8, 2020

# Preliminaries
library(here)
library(SuperLearner)
# # XGBoost
# install.packages("drat", repos="https://cran.rstudio.com")
# drat:::addRepo("dmlc")
# install.packages("xgboost", repos="http://dmlc.ml/drat/", type = "source") # not available for R 4.0.0

# Get data and other dependencies ####
if (!('cohort' %in% ls(envir = .GlobalEnv))) {
source(here::here('get-data.R'))
}
setorder(cohort_analytic, studyno, year)

# Drop unnecessary variables and take subset for fast testing
X_cat.names <- c(
	"Years_since_hire", "Age",
	"Race", "Plant", "Sex", "Year_of_hire",
	"Time_spent_machining", "Time_spent_off", "Cumulative_time_off",
	"Cumulative_soluble_exposure",
	"Cumulative_straight_exposure",
	"Cumulative_synthetic_exposure",
	"Employment_status")
X.names <- c(
	"sincehire.years", "age",
	"yin16",
	"machining", "off", "cumulative_off",
	"cum_soluble",
	"cum_straight",
	"cum_synthetic"
)
col.names <- c(
	"studyno", "year",
	"age.year1", "age.year2", "year1", "year2",
	"y",
	X_cat.names, X.names
)

set.seed(236)
studyno.sample <- sample(unique(cohort_analytic$studyno), 1500)
cohort_select <- cohort_analytic[
	studyno %in% studyno.sample
	, col.names, with = F]

# Check data characteristics
str(cohort_select)

# # One-hot encoding? ####
# col.is_long_factor <- sapply(cohort_select, function(x) {is.factor(x) & length(levels(x)) > 2})
# one_hot <- sapply(which(col.is_long_factor), function(i) {
# 	one_hot <- as.data.table(model.matrix(~ . - 1, cohort_select[, i, with = F]))
# 	return(one_hot)
# })
# 
# # cbind it up!
# one_hot.bind <- function(x) {
# 	x[[2]] <- cbind(x[[1]], x[[2]])
# 	x <- x[-1]
# 	if (length(x) > 1) {
# 		one_hot.bind(x)
# 	} else {
# 		return(x[[1]])
# 	}}
# 
# cohort_select <- cbind(cohort_select[,!col.is_long_factor, with = F], one_hot.bind(one_hot))

# Fit ####
sl.library <- c("sl", "mean",
								"lm", "glm", # categorical
								"glm.interaction", # continuous
								"step", "step.interaction", # categorical
								"glmnet", "biglasso", # categorical
								# "ridge", # Not possible for binomial
								"gam", # continuous
								# "loess", # Not possible for binomial
								# "earth", # continuous # very long computation
								"ranger", # categorical
								"knn", # continuous, else too many ties
								# "ksvm", # continous, else vector memory exausted
								# "svm", # Infeasible
								"rpart", "rpartPrune", # categorical
								# "xgboost" # categorical # very long computation time
								)
source(here::here("sl_library.R"))
# Set up parallel computation - no windows
num_cores <- RhpcBLASctl::get_num_cores()
options(mc.cores = num_cores - 1)
set.seed(236, "L'Ecuyer-CMRG")
sl <- CV.SuperLearner(
	Y = cohort_select$y,
	X = as.data.frame(cohort_select[,c(X_cat.names, X.names), with = F]),
	family = "binomial",
	SL.library = paste0("SL.", sl.library[-1], "2"),
	cvControl = list(V = 10),
	parallel = "multicore"
	)

saveRDS(sl, file = to_drive_D(here::here("resources", "sl.rdata")))
# sl <- readRDS(file = to_drive_D(here::here("resources", "sl.rdata")))

# Get predictions
sl.predict <- predict(sl, as.data.frame(cohort_select[,c(X_cat.names, X.names), with = F]), onlySL = T)

sl.predict <- data.table(
	studyno = cohort_select$studyno,
	year = cohort_select$year,
	sl = as.vector(sl.predict$pred),
	sl.predict$library.predict)

names(sl.predict)[-(1:3)] <- sl.library[-1]

# box_save(sl.predict,
# 				 dir_id = 117568282871,
# 				 file_name = paste0("sl.predict.rdata"),
# 				 description = "Predicted values from Super Learner")

# Probability of being alive ####
# (1 - Pr(subject died at t = 1)) * (1 - Pr(subject died at t = 2)) * ... * (1 - Pr(subject died at t = t))
cohort_select[, (paste0(sl.library, ".fitted")):=(
	lapply(sl.predict[,-(1:2)], function(x) {x})
)]
cohort_select[, (paste0(sl.library, ".prob")):=(
	lapply(sl.predict[,-(1:2)], function(x) {
		cumprod(1 - x)})
)]
cohort_select[,`:=`(
	I = 1:(.N),
	N = .N,
	mortality = as.numeric(max(y) > 0)),
	by = .(studyno)]
prob.tab <- cohort_select[I == N, c(
	"studyno",
	paste0(sl.library, ".prob"),
	"mortality",
	"Age", "Year_of_hire", "Race", "Sex",
	"Cumulative_time_off",
	"Cumulative_straight_exposure",
	"Cumulative_soluble_exposure",
	"Cumulative_synthetic_exposure",
	"Employment_status"
), with = F]
names(prob.tab) <- gsub(" exposure$", "", gsub("_", " ", names(prob.tab)))

# Probabilities by stratum ####
for (covariate in c("Age", "Employment status", "Race", "Sex", "Year of hire")) {
	tmp.tab <- prob.tab[mortality == 0, .(
		`sl` = mean(sl.prob),
		`mean` = mean(mean.prob),
		`glm` = mean(glm.prob),
		`ranger` = mean(ranger.prob)
	), by = .(get(covariate))][order(get),]
	names(tmp.tab)[1] <- covariate
	# box_save(
	# 	tmp.tab,
	# 	dir_id = 117612329554,
	# 	file_name = paste0("survival_by-", gsub(" ", "-" , covariate), ".rdata"),
	# 	description = paste("Average survival probability stratified by", tolower(covariate))
	# )
	saveRDS(
		tmp.tab,
		file = here::here("resources", paste0("survival_by-", gsub(" ", "-" , covariate), ".rds"))
	)
	assign(paste0(gsub(" ", "_" , tolower(covariate)), ".prob.tab"),
				 tmp.tab, envir = .GlobalEnv)
}

# ROC ####
library(pROC)
for (method in sl.library) {
	assign(paste0(method, ".roc"),
				 roc(prob.tab[,.(mortality, p = 1 - get(paste0(method, ".prob")))],
				 		mortality, p, ci = T))}

# Plot
roc.ggtab <- rbindlist(lapply(sl.library, function (method) {
	data.frame(
		Sensitivity = get(paste0(method, ".roc"))$sensitivities,
		Specificity = get(paste0(method, ".roc"))$specificities,
		AUC = as.character(round(as.numeric(get(paste0(method, ".roc"))$auc), 3)),
		threshold = get(paste0(method, ".roc"))$thresholds,
		method = paste(method))
}))

setDT(roc.ggtab)

roc.ggtab[,`:=`(AUC = factor(AUC, sort(as.numeric(unique(AUC)), T)))]
roc.ggtab[,`:=`(method = factor(
	method, roc.ggtab[,.(method = method[1]), by = .(AUC)][order(AUC),]$method
	))]

box_save(
	roc.ggtab,
	dir_id = 117612329554,
	file_name = "roc.rdata",
	description = "ROC results by method"
)

box_write(
	roc.ggtab,
	dir_id = 117612329554,
	file_name = "roc.csv",
	description = "ROC results by method"
)

roc.threshold <- roc.ggtab[,.(
	threshold.min = min(threshold[is.finite(threshold)]),
	threshold.Q2 = quantile(threshold[is.finite(threshold)], 0.25),
	threshold.med = median(threshold[is.finite(threshold)]),
	threshold.Q3 = quantile(threshold[is.finite(threshold)], 0.75),
	threshold.max = max(threshold[is.finite(threshold)])
	), by = .(method)]

saveRDS(roc.threshold,
				file = here::here("resources",
													"roc.threshold.rds"))

# Thin out number of lines for ggplot
n_i <- 40
roc.ggtab[,`:=`(I = c(rep(1:n_i, length(Specificity) %/% n_i),
											seq(1, length.out = length(Specificity) - length(Specificity) %/% n_i * n_i))
)]

roc.ggtab[I == 1 | method == "mean"] %>% ggplot(aes(
	x = Specificity, y = Sensitivity, color = method
)) + geom_step(size = 0.5) +
	geom_segment(x = -1, y = 0, xend = 0, yend = 1, color = "black", linetype = 2, size = 0.5) +
	geom_rect(aes(xmin = 0.5, xmax = 0.5, ymin = 0.5, ymax = 0.5,
								fill = AUC), alpha = 0) +
	guides(fill = guide_legend(
		override.aes = list(alpha = 1))) +
	scale_x_reverse() + mytheme -> roc.ggplot

roc.ggplot

sapply(sl.library, function(method) {
	get(paste0(method, ".roc"))$auc})

# # Compile plot
# library(tikzDevice)
# tikz(file = here::here("reports/survival to 1985/resources", "sl.roc.tex"),
# 		 height = 3, width = 4.75, standAlone = T)
# roc.ggplot
# dev.off()
# 
# lualatex(pattern = "^sl\\.roc\\.tex",
# 				 directory = here::here("reports/survival to 1985/resources"),
# 				 break_after = 60)


# Look at the weights ####
weight.ggtab <- melt(prob.tab,
		 id.vars = 1,
		 measure = 2:5)[,.(
		 	studyno,
		 	method = substr(variable, 1,  unlist(gregexpr("\\.", variable)) - 1),
		 	p = value
		 )]

# Truncate?
weight.ggtab[,`:=`(
	p = ifelse(p > quantile(p, 0.99), quantile(p, 0.99), p)
), by = .(method)]

weight.ggtab[,`:=`(
	p.mean = mean(p)
), by = .(method)]

weight.ggtab[,`:=`(method = relevel(factor(method), "sl"))]

weight.ggtab %>% ggplot(aes(
	x = (1 - p.mean)/(1 - p)
)) + geom_histogram(
	# aes(y = ..density..),
	bins = 500) +
	# geom_density(color = "red") +
	# coord_cartesian(xlim = c(0, 100)) +
	facet_wrap(. ~ method) + mytheme.web
