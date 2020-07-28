# survival-sl.R
# Estimate survival in each year up to 1985 using SuperLearner
# Kevin Chen
# July 8, 2020

# Preliminaries
library(here)
library(SuperLearner)
# sl3, the Coyle, Nejazi, Malencia, Sofrygin implementation
# remotes::install_github("tlverse/sl3")
library(sl3)
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
# Variables to use
X.names <- c(
	"Duration_of_employment",
	"Calendar_year",
	"Age",
	"Race", "Plant", "Sex",
	"Cumulative_time_off",
	"Year_of_hire",
	"Cumulative_soluble_exposure",
	"Cumulative_straight_exposure",
	"Cumulative_synthetic_exposure",
	"Employment_status")
X_continuous.names <- c(
	"employment.years", "year",
	"age",
	"cumulative_off",
	"yin16",
	"cum_soluble",
	"cum_straight",
	"cum_synthetic",
	"Employment_status"
)
col.names <- unique(c(
	"studyno", "year", "y", "Age",
	X.names, X_continuous.names
))

set.seed(242)
studyno.sample <- sample(unique(cohort_analytic$studyno), 2000)
cohort_select <- cohort_analytic[
	# studyno %in% studyno.sample
	, col.names, with = F]

cohort_select$Calendar_year %>% is.factor

# Get covariate levels
sapply(X.names, function(x) {
	table(cohort_select[y == 1, x, with = F])
}, simplify = F)

# Get summaries
sapply(X_continuous.names, function(x) {
	summary(cohort_select[y == 1 & {
		if (grepl("cum", x)) {
			get(x) != 0
		} else {T}}, x, with = F])
}, simplify = F)

# Choose library
sl.library <- c("sl", "mean",
								"glm",
								# "glm.interaction",
								# "step",
								# "step.interaction",
								"glmnet",
								"gam",
								"ranger",
								# "knn", # breaks the cv
								# "rpart", # breaks the cv,
								"xgboost",
								NULL
)

# sl3 pipeline ####
# Make sure predict.gam returns "response"
Lrnr_gam$private_methods$.predict <- function (task) {
	predictions <- stats::predict(private$.fit_object, newdata = task$X, type = "response")
	predictions <- as.numeric(predictions)
	return(predictions)
}

# Instantiate task
task <- make_sl3_Task(as.data.frame(cohort_select),
											covariates = X.names,
											outcome = "y",
											outcome_type = "binomial")

# Make learners

# Instantiate
sapply(sl.library[-1], function(x = "rpart") {
	assign(paste0("lrnr_", x),
				 {if (sum(grepl(x, sl3_list_learners("binomial"))) > 0) {
				 	make_learner(get(paste0("Lrnr_", x)))
				 } else {
				 	# Use SuperLearner package, if needed
				 	make_learner(Lrnr_pkg_SuperLearner, paste0("SL.", x))}
				 },
				 envir = .GlobalEnv)
})

# # Run and check system time
# sapply(sl.library[-1], function(x) {
# 	system.time(
# 	assign(paste0("lrnr_", x, "_fit"),
# 				 get(paste0("lrnr_", x))$train(task),
# 				 envir = .GlobalEnv),
# 	F)
# })

# No screening for now, just stacking
stack <- make_learner(Stack,
											sapply(sl.library[-1], function(x) {
												get(paste0("lrnr_", x))
											}))
stack_fit <- stack$train(task)
stack_preds <- stack_fit$predict(task)
# Save stack results
saveRDS(stack_preds, file = to_drive_D(here::here(
	"resources",
	paste0("stack_preds", ifelse(length(table(sapply(probs, length))) == 1, "_by-deciles", ""), ".rds"))))

# Cross-validation
cv_stack <- Lrnr_cv$new(stack)
cv_fit <- cv_stack$train(task)
cv_preds <- cv_fit$predict(task)
# Save CV stack results
saveRDS(cv_preds, file = to_drive_D(here::here(
	"resources",
	paste0("cv_preds", ifelse(length(table(sapply(probs, length))) == 1, "_by-deciles", ""), ".rds"))))
cv_risks <- cv_fit$cv_risk(loss_squared_error)
saveRDS(cv_risks, file = to_drive_D(here::here(
	"resources",
	paste0("cv_risks", ifelse(length(table(sapply(probs, length))) == 1, "_by-deciles", ""), ".rds"))))

# Super Learner
metalearner <- make_learner(Lrnr_nnls)
cv_task <- cv_fit$chain(task)
ml_fit <- metalearner$train(cv_task)

sl_pipeline <- make_learner(Pipeline, stack_fit, ml_fit)
sl_weights <- sl_pipeline$learner_fits
# Save learner fits #
saveRDS(sl_weights, file = to_drive_D(here::here(
	"resources",
	paste0("sl_weights", ifelse(length(table(sapply(probs, length))) == 1, "_by-deciles", ""), ".rds"))))
sl_preds <- sl_pipeline$predict()
# Save sl results #
saveRDS(sl_preds, file = to_drive_D(here::here(
	"resources",
	paste0("sl_preds", ifelse(length(table(sapply(probs, length))) == 1, "_by-deciles", ""), ".rds"))))

# # Load sl3 results
# sl_preds <- readRDS(file = to_drive_D(here::here("resources", paste0("sl_preds", ifelse(length(table(sapply(probs, length))) == 1, "_by-deciles", ""), ".rds"))))
# sl_weights <- readRDS(file = to_drive_D(here::here("resources", paste0("sl_weights", ifelse(length(table(sapply(probs, length))) == 1, "_by-deciles", ""), ".rds"))))
# stack_preds <- readRDS(file = to_drive_D(here::here("resources", paste0("stack_preds", ifelse(length(table(sapply(probs, length))) == 1, "_by-deciles", ""), ".rds"))))
# cv_preds <- readRDS(file = to_drive_D(here::here("resources", paste0("cv_preds", ifelse(length(table(sapply(probs, length))) == 1, "_by-deciles", ""), ".rds"))))
# cv_risks <- readRDS(file = to_drive_D(here::here("resources", paste0("cv_risks", ifelse(length(table(sapply(probs, length))) == 1, "_by-deciles", ""), ".rds"))))

sl.predict <- data.table(
	sl = sl_preds,
	sapply(stack_preds, function(x) {x})
)

# # Basic implementation ####
# # Set up parallel computation - no windows
# num_cores <- RhpcBLASctl::get_num_cores()
# options(mc.cores = num_cores - 1)
# set.seed(236, "L'Ecuyer-CMRG")
# sl <- SuperLearner(
# 	Y = cohort_select$y,
# 	X = as.data.frame(cohort_select[,c(X.names), with = F]),
# 	family = "binomial",
# 	SL.library = paste0("SL.", sl.library[-1][-3]),
# 	cvControl = list(V = 10)
# 	# parallel = "multicore"
# 	)
# # # Save SuperLearner results ####
# # saveRDS(sl, file = to_drive_D(here::here("resources", "sl.rds")))
# # sl <- readRDS(file = to_drive_D(here::here("resources", "sl.rds")))
# 
# # Get predictions
# sl.predict <- as.data.table(
# 	cbind(studyno = cohort_select$studyno,
# 	year = cohort_select$year,
# 	sl$library.predict)
# 	)
# names(sl.predict) <- gsub("2_All$|^SL.", "", names(sl.predict))

# Save SL fitted values to Box
box_save(sl.predict,
				 dir_id = 117568282871,
				 file_name = paste0("sl.predict", ifelse(length(table(sapply(probs, length))) == 1, "_by-deciles", ""), ".rdata"),
				 description = "Predicted values from Super Learner")

# Probability of being alive ####
# (1 - Pr(subject died at t = 1)) * (1 - Pr(subject died at t = 2)) * ... * (1 - Pr(subject died at t = t))
# Fitted values
cohort_select[, (paste0(sl.library, ".fitted")):=(
	lapply(sl.predict, function(x) {x})
)]
# Cumulative products of the complements (survival)
cohort_select[, (paste0(sl.library, ".prob")):=(
	lapply(paste0(sl.library, ".fitted"), function(preds) {
		cumprod(1 - get(preds))})
), by = .(studyno)]
# Get ready for aggregating across year
cohort_select[,`:=`(
	I = 1:(.N),
	N = .N,
	mortality = as.numeric(max(y) > 0)),
	by = .(studyno)]

# Individual-specific
prob.tab <- cohort_select[I == N, c(
	"studyno", "year",
	paste0(sl.library, ".prob"),
	"mortality",
	"Age", "Year_of_hire", "Race", "Plant", "Sex",
	"Cumulative_time_off",
	"Cumulative_straight_exposure",
	"Cumulative_soluble_exposure",
	"Cumulative_synthetic_exposure",
	"Employment_status"
), with = F]
names(prob.tab) <- gsub(" exposure$", "", gsub("_", " ", names(prob.tab)))

# Probabilities by stratum ####
for (covariate in c("Age", "Employment status", "Race", "Sex", "Year of hire",
										paste("Cumulative", c("straight", "soluble", "synthetic")))) {
	tmp.tab <- prob.tab[mortality == 0, .(
		`sl` = mean(sl.prob),
		`mean` = mean(mean.prob),
		`glm` = mean(glm.prob),
		`glmnet` = mean(glmnet.prob),
		`gam` = mean(gam.prob),
		`ranger` = mean(ranger.prob),
		`xgboost` = mean(xgboost.prob)
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
		file = here::here("resources", paste0("survival_by-", gsub(" ", "-" , covariate), ifelse(length(table(sapply(probs, length))) == 1, "_by-deciles", ""), ".rds"))
	)
	assign(paste0(gsub(" ", "_" , tolower(covariate)), ".prob.tab"),
				 tmp.tab, envir = .GlobalEnv)
	print(tmp.tab)
}

# ROC ####
library(pROC)
for (method in sl.library) {
	assign(paste0(method, ".roc"),
				 roc(prob.tab[,.(mortality, p = 1 - get(paste0(method, ".prob")))],
				 		mortality, p, ci = T))
}

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
pair <- roc.ggtab[,.(AUC = paste0(AUC[1])), by = .(method)][order(AUC),]

roc.ggtab[,`:=`(
	method = factor(method, pair$method),
	AUC = factor(method, pair$method, pair$AUC)
)]

box_save(
	roc.ggtab,
	dir_id = 117612329554,
	file_name = paste0("roc", ifelse(length(table(sapply(probs, length))) == 1, "_by-deciles", ""), ".rdata"),
	description = "ROC results by method"
)

box_write(
	roc.ggtab,
	dir_id = 117612329554,
	file_name = paste0("roc", ifelse(length(table(sapply(probs, length))) == 1, "_by-deciles", ""), ".csv"),
	description = "ROC results by method"
)

roc.threshold <- roc.ggtab[, .(
	threshold.min = min(threshold[is.finite(threshold)]),
	threshold.Q2 = quantile(threshold[is.finite(threshold)], 0.25),
	threshold.med = median(threshold[is.finite(threshold)]),
	threshold.Q3 = quantile(threshold[is.finite(threshold)], 0.75),
	threshold.max = max(threshold[is.finite(threshold)])
), by = .(method)]

saveRDS(roc.threshold, file = here::here("resources", paste0("roc.threshold", ifelse(length(table(sapply(probs, length))) == 1, "_by-deciles", ""), ".rds")))

# Thin out number of lines for ggplot
n_i <- 40
roc.ggtab[,`:=`(I = c(
	rep(1:n_i, .N %/% n_i),
	seq(1, length.out = .N - .N %/% n_i * n_i))
), by = .(method)]

roc.ggtab[I == 1 | method == "mean"] %>% ggplot(aes(
	x = Specificity, y = Sensitivity, color = method
)) + geom_step(size = 0.5) +
	geom_segment(x = -1, y = 0, xend = 0, yend = 1, color = "black", linetype = 2, size = 0.5) +
	geom_rect(aes(xmin = 0.5, xmax = 0.5, ymin = 0.5, ymax = 0.5, fill = method), alpha = 0) +
	scale_fill_discrete(name = "AUC", labels = pair$AUC) +
	guides(fill = guide_legend(
		override.aes = list(alpha = 1))) +
	scale_x_reverse() + mytheme +
	theme(
		legend.box = "horizontal"
	) -> roc.ggplot

roc.ggplot

sort(sapply(sl.library, function(method) {
	get(paste0(method, ".roc"))$auc}))

# Compile plot
library(tikzDevice)
tikz(file = here::here("reports/survival to 1985/resources", paste0("sl.roc", ifelse(length(table(sapply(probs, length))) == 1, "_by-deciles", ""), ".tex")),
		 height = 3, width = 5, standAlone = T)
roc.ggplot
dev.off()

lualatex(pattern = paste0("^sl\\.roc", ifelse(length(table(sapply(probs, length))) == 1, "_by-deciles", ""), "\\.tex"),
				 directory = here::here("reports/survival to 1985/resources"),
				 break_after = 120)

# # Look at the weights ####
# weight.ggtab <- melt(prob.tab,
# 		 id.vars = 1,
# 		 measure = 2:5)[,.(
# 		 	studyno,
# 		 	method = substr(variable, 1,  unlist(gregexpr("\\.", variable)) - 1),
# 		 	p = value
# 		 )]
# 
# # Truncate?
# weight.ggtab[,`:=`(
# 	p = ifelse(p > quantile(p, 0.99), quantile(p, 0.99), p)
# ), by = .(method)]
# 
# weight.ggtab[,`:=`(
# 	p.mean = mean(p)
# ), by = .(method)]
# 
# weight.ggtab[,`:=`(method = relevel(factor(method), "sl"))]
# 
# weight.ggtab %>% ggplot(aes(
# 	x = (1 - p.mean)/(1 - p)
# )) + geom_histogram(
# 	# aes(y = ..density..),
# 	bins = 500) +
# 	# geom_density(color = "red") +
# 	# coord_cartesian(xlim = c(0, 100)) +
# 	facet_wrap(. ~ method) + mytheme.web
