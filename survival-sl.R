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

# Fix theme elements
mytheme <- mytheme +
	theme(strip.text = element_text(size = 7, margin = margin(1.5, 1.5, 1.5, 1.5, "pt"))
				)

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
	X.names, paste0(X.names[X.names %in% new_cat.cols], "_deciles"),
	X_continuous.names
))

set.seed(242)
studyno.sample <- sample(unique(cohort_analytic$studyno), 1500)
cohort_select <- cohort_analytic[
	# studyno %in% studyno.sample
	, col.names, with = F]

# Get covariate levels
sapply(X.names, function(x) {
	table(cohort_select[y == 1, x, with = F])
}, simplify = F)
# Covariate levels (deciles)
sapply(X.names, function(x) {
	if (x %in% new_cat.cols) {
		x <- paste0(x, "_deciles")
	}
	table(cohort_select[y == 1, x, with = F])
}, simplify = F)

# Get summaries for continuous variables
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
run_sl3 <- F
# Make sure predict.gam returns "response"
Lrnr_gam$private_methods$.predict <- function (task) {
	predictions <- stats::predict(private$.fit_object, newdata = task$X, type = "response")
	predictions <- as.numeric(predictions)
	return(predictions)
}

# Run sl3 for both sets of cutpoints ####
for (i in length(probs)) { # loop over the different sets of cutpoints
	is.deciles <- length(table(sapply(c(probs[[i]], deciles = list(seq(0, 1, 0.1))), length))) == 1

	# Choose covariates
	X <- X.names
	if (is.deciles) {
		X <- sapply(X, function(x) {
			if (x %in% new_cat.cols) {paste0(x, "_deciles")} else {x}}, USE.NAMES = F)
	}

	if (run_sl3) {
		# Instantiate task
		task <- make_sl3_Task(as.data.frame(cohort_select),
													covariates = X,
													outcome = "y",
													outcome_type = "binomial")

		# Instantiate learners
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
		saveRDS(stack_preds, file = here::here(
			"resources/sl", paste0("stack_preds", ifelse(is.deciles, "_by-deciles", ""), ".rds")))

		# Cross-validation
		cv_stack <- Lrnr_cv$new(stack)
		cv_fit <- cv_stack$train(task)
		cv_preds <- cv_fit$predict(task)
		# Save CV stack results
		saveRDS(cv_preds, file = here::here(
			"resources/sl", paste0("cv_preds", ifelse(is.deciles, "_by-deciles", ""), ".rds")))
		cv_risks <- cv_fit$cv_risk(loss_squared_error)
		saveRDS(cv_risks, file = here::here(
			"resources/sl", paste0("cv_risks", ifelse(is.deciles, "_by-deciles", ""), ".rds")))

		# Super Learner
		metalearner <- make_learner(Lrnr_nnls)
		cv_task <- cv_fit$chain(task)
		ml_fit <- metalearner$train(cv_task)

		sl_pipeline <- make_learner(Pipeline, stack_fit, ml_fit)
		sl_weights <- sl_pipeline$learner_fits
		# Save learner fits #
		saveRDS(sl_weights, file = to_drive_D(here::here(
			"resources/sl", paste0("sl_weights", ifelse(is.deciles, "_by-deciles", ""), ".rds"))))
		sl_preds <- sl_pipeline$predict()
		# Save sl results #
		saveRDS(sl_preds, file = here::here(
			"resources/sl", paste0("sl_preds", ifelse(is.deciles, "_by-deciles", ""), ".rds")))

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
		# 	X = as.data.frame(cohort_select[,c(X), with = F]),
		# 	family = "binomial",
		# 	SL.library = paste0("SL.", sl.library[-1][-3]),
		# 	cvControl = list(V = 10)
		# 	# parallel = "multicore"
		# 	)
		# # # Save SuperLearner results ####
		# # saveRDS(sl, file = to_drive_D(here::here("resources/sl", "sl.rds")))
		# # sl <- readRDS(file = to_drive_D(here::here("resources/sl", "sl.rds")))
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
						 file_name = paste0("sl.predict", ifelse(is.deciles, "_by-deciles", ""), ".rdata"),
						 description = "Predicted values from Super Learner")
	} else {
		# Load sl3 results (rather than run)
		sl_preds <- readRDS(file = here::here("resources/sl", paste0("sl_preds", ifelse(is.deciles, "_by-deciles", ""), ".rds")))
		sl_weights <- readRDS(file = to_drive_D(here::here("resources/sl", paste0("sl_weights", ifelse(is.deciles, "_by-deciles", ""), ".rds"))))
		stack_preds <- readRDS(file = here::here("resources/sl", paste0("stack_preds", ifelse(is.deciles, "_by-deciles", ""), ".rds")))
		cv_preds <- readRDS(file = here::here("resources/sl", paste0("cv_preds", ifelse(is.deciles, "_by-deciles", ""), ".rds")))
		cv_risks <- readRDS(file = here::here("resources/sl", paste0("cv_risks", ifelse(is.deciles, "_by-deciles", ""), ".rds")))
		sl.predict <- box_read(ifelse(is.deciles, 696590655167, 691257036647))
	}

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
	prob.tab <- cohort_select[I == N, sapply(c(
		"studyno", "year",
		paste0(sl.library, ".prob"),
		"mortality",
		"Age", "Year_of_hire", "Race", "Plant", "Sex",
		"Cumulative_time_off",
		"Cumulative_straight_exposure",
		"Cumulative_soluble_exposure",
		"Cumulative_synthetic_exposure",
		"Employment_status"
	), function(x) {
		if (x %in% new_cat.cols & is.deciles) {
		paste0(x, "_deciles")} else {x}
		}), with = F]
	names(prob.tab) <- gsub(" exposure| deciles", "", gsub("_", " ", names(prob.tab)))

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
		saveRDS(
			tmp.tab,
			file = here::here("resources/survival by covariate", paste0(gsub(" ", "-" , covariate), ifelse(is.deciles, "_by-deciles", ""), ".rds"))
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
		file_name = paste0("roc", ifelse(is.deciles, "_by-deciles", ""), ".rdata"),
		description = "ROC results by method"
	)

	box_write(
		roc.ggtab,
		dir_id = 117612329554,
		file_name = paste0("roc", ifelse(is.deciles, "_by-deciles", ""), ".csv"),
		description = "ROC results by method"
	)

	roc.threshold <- roc.ggtab[, .(
		threshold.min = min(threshold[is.finite(threshold)]),
		threshold.Q2 = quantile(threshold[is.finite(threshold)], 0.25),
		threshold.med = median(threshold[is.finite(threshold)]),
		threshold.Q3 = quantile(threshold[is.finite(threshold)], 0.75),
		threshold.max = max(threshold[is.finite(threshold)])
	), by = .(method)]

	saveRDS(roc.threshold, file = here::here("resources/roc", paste0("roc.threshold", ifelse(is.deciles, "_by-deciles", ""), ".rds")))

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
		guides(
			color = guide_legend(order = 1),
			fill = guide_legend(
				override.aes = list(alpha = 1),
				order = 2)) +
		scale_x_reverse() + mytheme +
		theme(
			legend.box = "horizontal"
		) -> roc.ggplot


	roc.sort <- sort(sapply(sl.library, function(method) {
		get(paste0(method, ".roc"))$auc}))
	roc.sort

	# Compile plot
	library(tikzDevice)
	# quartz(height = 3, width = 5)
	tikz(file = here::here("reports/survival to 1985/resources", paste0("sl.roc", ifelse(is.deciles, "_by-deciles", ""), ".tex")),
			 height = 3, width = 5, standAlone = T)
	roc.ggplot
	dev.off()

	lualatex(pattern = paste0("^sl\\.roc", ifelse(is.deciles, "_by-deciles", ""), "\\.tex"),
					 directory = here::here("reports/survival to 1985/resources"),
					 break_after = 120)

	# Look at the weights ####
	weight.ggtab <- melt(prob.tab,
											 id.vars = 1,
											 measure = which(grepl("\\.prob", names(prob.tab))))[,.(
											 	studyno,
											 	method = substr(variable, 1,  unlist(gregexpr("\\.", variable)) - 1),
											 	p = value
											 )]

	weight.ggtab[,`:=`(method = factor(factor(method), levels = names(roc.sort)))]

	weight.ggtab %>% group_by(method) %>%
		mutate(p.mean = mean(p)) %>%
		mutate(weight = (1 - p.mean)/(1 - p)) %>%
		mutate(
			weight = ifelse(weight > quantile(weight, 0.95), quantile(weight, 0.95), weight)
			) %>% ggplot(aes(
		x = weight
	)) + geom_histogram(
		# aes(y = ..density..),
		bins = 500) +
		# geom_density(color = "red") +
		# coord_cartesian(xlim = c(0, 100)) +
		facet_wrap(. ~ method, scales = "free") + mytheme -> weight.ggplot

	# quartz(height = 4, width = 5)
	tikz(file = here::here("reports/survival to 1985/resources", paste0("weights", ifelse(is.deciles, "_by-deciles", ""), ".tex")),
			 height = 4, width = 5, standAlone = T)
	weight.ggplot
	dev.off()

	lualatex(pattern = paste0("weights", ifelse(is.deciles, "_by-deciles", ""), "\\.tex"),
					 directory = here::here("reports/survival to 1985/resources"),
					 break_after = 120)

	weight.ggtab %>% group_by(method) %>% mutate(
		p = ifelse(p > quantile(p, 0.95), quantile(p, 0.95), p)
	) %>% ggplot(aes(
		x = p
	)) + geom_histogram(
		# aes(y = ..density..),
		bins = 500) +
		# geom_density(color = "red") +
		# coord_cartesian(xlim = c(0, 100)) +
		facet_wrap(. ~ method, scales = "free") + mytheme -> prob.ggplot

	# quartz(height = 4, width = 5)
	tikz(file = here::here("reports/survival to 1985/resources", paste0("probs", ifelse(is.deciles, "_by-deciles", ""), ".tex")),
			 height = 4, width = 5, standAlone = T)
	prob.ggplot
	dev.off()

	lualatex(pattern = paste0("probs", ifelse(is.deciles, "_by-deciles", ""), "\\.tex"),
					 directory = here::here("reports/survival to 1985/resources"),
					 break_after = 120)

} # End loop over different cutpoints