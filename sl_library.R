# SL library to use different variables for X
# Kevin Chen
# July 14

# # Print available functions to console
# invisible(sapply(sl.library[-1], function(x) {
#     cat(paste0("SL.", x, "2 <-", "\n"))
#     print(get(paste0("SL.", x)))
#     if (paste0("predict.SL.", x) %in% ls("package:SuperLearner")) {
#     	cat(paste0("predict.SL.", x, "2 <-", "\n"))
#     print(get(paste0("predict.SL.", x)))}
# }))

SL.mean2 <- function (Y, X, newX, family, obsWeights, id, ...) {
	meanY <- weighted.mean(Y, w = obsWeights)
	pred <- rep.int(meanY, times = nrow(newX))
	fit <- list(object = meanY)
	out <- list(pred = pred, fit = fit)
	class(out$fit) <- c("SL.mean")
	return(out)
}


predict.SL.mean2 <- function (object, newdata, family, X = NULL, Y = NULL, ...) {
	pred <- rep.int(object$object, times = nrow(newdata))
	return(pred)
}


SL.lm2 <- function (Y, X, newX, family, obsWeights, model = TRUE, ...) {
	X = X[,which(names(X) %in% X_cat.names)]
	
	if (is.matrix(X)) {
		X = as.data.frame(X)
	}
	fit <- stats::lm(Y ~ ., data = X, weights = obsWeights, model = model)
	if (is.matrix(newX)) {
		newX = as.data.frame(newX)
	}
	pred <- predict(fit, newdata = newX, type = "response")
	if (family$family == "binomial") {
		pred = pmin(pmax(pred, 0), 1)
	}
	fit <- list(object = fit, family = family)
	class(fit) <- "SL.lm"
	out <- list(pred = pred, fit = fit)
	return(out)
}


SL.glm2 <- function (Y, X, newX, family, obsWeights, model = TRUE, ...) {
	X = X[,which(names(X) %in% X_cat.names)]
	
	if (is.matrix(X)) {
		X = as.data.frame(X)
	}
	fit.glm <- glm(Y ~ ., data = X, family = family, weights = obsWeights, 
								 model = model)
	if (is.matrix(newX)) {
		newX = as.data.frame(newX)
	}
	pred <- predict(fit.glm, newdata = newX, type = "response")
	fit <- list(object = fit.glm)
	class(fit) <- "SL.glm"
	out <- list(pred = pred, fit = fit)
	return(out)
}


SL.glm.interaction2 <- function (Y, X, newX, family, obsWeights, ...) {
	X = X[,which(names(X) %in% X.names)]
	
	if (is.matrix(X)) {
		X = as.data.frame(X)
	}
	fit.glm <- glm(Y ~ .^2, data = X, family = family, weights = obsWeights)
	if (is.matrix(newX)) {
		newX = as.data.frame(newX)
	}
	pred <- predict(fit.glm, newdata = newX, type = "response")
	fit <- list(object = fit.glm)
	class(fit) <- "SL.glm"
	out <- list(pred = pred, fit = fit)
	return(out)
}


SL.step2 <- function (Y, X, newX, family, direction = "both", trace = 0, 
											k = 2, ...) {
	X = X[,which(names(X) %in% X_cat.names)]
	
	fit.glm <- glm(Y ~ ., data = X, family = family)
	fit.step <- step(fit.glm, direction = direction, trace = trace, 
									 k = k)
	pred <- predict(fit.step, newdata = newX, type = "response")
	fit <- list(object = fit.step)
	out <- list(pred = pred, fit = fit)
	class(out$fit) <- c("SL.step")
	return(out)
}

SL.step.interaction2 <- function (Y, X, newX, family, direction = "both", trace = 0, 
																	k = 2, ...) {
	X = X[,which(names(X) %in% X.names)]
	
	fit.glm <- glm(Y ~ ., data = X, family = family)
	fit.step <- step(fit.glm, scope = Y ~ .^2, direction = direction, 
									 trace = trace, k = k)
	pred <- predict(fit.step, newdata = newX, type = "response")
	fit <- list(object = fit.step)
	out <- list(pred = pred, fit = fit)
	class(out$fit) <- c("SL.step")
	return(out)
}


SL.glmnet2 <- function (Y, X, newX, family = "binomial", obsWeights, id, alpha = 1, nfolds = 10, 
												nlambda = 100, useMin = TRUE, loss = "deviance", ...) {
	X = X[,which(names(X) %in% X_cat.names)]
	
	require("glmnet")
	if (!is.matrix(X)) {
		X <- model.matrix(~-1 + ., X)
		newX <- model.matrix(~-1 + ., newX[,which(names(newX) %in% X_cat.names)])
	}
	fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, 
														 lambda = NULL, type.measure = loss, nfolds = nfolds, 
														 family = family$family, alpha = alpha, nlambda = nlambda, 
														 ...)
	pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin, 
																																		"lambda.min", "lambda.1se"))
	fit <- list(object = fitCV, useMin = useMin)
	class(fit) <- "SL.glmnet"
	out <- list(pred = pred, fit = fit)
	return(out)
}

predict.SL.glmnet <- function (object, newdata, remove_extra_cols = T, add_missing_cols = T, 
															 ...) {
	require("glmnet")
	if (!is.matrix(newdata)) {
		newdata <- model.matrix(~-1 + ., newdata[,which(names(newdata) %in% X_cat.names)])
	}
	original_cols = rownames(object$object$glmnet.fit$beta)
	if (remove_extra_cols) {
		extra_cols = setdiff(colnames(newdata), original_cols)
		if (length(extra_cols) > 0) {
			warning(paste("Removing extra columns in prediction data:", 
										paste(extra_cols, collapse = ", ")))
			newdata = newdata[, !colnames(newdata) %in% extra_cols, 
												drop = FALSE]
		}
	}
	if (add_missing_cols) {
		missing_cols = setdiff(original_cols, colnames(newdata))
		if (length(missing_cols) > 0) {
			warning(paste("Adding missing columns in prediction data:", 
										paste(missing_cols, collapse = ", ")))
			new_cols = matrix(0, nrow = nrow(newdata), ncol = length(missing_cols))
			colnames(new_cols) = missing_cols
			newdata = cbind(newdata, new_cols)
			newdata = newdata[, original_cols, drop = FALSE]
		}
	}
	pred <- predict(object$object, newx = newdata, type = "response", 
									s = ifelse(object$useMin, "lambda.min", "lambda.1se"))
	return(pred)
}

SL.biglasso2 <- function (Y, X, newX, family, obsWeights, penalty = "lasso", 
													alg.logistic = "Newton", screen = "SSR", alpha = 1, nlambda = 100, 
													eval.metric = "default", ncores = 1, nfolds = 5, ...) {
	X = X[,which(names(X) %in% X_cat.names)]
	
	require("biglasso")
	require("bigmemory")
	if (!is.matrix(X)) {
		X = model.matrix(~., X)
		X = X[, -1]
	}
	X = bigmemory::as.big.matrix(X)
	fit = biglasso::cv.biglasso(X, Y, family = family$family, 
															penalty = penalty, alg.logistic = alg.logistic, screen = screen, 
															eval.metric = eval.metric, ncores = ncores, alpha = alpha, 
															nfolds = nfolds, nlambda = nlambda)
	if (!is.matrix(newX)) {
		newX = model.matrix(~., newX)
		newX = newX[, -1]
	}
	newX = bigmemory::as.big.matrix(newX)
	pred <- predict(fit, newX, type = "response")
	fit <- list(object = fit)
	class(fit) <- c("SL.biglasso")
	out <- list(pred = as.vector(pred), fit = fit)
	return(out)
}


SL.gam2 <- function (Y, X, newX, family, obsWeights, deg.gam = 2, cts.num = 4, 
										 ...) {
	X = X[,which(names(X) %in% X.names)]
	
	if (!require("gam")) {
		stop("SL.gam requires the gam package, but it isn't available")
	}
	if ("mgcv" %in% loadedNamespaces()) 
		warning("mgcv and gam packages are both in use. You might see an error because both packages use the same function names.")
	cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
	if (sum(!cts.x) > 0) {
		gam.model <- as.formula(paste("Y~", paste(paste("s(", 
																										colnames(X[, cts.x, drop = FALSE]), ",", deg.gam, 
																										")", sep = ""), collapse = "+"), "+", paste(colnames(X[, 
																																																					 !cts.x, drop = FALSE]), collapse = "+")))
	}
	else {
		gam.model <- as.formula(paste("Y~", paste(paste("s(", 
																										colnames(X[, cts.x, drop = FALSE]), ",", deg.gam, 
																										")", sep = ""), collapse = "+")))
	}
	if (sum(!cts.x) == length(cts.x)) {
		gam.model <- as.formula(paste("Y~", paste(colnames(X), 
																							collapse = "+"), sep = ""))
	}
	fit.gam <- gam::gam(gam.model, data = X, family = family, 
											control = gam::gam.control(maxit = 50, bf.maxit = 50), 
											weights = obsWeights)
	if (packageVersion("gam") >= 1.15) {
		pred <- gam::predict.Gam(fit.gam, newdata = newX, type = "response")
	}
	else {
		stop("This SL.gam wrapper requires gam version >= 1.15, please update the gam package with 'update.packages('gam')'")
	}
	fit <- list(object = fit.gam)
	out <- list(pred = pred, fit = fit)
	class(out$fit) <- c("SL.gam")
	return(out)
}


SL.loess2 <- function (Y, X, newX, family, obsWeights, span = 0.75, l.family = "gaussian", 
											 ...) {
	X = X[,which(names(X) %in% X.names)]
	
	if (family$family == "gaussian") {
		fit.loess <- loess(as.formula(paste("Y~", names(X))), 
											 data = X, family = l.family, span = span, control = loess.control(surface = "direct"), 
											 weights = obsWeights)
	}
	if (family$family == "binomial") {
		stop("family = binomial() not currently implemented for SL.loess")
	}
	pred <- predict(fit.loess, newdata = newX)
	fit <- list(object = fit.loess)
	out <- list(pred = pred, fit = fit)
	class(out$fit) <- c("SL.loess")
	return(out)
}


SL.earth2 <- function (Y, X, newX, family, obsWeights, id, degree = 2, penalty = 3, 
											 nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", nfold = 0, 
											 ncross = 1, minspan = 0, endspan = 0, ...) {
	X = X[,which(names(X) %in% X.names)]
	
	require("earth")
	if (family$family == "gaussian") {
		fit.earth <- earth::earth(x = X, y = Y, degree = degree, 
															nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
															ncross = ncross, minspan = minspan, endspan = endspan)
	}
	if (family$family == "binomial") {
		fit.earth <- earth::earth(x = X, y = Y, degree = degree, 
															nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
															ncross = ncross, minspan = minspan, endspan = endspan, 
															glm = list(family = binomial))
	}
	pred <- predict(fit.earth, newdata = newX, type = "response")
	fit <- list(object = fit.earth)
	out <- list(pred = pred, fit = fit)
	class(out$fit) <- c("SL.earth")
	return(out)
}


SL.ranger2 <- function (Y, X, newX, family, obsWeights, num.trees = 500, mtry = floor(sqrt(ncol(X))), 
												write.forest = TRUE, probability = family$family == "binomial", 
												min.node.size = ifelse(family$family == "gaussian", 5, 1), 
												replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632), 
												num.threads = 1, verbose = T, ...) {
	X = X[,which(names(X) %in% X_cat.names)]
	
	require("ranger")
	if (family$family == "binomial") {
		Y = as.factor(Y)
	}
	if (is.matrix(X)) {
		X = data.frame(X)
	}
	fit <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, X), 
												num.trees = num.trees, mtry = mtry, min.node.size = min.node.size, 
												replace = replace, sample.fraction = sample.fraction, 
												case.weights = obsWeights, write.forest = write.forest, 
												probability = probability, num.threads = num.threads, 
												verbose = verbose)
	pred <- predict(fit, data = newX)$predictions
	if (family$family == "binomial") {
		pred = pred[, "1"]
	}
	fit <- list(object = fit, verbose = verbose)
	class(fit) <- c("SL.ranger")
	out <- list(pred = pred, fit = fit)
	return(out)
}


SL.knn2 <- function (Y, X, newX, family, k = 10, ...) {
	X = X[,which(names(X) %in% X.names)]
	newX = newX[,which(names(newX) %in% X.names)]
	
	if (!is.matrix(X)) {
		X = model.matrix(~ -1 + ., data = X)
	}
	
	if (!is.matrix(newX)) {
		newX = model.matrix(~ -1 + ., data = newX)
	}
	
	require("class")
	if (family$family == "gaussian") {
		stop("SL.knn only available for family = binomial()")
	}
	fit.knn <- class::knn(train = X, test = newX, k = k, cl = Y, 
												prob = TRUE)
	pred <- (as.numeric(fit.knn) - 1) * attr(fit.knn, "prob") + 
		(1 - (as.numeric(fit.knn) - 1)) * (1 - attr(fit.knn, 
																								"prob"))
	fit <- list(k = k)
	out <- list(pred = pred, fit = fit)
	class(out$fit) <- c("SL.knn")
	return(out)
}

predict.SL.knn <- function (object, newdata, X, Y, ...) {
	X = X[,which(names(X) %in% X.names)]
	newdata = newdata[,which(names(newdata) %in% X.names)]
	
	if (!is.matrix(X)) {
		X = model.matrix(~ -1 + ., data = X)
	}
	
	if (!is.matrix(newdata)) {
		newdata = model.matrix(~ -1 + ., data = newdata)
	}
	
	require("class")
	fit.knn <- class::knn(train = X, test = newdata, k = object$k, 
												cl = Y, prob = TRUE)
	pred <- (as.numeric(fit.knn) - 1) * attr(fit.knn, "prob") + 
		(1 - (as.numeric(fit.knn) - 1)) * (1 - attr(fit.knn, 
																								"prob"))
	return(pred)
}

SL.ksvm2 <- function (Y, X, newX, family, type = NULL, kernel = "rbfdot", 
											kpar = "automatic", scaled = T, C = 1, nu = 0.2, epsilon = 0.1, 
											cross = 0, prob.model = family$family == "binomial", class.weights = NULL, 
											cache = 40, tol = 0.001, shrinking = T, ...) {
	X = X[,which(names(X) %in% X.names)]
	
	require("kernlab")
	if (!is.matrix(X)) {
		X = model.matrix(~., data = X)
		X = X[, -1]
	}
	
	if (family$family == "binomial") {
		Y = as.factor(Y)
		predict_type = "probabilities"
	}
	else {
		predict_type = "response"
	}
	model = kernlab::ksvm(X, Y, scaled = scaled, type = type, 
												kernel = kernel, kpar = kpar, C = C, nu = nu, epsilon = epsilon, 
												prob.model = prob.model, class.weights = class.weights)
	newX = newX[,which(names(newX) %in% X.names)]
	if (!is.matrix(newX)) {
		newX = model.matrix(~., data = newX)
		newX = newX[, -1, drop = FALSE]
	}
	pred = kernlab::predict(model, newX, predict_type)
	if (family$family == "binomial") {
		pred = pred[, 2]
	}
	fit = list(object = model, family = family)
	out = list(pred = pred, fit = fit)
	class(out$fit) = "SL.ksvm"
	return(out)
}

predict.SL.ksvm <- function (object, newdata, family, coupler = "minpair", ...) {
	newdata = newX[,which(names(newdata) %in% X.names)]
	
	require("kernlab")
	if (!is.matrix(newdata)) {
		newdata = model.matrix(~., data = newdata)
		newdata = newdata[, -1, drop = FALSE]
	}
	if (family$family == "binomial") {
		predict_type = "probabilities"
	}
	else {
		predict_type = "response"
	}
	pred = kernlab::predict(object$object, newdata, predict_type, 
													coupler = coupler)
	if (family$family == "binomial") {
		pred = pred[, 2]
	}
	return(pred)
}

SL.svm2 <- function (Y, X, newX, family, type.reg = "nu-regression", type.class = "nu-classification", 
										 kernel = "radial", nu = 0.5, degree = 3, cost = 1, coef0 = 0, 
										 ...) {
	X = X[,which(names(X) %in% X_cat.names)]
	newX = newX[,which(names(newX) %in% X_cat.names)]
	
	if (!is.matrix(X)) {
		X = model.matrix(~ -1 + ., data = X)
	}
	
	if (!is.matrix(newX)) {
		newX = model.matrix(~ -1 + ., data = newX)
	}
	
	require("e1071")
	if (family$family == "gaussian") {
		fit.svm <- e1071::svm(y = Y, x = X, nu = nu, type = type.reg, 
													fitted = FALSE, kernel = kernel, degree = degree, 
													cost = cost, coef0 = coef0)
		pred <- predict(fit.svm, newdata = newX)
		fit <- list(object = fit.svm)
	}
	if (family$family == "binomial") {
		fit.svm <- e1071::svm(y = as.factor(Y), x = X, nu = nu, 
													type = type.class, fitted = FALSE, probability = TRUE, 
													kernel = kernel, degree = degree, cost = cost, coef0 = coef0)
		pred <- attr(predict(fit.svm, newdata = newX, probability = TRUE), 
								 "prob")[, "1"]
		fit <- list(object = fit.svm)
	}
	out <- list(pred = pred, fit = fit)
	class(out$fit) <- c("SL.svm")
	return(out)
}

predict.SL.svm <- function (object, newdata, family, ...) {
	newdata = newdata[,which(names(newdata) %in% X_cat.names)]
	if (!is.matrix(newX)) {
		newdata = model.matrix(~ -1 + ., data = newdata)
	}
	require("e1071")
	if (family$family == "gaussian") {
		pred <- predict(object$object, newdata = newdata)
	}
	if (family$family == "binomial") {
		pred <- attr(predict(object$object, newdata = newdata, 
												 probability = TRUE), "prob")[, "1"]
	}
	return(pred)
}


SL.rpart2 <- function (Y, X, newX, family, obsWeights, cp = 0.01, minsplit = 20, 
											 xval = 0L, maxdepth = 30, minbucket = round(minsplit/3), 
											 ...) {
	X = X[,which(names(X) %in% X_cat.names)]
	
	require("rpart")
	if (family$family == "gaussian") {
		fit.rpart <- rpart::rpart(Y ~ ., data = data.frame(Y, 
																											 X), control = rpart::rpart.control(cp = cp, minsplit = minsplit, 
																											 																	 xval = xval, maxdepth = maxdepth, minbucket = minbucket), 
															method = "anova", weights = obsWeights)
		pred <- predict(fit.rpart, newdata = newX)
	}
	if (family$family == "binomial") {
		fit.rpart <- rpart::rpart(Y ~ ., data = data.frame(Y, 
																											 X), control = rpart::rpart.control(cp = cp, minsplit = minsplit, 
																											 																	 xval = xval, maxdepth = maxdepth, minbucket = minbucket), 
															method = "class", weights = obsWeights)
		pred <- predict(fit.rpart, newdata = newX)[, 2]
	}
	fit <- list(object = fit.rpart)
	out <- list(pred = pred, fit = fit)
	class(out$fit) <- c("SL.rpart")
	return(out)
}


SL.rpartPrune2 <- function (Y, X, newX, family, obsWeights, cp = 0.001, minsplit = 20, 
														xval = 10, maxdepth = 20, minbucket = 5, ...) {
	X = X[,which(names(X) %in% X_cat.names)]
	
	require("rpart")
	if (family$family == "gaussian") {
		fit.rpart <- rpart::rpart(Y ~ ., data = data.frame(Y, 
																											 X), control = rpart::rpart.control(cp = cp, minsplit = minsplit, 
																											 																	 xval = xval, maxdepth = maxdepth, minbucket = minbucket), 
															method = "anova", weights = obsWeights)
		CP <- fit.rpart$cptable[which.min(fit.rpart$cptable[, 
																												"xerror"]), "CP"]
		fitPrune <- rpart::prune(fit.rpart, cp = CP)
		pred <- predict(fitPrune, newdata = newX)
	}
	if (family$family == "binomial") {
		fit.rpart <- rpart::rpart(Y ~ ., data = data.frame(Y, 
																											 X), control = rpart::rpart.control(cp = cp, minsplit = minsplit, 
																											 																	 xval = xval, maxdepth = maxdepth, minbucket = minbucket), 
															method = "class", weights = obsWeights)
		CP <- fit.rpart$cptable[which.min(fit.rpart$cptable[, 
																												"xerror"]), "CP"]
		fitPrune <- rpart::prune(fit.rpart, cp = CP)
		pred <- predict(fitPrune, newdata = newX)[, 2]
	}
	fit <- list(object = fitPrune, fit = fit.rpart, cp = CP)
	out <- list(pred = pred, fit = fit)
	class(out$fit) <- c("SL.rpart")
	return(out)
}


SL.xgboost2 <- function (Y, X, newX, family, obsWeights, id, ntrees = 1000, 
												 max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(), 
												 nthread = 1, verbose = 0, save_period = NULL, ...) {
	X = X[,which(names(X) %in% X_cat.names)]
	newX = newX[,which(names(newX) %in% X_cat.names)]
	
	require("xgboost")
	if (packageVersion("xgboost") < 0.6) 
		stop("SL.xgboost requires xgboost version >= 0.6, try help('SL.xgboost') for details")
	if (!is.matrix(X)) {
		X = model.matrix(~. - 1, X)
	}
	xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
	if (family$family == "gaussian") {
		model = xgboost::xgboost(data = xgmat, objective = "reg:linear", 
														 nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
														 eta = shrinkage, verbose = verbose, nthread = nthread, 
														 params = params, save_period = save_period)
	}
	if (family$family == "binomial") {
		model = xgboost::xgboost(data = xgmat, objective = "binary:logistic", 
														 nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
														 eta = shrinkage, verbose = verbose, nthread = nthread, 
														 params = params, save_period = save_period)
	}
	if (family$family == "multinomial") {
		model = xgboost::xgboost(data = xgmat, objective = "multi:softmax", 
														 nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
														 eta = shrinkage, verbose = verbose, num_class = length(unique(Y)), 
														 nthread = nthread, params = params, save_period = save_period)
	}
	if (!is.matrix(newX)) {
		newX = model.matrix(~. - 1, newX)
	}
	pred = predict(model, newdata = newX)
	fit = list(object = model)
	class(fit) = c("SL.xgboost")
	out = list(pred = pred, fit = fit)
	return(out)
}

predict.SL.xgboost <- function (object, newdata, family, ...) {
	newdata = newdata[,which(names(newdata) %in% X_cat.names)]
	if (!is.matrix(newdata)) {
		newdata = model.matrix(~ -1 + ., newdata)
	}
		require("xgboost")
	if (packageVersion("xgboost") < 0.6) 
		stop("SL.xgboost requires xgboost version >= 0.6, try help('SL.xgboost') for details")
	if (!is.matrix(newdata)) {
		newdata = model.matrix(~. - 1, newdata)
	}
	pred = predict(object$object, newdata = newdata)
	return(pred)
}