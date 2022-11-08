# mortality.R
# Left truncation weighting - does it "work" for mortality?
# December 7, 2020
# Kevin Chen

# Aims ####
# 1. Restrict FU to start in 1985
# 2. Compute left-truncation weights
# 3. Run weighted analysis

library(here)

og.dir <- getwd()
detach("package:here", unload = T)
setwd("../gm-cancer-mort")
source(here::here("modeling.R"))
cohort_analytic_21 <- cohort_analytic[year <= 1999]
detach("package:here", unload = T)
setwd(og.dir)
library(here)

# Get analytic data for other exposure lags
sapply(c(5, 10), function(exposure.lag) {
	if (!(paste0('cohort_analytic', "_", exposure.lag) %in% ls(envir = .GlobalEnv))) {
		cohort_analytic <- get.cohort_analytic(
			outcome_type = outcome.type,
			exposure.lag = exposure.lag,
			deathage.max = NULL
		)
		setorder(cohort_analytic, studyno, year)
		cohort_analytic <- cohort_analytic[year >= 1941 & (
			year(yin) < 1938 | year >= year(yin + 365.25 * 3)
		)]

		# Remove people with unknown cause of death
		cohort_analytic <- cohort_analytic[!(status15 == 6 & (is.na(icd)))]

		# Limit FU
		cohort_analytic <- cohort_analytic[year <= 1999]
	}

	assign(paste0('cohort_analytic', "_", exposure.lag),
				 cohort_analytic, envir = .GlobalEnv)
})

# Get left-truncation weights
weights.year.min <- 1941
if (weights.year.min == 1941) {
	# Survival: 1941 to 1985
	survival_to_85 <- box_read(674612305081)}
if (weights.year.min == 1970) {
	# Survival: 1970 to 1985
	survival_to_85 <- box_read(671315790911)
}
setDT(survival_to_85)

# Modeling helper function
run.mod <- function(outcome = "Non-Hodgkin lymphoma",
										cohort_analytic.name = "cohort_analytic_5",
										new_data = T,
										run_model = F,
										weight_type = c("none", "p", "sw"),
										full_FU = F
) {

	if (!(paste0(gsub(" ", "_", outcome ), ".cohort_prepped") %in% ls(envir = .GlobalEnv)) |
			new_data == T) {
		get.mod(
			cohort_py = get(cohort_analytic.name, envir = .GlobalEnv),
			outcome = outcome,
			probs = probs,
			outcome_type = outcome.type,
			run_coxph = F
		)}

	tmp.cohort_prepped <- get(paste0(
		gsub(" ", "_", outcome ), ".cohort_prepped"
	))

	# Arbitrary MWF categories for straight
	tmp.cohort_prepped[,`:=`(
		Straight.cat = cut(Straight, unique(c(-Inf, 0, 2, max(Straight))))
	)]

	if (!full_FU) {

		tmp.cohort_prepped <- tmp.cohort_prepped[Year >= 1985]

		tmp.cohort_prepped <- merge(
			tmp.cohort_prepped,
			survival_to_85[dead_by_85 == 0, .(
				studyno,
				p = big.gam.prob)],
			on = "studyno",
			all.x = T)

		# If hired in 1985, go ahead and give them full weight
		tmp.cohort_prepped[is.na(p) & year(`Year of hire`) + 3 > 1985, p := 1]
		tmp.cohort_prepped <- tmp.cohort_prepped[!is.na(p)]

		overall.p <- 1 - mean(tmp.cohort_prepped[,.(p = p[1]), by = .(studyno)]$p)

		tmp.cohort_prepped[,`:=`(
			sw = overall.p/(1 - p)
		), by = .(studyno)]

		# Truncate at 90th percentile
		sw.max <- quantile(tmp.cohort_prepped[, .(sw = sw[1]), by = studyno]$sw,
											 .99)
		tmp.cohort_prepped[sw > sw.max,`:=`(
			sw = sw.max
		), by = .(studyno)]

		tmp.cohort_prepped[,`:=`(
			sw = c(sw[1], rep(1, .N - 1))
		), by = .(studyno)]
	} else {
		weight_type = "none"
	}


	# Run model
	basic.formula <- 	paste(
		"Surv(age.year1, age.year2, event) ~",
		if (n_distinct(tmp.cohort_prepped$Straight.cat) > 1) {
			"Straight.cat +"
		},
		if (n_distinct(tmp.cohort_prepped$Soluble.cat) > 1) {
			"Soluble5.cat +"
		},
		if (n_distinct(tmp.cohort_prepped$Synthetic.cat) > 1) {
			"Synthetic.cat +"
		},
		if (min(table(as.numeric(tmp.cohort_prepped[event == 1]$Year.cat))) > 1 &
				length(unique(tmp.cohort_prepped[event == 1]$Year.cat)) > 1) {
			"Year.cat +"
		},
		if (min(table(as.numeric(tmp.cohort_prepped[event == 1]$`Year of hire.cat`))) > 1 &
				length(unique(tmp.cohort_prepped[event == 1]$`Year of hire.cat`)) > 1) {
			"`Year of hire.cat` +"
		},
		if (min(table(tmp.cohort_prepped[event == 1]$`Race`)) > 10 &
				length(unique(tmp.cohort_prepped[event == 1]$`Race`)) > 1) {
			"Race +"
		},
		if (min(table(tmp.cohort_prepped[event == 1]$Sex)) > 10 &
				length(unique(tmp.cohort_prepped[event == 1]$Sex)) > 1) {
			"Sex +"
		},
		'Plant'
	)

	# Cox models
	sapply(weight_type, function(weight.i = "none") {

		if (weight.i != "none") {
			tmp.cohort_prepped$weight <- tmp.cohort_prepped[, weight.i, with = F]
		}

		if (run_model) {
			tmp.coxph <- coxph(as.formula(basic.formula),
												 data = tmp.cohort_prepped,
												 ties = 'efron',
												 cluster = studyno,
												 weights = if (weight.i != "none") {weight
												 } else {NULL})

			mod.dir <- here::here("resources/mortality")
			mod.name <- paste0(gsub(" ", "_", outcome), "_", weight.i,
												 ifelse(full_FU,  "_full-FU", ""),
												 ".coxph.rds")

			# Save model
			saveRDS(tmp.coxph, file = gsub("//", "/", paste0(mod.dir, "/", mod.name)))
		} else {
			tmp.coxph <- readRDS(file = gsub("//", "/", paste0(mod.dir, "/", mod.name)))
		}

		print(summary(tmp.coxph)$coefficients[,c(2,4,5)])

		cbind(summary(tmp.coxph)$coefficients[, -c(1, 3, 4)],
					summary(tmp.coxph)$conf.int[, 3:4]) %>%
			as.data.frame() %>% cbind(
				name = substr(rownames(.), 1, unlist(
					gregexpr(
						'\\.cat`\\(|\\.cat\\(|Not |Female$|2$|3$|Hydra|Sagin',
						rownames(.)
					)
				) - 1),
				level = gsub(',', ", ", substring(rownames(.), unlist(
					gregexpr('t`\\(|t\\(|eNot|xF|t2|t3|tHydra|tSagin', rownames(.))
				) + 1)),
				n = NaN,
				.
			) %>% mutate(" " = ifelse(
				`Pr(>|z|)` < 0.05, "$*$",
				ifelse(`Pr(>|z|)` < 0.1, "$\\cdot$", NA))) %>% apply(
					2, function(x) {
						gsub("`", "", as.character(x))
					}) %>% as.data.frame -> coxph.tab

		coxph.tab[, 3:(ncol(coxph.tab) - 1)] <- apply(coxph.tab[, 3:(ncol(coxph.tab) - 1)], 2, as.numeric)

		mwf.vec <- NULL
		messy_ref <- '_sol5'

		coxph.tab <- coxph.tab[grep(
			ifelse(is.null(mwf.vec), "^strai|^sol|^synth", paste0('^', mwf)),
			coxph.tab$name, ignore.case = T),]

		coxph.tab$name <- with(coxph.tab, gsub("500|10|5", "", name))

		if (is.null(mwf.vec)) {
			coxph.tab$n <- as.vector(unlist(
				sapply(tmp.cohort_prepped[event == 1, .(
					Straight.cat,
					get(paste0('Soluble',
										 ifelse(
										 	!is.null(messy_ref),
										 	substring(messy_ref, 5), ''), '.cat')),
					Synthetic.cat)], function(x) {
						table(x)[-1]
					})
			))
		} else {
			coxph.tab$n <- as.vector(unlist(
				sapply(tmp.cohort_prepped[
					event == 1,
					paste0(mwf, ifelse(
						paste0(mwf) == "Soluble" & !is.null(messy_ref),
						substring(messy_ref, 5), ""), '.cat'),
					with = F], function(x) {
						table(x)[-1]
					})
			))
		}

		setDT(coxph.tab)

		coxph.tab[!is.na(`lower .95`), `95\\% CI` := paste0(
			"(",
			formatC(`lower .95`, digits = 2, format = 'f'),
			", ",
			formatC(`upper .95`, digits = 2, format = 'f'),
			")"
		)]

		referent.tab <- data.table()
		referent.tab[,	(names(coxph.tab)) := list(
			# name
			if (is.null(mwf.vec)) {
				c("Straight", "Soluble", "Synthetic")
			} else {paste(mwf)},
			# level
			{if (is.null(mwf.vec)) {sapply(tmp.cohort_prepped[event == 1, .(
				Straight.cat,
				get(paste0('Soluble', ifelse(
					!is.null(messy_ref),
					substring(messy_ref, 5), ""), '.cat')),
				Synthetic.cat)],
				function(x) {
					names(table(x)[1])
				})} else {
					sapply(tmp.cohort_prepped[
						event == 1, paste0(mwf, ifelse(
							!is.null(messy_ref) & mwf == "Soluble",
							substring(messy_ref, 5), ""), '.cat'),
						with = F],
						function(x) {
							names(table(x)[1])
						})
				}},
			# n
			{if (is.null(mwf.vec)) {sapply(tmp.cohort_prepped[event == 1, .(
				Straight.cat,
				get(paste0('Soluble', ifelse(
					!is.null(messy_ref),
					substring(messy_ref, 5), ""), '.cat')),
				Synthetic.cat)],
				function(x) {
					table(x)[1]
				})} else {
					sapply(
						tmp.cohort_prepped[
							event == 1,
							paste0(mwf, ifelse(
								!is.null(messy_ref) & mwf == "Soluble",
								substring(messy_ref, 5), ""), '.cat'), with = F],
						function(x) {
							table(x)[1]
						})
				}
			},
			# exp(coef)
			NA,
			# z
			NA,
			# Pr()
			NA,
			# lower .95
			NA,
			# upper .95
			NA,
			# " "
			NA,
			# CI
			NA
		)]

		referent.tab[, level := gsub(",", ", ", level)]

		coxph.tab <- rbindlist(list(coxph.tab, referent.tab))

		coxph.tab[, `:=`(
			exposure.lower = as.numeric(substr(
				level, 2, unlist(regexpr(",", level)) - 1)),
			exposure.upper = as.numeric(ifelse(grepl("]", level), substr(
				level, unlist(regexpr(", ", level)) + 2,
				unlist(regexpr("]", level)) - 1), NA))
		)]

		setorder(coxph.tab, name, exposure.lower)

		coxph.tab[, `:=`(
			level =
				if (.N > 2) {
					c(
						ifelse(exposure.upper[1] == 0, '$0$',
									 paste0("$0$ to $",
									 			 exposure.upper[1], "$"
									 )),
						paste0(paste0(
							"$>",
							{if (exposure.lower[2] == .05) {
								sapply(1:(.N - 2),
											 function(i) {
											 	ifelse(
											 		i == 1,
											 		exposure.lower[-c(1, .N)][i],
											 		round(exposure.lower[-c(1, .N)][i], 1))

											 })
							} else {
								round(exposure.lower[-c(1, .N)], 1)
							}}),
							"$ to $",
							round(
								exposure.upper[-c(1, .N)], 1), "$ "
						),
						paste0("$>", round(exposure.lower[.N], 1), "$")
					)} else {c(
						ifelse(exposure.upper[1] == 0, '$0$',
									 paste0("$0$ to $",
									 			 exposure.upper[1], "$ "
									 )), paste0("$>", round(exposure.lower[.N], 1), "$"))}
		), by = .(name)]

		# make column for mwf tpe
		coxph.tab[, mwf := name]
		# make name column pretty
		coxph.tab[duplicated(name), name := NA]

		# Order rows
		setorder(coxph.tab, mwf, exposure.lower)

		coxph.tab <- coxph.tab[,-c('exposure.lower', 'exposure.upper')]

		# Add column for units
		coxph.tab <- cbind(coxph.tab[,1:2],
											 " " = 'mg/m$^3\\cdot$years',
											 coxph.tab[,3:ncol(coxph.tab)])

		coxph.tab$`exp(coef)` <- formatC(round(coxph.tab$`exp(coef)`, digits = 2), format = 'f', digits = 2)
		coxph.tab$`exp(coef)`[grep('NA', coxph.tab$`exp(coef)`)] <- NA
		coxph.tab$`Pr(>|z|)` <- formatC(round(coxph.tab$`Pr(>|z|)`, digits = 2), format = 'f', digits = 2)

		coxph.tab$`Pr(>|z|)` <- gsub("1.00", "$0.99$", coxph.tab$`Pr(>|z|)`)

		coxph.tab$`Pr(>|z|)`[grep('NA', coxph.tab$`Pr(>|z|)`)] <- NA

		coxph.tab[, `:=`(outcome = paste(outcome))]

		names(coxph.tab) <- c('name', 'exposure',
													"   ", '$n$', 'HR', '$z$',
													'$p$', 'lower .95', 'upper .95',
													"    ", '95\\% CI', 'mwf', 'outcome')

		# Save coefficient table
		tab.dir <- here::here(paste0("resources/mortality/coefficient tables/FU to ",
																 max(tmp.cohort_prepped$Year)))
		if (max(tmp.cohort_prepped$Year) != 2015) {
			tab.dir <- paste(tab.dir, paste0("lag ", gsub("[a-z]|_", "", cohort_analytic.name)),
											 sep = "/")
		}
		tab.name <- paste0(gsub(" ", "_", outcome), "_", weight.i,
											 ifelse(full_FU,  "_full-FU", ""),
											 ".coef.rds")
		saveRDS(coxph.tab, file = gsub("//", "/", paste0(tab.dir, "/", tab.name)))

		return(coxph.tab)
	}, simplify = F)
}

# Run model
for (j in c(21, 10, 5)) {
	for (x in c(outcome.selected[c(2, 4, 11)], "Non-Hodgkin lymphoma")) {
		for (FU in c(T, F)) {
			run.mod(x,
							cohort_analytic.name = paste0("cohort_analytic_", j),
							run_model = T, full_FU = FU)
		}
	}
}

# Print results
print_tab <- function(
	outcomes = c(outcome.selected[c(2, 4, 11)], "Non-Hodgkin lymphoma"),
	scale = 0.7,
	FU = 1999,
	exposure.lag = 21) {
	sapply(outcomes, function(x = "Non-hodgkin lymphoma") {
		sapply(exposure.lag, function(lag) {
			# x <- outcome.selected[c(2, 4, 11)][1]

			tab.dir <- here::here(paste0("resources/mortality/coefficient tables",
																	 "/FU to ", FU,
																	 "/lag ", lag))

			original_coxph.tab <- readRDS(
				file = paste0(tab.dir, "/",
											paste0(gsub(" ", "_", x), "_none_full-FU.coef.rds")))

			original_coxph.tab[!is.na(name), name := paste0("\\multicolumn{3}{l}{", name, "}\\\\\n")]

			# Clean up names
			names(original_coxph.tab)[grep("name", names(original_coxph.tab))] <- 'Cumulative exposure'
			names(original_coxph.tab)[grep("exposure", names(original_coxph.tab))] <- '\\vspace{0em}'

			coxph.tab <- lapply(c("none", "p", "sw"), function(weight.i = "p") {
				tab.name <- paste0(gsub(" ", "_", x), "_", weight.i, ".coef.rds")
				coxph.tab <- readRDS(file = gsub("//", "/", paste0(tab.dir, "/", tab.name)))
				return(coxph.tab)
			})

			og_results <- as.data.table(as.data.frame(original_coxph.tab))
			og_results <- og_results[,c(1:2, 4:5, 11)]

			order_mwf <- function(name = og_results[,1]) {
				as.vector(unlist(sapply(c("Str", "Sol", "Syn"), function(mwf = "Str") {
					grep(mwf, unlist(zoo::na.locf(name)))
				})))}

			cbind(
				og_results[order_mwf(og_results[,1])],
				" " = NA,
				coxph.tab[[1]][,c(4:5, 11)][order_mwf(coxph.tab[[1]][,1])],
				" " = NA,
				coxph.tab[[2]][,c(5, 11)][order_mwf(coxph.tab[[2]][,1])],
				" " = NA,
				coxph.tab[[3]][,c(5, 11)][order_mwf(coxph.tab[[3]][,1])]
			) %>% xtable %>% print(
				hline.after = c(0, nrow(.)),
				add.to.row = list(
					pos = list(-1),
					command = c(
						paste0(
							"\\toprule ",
							"& & \\multicolumn{3}{c}{Full follow-up}",
							"& & \\multicolumn{3}{c}{1985 cohort}",
							"& & \\multicolumn{2}{c}{Weighted by surv.}",
							"& & \\multicolumn{2}{c}{St. inverse weight}",
							"\\\\\n",
							"\\cline{3-5} \\cline{7-9} \\cline{11-12} \\cline{14-15}",
							"\n")
					)
				)) -> tab.text

			return(paste0(
				paste0(
					"# \n\n",
					"## ", x, " (exposure lagged ", lag, " years)",
					"\n\n",
					"\\begin{table}\\linespread{1.15}\n",
					"\\begin{adjustbox}{scale = ", scale, "}\n",
					"%\\caption{",
					"Adjusted HR estimates for ", tolower(x), ".",
					"}\n"),
				tab.text,
				"\\end{adjustbox}",
				"\\end{table}\n\n"
			))
		}, USE.NAMES = F)
	}, USE.NAMES = F)
}

print_tab(exposure.lag = c(21, 10, 5)) %>% as.vector %>% gsub("\t", "", .) %>% clipr::write_clip()

# # Age at entry
# cohort_prepped <- cohort_analytic[
# 	immortal == 0 & nohist == 0 & wh == 1 & right.censored != 1,
# 	.(`Year of cohort entry` = year(year1[1]),
# 		`Age at cohort entry` = age.year1[1]/365,
# 		`Age in 1985` = age.year1[year == 1985]/365
# 		), by = .(studyno)]
#
# rbindlist(
# 	list(cohort_prepped[,
# 			 		.(studyno,
# 			 			`Year of cohort entry`,
# 			 			`Age at cohort entry`,
# 			 			population = "Full cohort",
# 			 			FU = "From 1941"
# 			 			)],
# 			 cohort_prepped[studyno %in% cohort_analytic[year >= 1985]$studyno,
# 			 		.(`Year of cohort entry` = max(1985, `Year of cohort entry`),
# 			 			`Age at cohort entry` = max(`Age in 1985`, `Age at cohort entry`),
# 			 			population = "Full cohort",
# 			 			FU = "From 1985"
# 			 			), by = .(studyno)],
# 			 rbindlist(
# lapply(c(outcome.selected[c(2, 4, 11)], "Non-Hodgkin lymphoma"),
# 			 function(outcome) {
# 			 	rbindlist(list(cohort_prepped[
# 			 		studyno %in% cohort_analytic[get(outcome) == 1]$studyno,
# 			 		.(studyno,
# 			 			`Year of cohort entry`,
# 			 			`Age at cohort entry`,
# 			 			population = outcome,
# 			 			FU = "From 1941")],
# 			 	cohort_prepped[
# 			 		studyno %in% cohort_analytic[
# 			 			get(outcome) == 1 & year >= 1985]$studyno,
# 			 		.(`Year of cohort entry` = max(1985, `Year of cohort entry`),
# 			 			`Age at cohort entry` = max(`Age in 1985`, `Age at cohort entry`),
# 			 			population = outcome,
# 			 			FU = "From 1985"),
# 			 		by = .(studyno)]))
# 			 })
# ))) -> ggtab
#
# tikz(here::here("reports/mortality sanity check/resources", "age-at-entry.tex"),
# 		 standAlone = T, width = 5 * 1.2, height = 3 * 1.2)
# ggtab %>% ggplot(aes(x = `Age at cohort entry`,
# 									 fill = FU,
# 									 color = FU)) +
# 	geom_density(alpha = 0.5) +
# 	facet_wrap(. ~ population, ncol = 3) +
# 	mytheme + labs(
# 		color = "Follow-up",
# 		fill = "Follow-up") +
# 	theme(legend.position = "bottom",
# 				legend.box.margin=margin(-10,5,5,5))
# dev.off()
# lualatex("age-at-entry\\.tex",
# 				 here::here("reports/mortality sanity check/resources"))