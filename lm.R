simple_h = function(lm){
	n = length(lm$model[[1]])
	mean = mean(lm$model[[2]])
	var = var(lm$model[[2]])
	h = 1 / n + (lm$model[[2]] - mean) ^ 2 / ((n - 1) * var)
	return (h)
}

simple_standard_residual = function(lm){
	residual = lm$residuals
	se = sigma(lm)
	h = simple_h(lm)
	sol = residual / (se * sqrt(1 - h))
	return (sol)
}

simple_outlier = function(lm){
	sr = simple_standard_residual(lm)
	sol = which(sr < -2 | sr > 2)
	names(sol) = NULL
	sol = data.frame("No." = sol, "Residual" = sr[sol])
	row.names(sol) = NULL
	return (sol)
}

simple_influential = function(lm){
	n = length(lm$model[[1]])
	h = simple_h(lm)
	sol = which(h > 6 / n)
	sol = data.frame("No." = sol, "h" = h[sol])
	row.names(sol) = NULL
	return (list(critical = 6 / n, solution = sol))
}

simple_residual_normality = function(lm, test = c("chisq", "shapiro")){
	test = match.arg(test)
	sr = simple_standard_residual(lm)
	if(test == "shapiro"){
		return (shapiro.test(sr)$p.value)
	}
}

simple_plot_scedasticity = function(lm, ...){
	fit = lm$fitted.values
	sr = simple_standard_residual(lm)
	plot(x = fit, y = sr, xlab = "Predicted", ylab = "Residual", main = "I'm Scatter", ...)
}

simple_independ = function(lm, test = c("run")){
	test = match.arg(test)
	sr = simple_standard_residual(lm)
	if(test == "run"){
		library(snpar)
		return (runs.test(sr, exact = FALSE)$p.value)
	}
}

simple_linear_regress = function(formula, unit = NULL, fit = FALSE, normality = "shapiro", independ = "run", lcol = "black", ...){
	lm = lm(formula)
	sr = simple_standard_residual(lm)
	fitted = lm$fitted.values
	h = simple_h(lm)
	outlier = simple_outlier(lm)
	outlier_vector = outlier$No.
	influent = simple_influential(lm)
	influent_vector = influent$No.
	name_y = deparse(formula[2])
	name_x = deparse(formula[3])
	name_y = substr(name_y, start = 1, stop = nchar(name_y) - 2)
	name_x = substr(name_x, start = 1, stop = nchar(name_x) - 2)
	xlab = paste("Independent Variable", name_x)
	ylab = paste("Dependent Variable", name_y)
	if(!is.null(unit)){
		if(!is.na(unit[1])){
			ylab = paste(ylab, " (" , unit[1], ")", sep = "")
		}
		if(length(unit) >= 2 && !is.na(unit[2])){
			xlab = paste(xlab, " (" , unit[2], ")", sep = "")
		}
	}
	plot(y = lm$model[[1]], x = lm$model[[2]], main = "Scatter before fit", xlab = xlab, ylab = ylab, ...)
	abline(lm, col = lcol)
	g1 = recordPlot()
	remove = union(outlier_vector, influent_vector)
	while(fit && length(remove) > 0){
		y = lm$model[[1]][-remove]
		x = lm$model[[2]][-remove]
		lm = lm(y ~ x)
		outlier2 = simple_outlier(lm)
		outlier_vector = outlier$No.
		influent2 = simple_influential(lm)
		influent_vector = influent$No.
		remove = union(outlier_vector, influent_vector)
		if(length(remove) == 0){
			names(lm$coefficients)[2] = name_x
			names(lm$effects)[2] = name_x
			dimnames(lm$qr$qr)[[2]][2] = name_x
			names(lm$model) = c(name_y, name_x)
		}
	}
	plot(y = lm$model[[1]], x = lm$model[[2]], main = "Scatter after fit", xlab = xlab, ylab = ylab, ...)
	abline(lm, col = lcol)
	g2 = recordPlot()
	hist(sr, main = "Histogram of Standard Residuals")
	g3 = recordPlot()
	plot(y = sr, x = fitted, main = "Plot of Residuals vs Predicted", xlab = "Predicted", ylab = "Residuals", ...)
	g4 = recordPlot()
	plot(y = sr, x = 1:length(sr), main = "Autocorrelation", ylab = "Residual", xlab = "Times")
	g_autocor = recordPlot()
	d = Durbin_Watson(sr)
	dev.off()
	library(lmtest)
	return (list(model = lm, standard_residuals = sr, h = h, outlier = outlier, influent = influent, p_value_for_normality = simple_residual_normality(lm, normality)
		, p_value_for_independ = simple_independ(lm, independ), graph_for_unfit = g1, graph_for_fit = g2, Hist = g3, Plot = g4, bptest = bptest(lm)
		, syntax = paste("y = ", lm$coefficients[2], "x + (", lm$coefficients[1], ")", sep = ""), anova = anova(lm), autocorrelation = g_autocor, d = d))
}

multiple_h = function(lm){
	n = length(lm$model[[1]])
	k = length(lm$model) - 1
	X = matrix(nrow = n, ncol = k + 1)
	X[, 1] = rep_len(1, n)
	for(i in 2:(k + 1)){
		X[, i] = lm$model[[i]]
	}
	H = X %*% solve(t(X) %*% X) %*% t(X)
	return (H)
}

multiple_standard_residual = function(lm){
	residual = lm$residuals
	se = sigma(lm)
	h = diag(multiple_h(lm))
	sol = residual / (se * sqrt(1 - h))
	return (sol)
}

multiple_outlier = function(lm){
	sr = multiple_standard_residual(lm)
	sol = which(sr < -2 | sr > 2)
	names(sol) = NULL
	sol = data.frame("No." = sol, "Residual" = sr[sol])
	row.names(sol) = NULL
	return (sol)
}

multiple_influential = function(lm){
	k = lm$rank - 1
	n = length(lm$model[[1]])
	h = diag(multiple_h(lm))
	sol1 = which(h > 3 * (k + 1) / n)
	sol1 = data.frame("No." = sol1, "h" = h[sol1])
	row.names(sol1) = NULL
	Dh = h / (1 - h) ^ 2 * (lm$model[[1]] - lm$fitted.values) ^ 2 / sigma(lm) ^ 2
	sol2 = which(Dh / (k - 1) > 1)
	sol2 = data.frame("No." = sol2, "D" = Dh[sol2] / (k - 1))
	sol3 = which(Dh / k > 1)
	sol3 = data.frame("No." = sol3, "D" = Dh[sol3] / k)
	return (list(h = list(critical = 3 * (k + 1) / n, solution = sol1), D_k_1 = list(critical = 1, solution = sol2), D_k = list(critical = 1, solution = sol3)))
}

multiple_residual_normality = function(lm, test = c("chisq", "shapiro", "lilliefors")){
	test = match.arg(test)
	sr = multiple_standard_residual(lm)
	if(test == "shapiro"){
		return (shapiro.test(sr)$p.value)
	}else if(test == "lilliefors"){
		library(nortest)
		return (lillie.test(sr)$p.value)
	}
}

multiple_plot_scedasticity = function(lm, ...){
	fit = lm$fitted.values
	sr = multiple_standard_residual(lm)
	plot(x = fit, y = sr, xlab = "Predicted", ylab = "Residual", main = "I'm Scatter", ...)
}

multiple_independ = function(lm, test = c("run")){
	test = match.arg(test)
	sr = multiple_standard_residual(lm)
	if(test == "run"){
		library(snpar)
		return (runs.test(sr, exact = FALSE)$p.value)
	}
}

multiple_anova = function(lm){
	sol = matrix(ncol = 5, nrow = 3, dimnames = list(c("Regression", "Residual", "Total"), c("DoF", "SQ", "MS", "F", "pvalue")))
	n = length(lm$model[[1]])
	k = lm$rank - 1
	sol[1, 1] = k
	sol[2, 1] = n - k - 1
	sol[3, 1] = n - 1
	SSE = sigma(lm) ^ 2 * (n - k - 1)
	SST = sum((lm$model[[1]] - mean(lm$model[[1]])) ^ 2)
	SSR = SST - SSE
	sol[1, 2] = SSR
	sol[2, 2] = SSE
	sol[3, 2] = SST
	MSR = SSR / k
	MSE = SSE / (n - k - 1)
	sol[1, 3] = MSR
	sol[2, 3] = MSE
	F = MSR / MSE
	sol[1, 4] = F
	sol[1, 5] = 1 - pf(F, k, n - k - 1)
	return(as.table(sol))
}

Durbin_Watson = function(sr){
	sr2 = sum(sr ^ 2)
	dif = sr
	dif[1] = 0
	for(i in 2:length(sr)){
		dif[i] = sr[i] - sr[i - 1]
	}
	return (sum(dif ^ 2) / sr2)
}

multiple_linear_regress = function(formula, unit = NULL, fit = FALSE, normality = "shapiro", independ = "run", lcol = "black", ...){
	lm = lm(formula)
	k = lm$rank - 1
	sr = multiple_standard_residual(lm)
	exp = lm$fitted.values
	h = diag(multiple_h(lm))
	outlier = multiple_outlier(lm)
	outlier_vector = outlier$No.
	influent = multiple_influential(lm)
	influent_vector = union(influent$h$No., influent$D_k_1$No.)
	name_y = deparse(formula[2])
	name_x = deparse(formula[3])
	name_x = substr(name_x, 2, nchar(name_x) - 3)
	name_x = strsplit(name_x, " + ", fixed = TRUE)
	name_x = name_x[[1]]
	name_y = substr(name_y, start = 1, stop = nchar(name_y) - 2)
	coref = matrix(nrow = k + 1, ncol = k + 1, dimnames = list(c(name_y, name_x), c(name_y, name_x)))
	for(i in 1:(k + 1)){
		for(j in i:(k + 1)){
			coref[j, i] = cor(lm$model[[i]], lm$model[[j]])
		}
	}
	xlab = paste("Independent Variable", name_x)
	ylab = paste("Dependent Variable", name_y)
	if(!is.null(unit)){
		if(!is.na(unit[1])){
			ylab = paste(ylab, " (" , unit[1], ")", sep = "")
		}
		tmp = min(length(unit) - 1, k)
		for(i in 1:tmp){
			if(!is.na(unit[1 + i])){
				xlab[i] = paste(xlab[i], " (" , unit[1 + i], ")", sep = "")
			}
		}
	}
	g_before = list()
	for(i in 1:k){
		plot(y = lm$model[[1]], x = lm$model[[1 + i]], main = "Scatter before fit", xlab = xlab[i], ylab = ylab, ...)
		abline(lm(lm$model[[1]] ~ lm$model[[1 + i]]), col = lcol)
		g_before[[i]] = recordPlot()
	}
	names(g_before) = name_x
	remove = union(outlier_vector, influent_vector)
	while(fit && length(remove) > 0){
		y = lm$model[[1]][-remove]
		x = list()
		str = "lm = lm(y ~ "
		for(i in 1:k){
			x[[i]] = lm$model[[1 + i]][-remove]
			if(i == 1){
				str = paste(str, "x[[", i, "]]", sep = "")
			}else{
				str = paste(str, " + x[[", i, "]]", sep = "")
			}
		}
		str = paste(str, ")", sep = "")
		eval(parse(text = str))
		outlier = multiple_outlier(lm)
		outlier_vector = outlier$No.
		influent = multiple_influential(lm)
		influent_vector = union(influent$h$No., influent$D_k_1$No.)
		remove = union(outlier_vector, influent_vector)
		if(length(remove) == 0){
			for(i in 1:k){
				names(lm$coefficients)[1 + i] = name_x[i]
				names(lm$effects)[1 + i] = name_x[i]
				dimnames(lm$qr$qr)[[2]][1 + i] = name_x[i]
			}
			names(lm$model) = c(name_y, name_x)
		}
	}
	g_after = list()
	for(i in 1:k){
		plot(y = lm$model[[1]], x = lm$model[[1 + i]], main = "Scatter after fit", xlab = xlab[i], ylab = ylab, ...)
		abline(lm(lm$model[[1]] ~ lm$model[[1 + i]]), col = lcol)
		g_after[[i]] = recordPlot()
	}
	names(g_before) = name_x
	hist(sr, main = "Histogram of Standard Residuals")
	g3 = recordPlot()
	plot(y = sr, x = exp, main = "Plot of Residuals vs Predicted", xlab = "Predicted", ylab = "Residuals", ...)
	g4 = recordPlot()
	dev.off()
	syntax = "y = "
	for(i in 1:k){
		if(i == 1){
			syntax = paste(syntax, lm$coefficients[1 + i], "x", i, sep = "")
		}else{
			syntax = paste(syntax, " + (", lm$coefficients[1 + i], ")x", i, sep = "")
		}
	}
	syntax = paste(syntax, " + (", lm$coefficients[1], ")", sep = "")
	plot(y = sr, x = 1:length(sr), main = "Autocorrelation", ylab = "Residual", xlab = "Times", ...)
	g_autocor = recordPlot()
	d = Durbin_Watson(sr)
	dev.off()
	library(lmtest)
	return (list(model = lm, standard_residuals = sr, h = h, outlier = outlier, influent = influent, p_value_for_normality = multiple_residual_normality(lm, normality)
		, p_value_for_independ = multiple_independ(lm, independ), graph_for_unfit = g_before, graph_for_fit = g_after, Hist = g3, Plot = g4, bptest = bptest(lm)
		, syntax = syntax, correlation = as.table(coref), anova = multiple_anova(lm), autocorrelation = g_autocor, d = d))
}

Scatt_with_category = function(y, x, f, color = c("red", "blue"), ...){
	kind_table = table(f)
	kind = length(kind_table)
	color = rep_len(color, kind)
	plot(x = x, y = y)
	t = ""
	for(i in 1:kind){
		index = which(f == names(kind_table)[i])
		points(x = x[index], y = y[index], col = color[i], ...)
		abline(lm(y[index] ~ x[index]), col = color[i])
		t = paste(t, names(kind_table)[i], " : ", color[i], "\n", sep = "")
	}
	g = recordPlot()
	sol = list(Graph = g, Text = t)
}
