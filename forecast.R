print.FORECAST = function(x, ...){
	sol = x$solution["Forecast", ]
	print(as.table(sol))
}

summary.FORECAST = function(x, newline = FALSE, ...){
	sol = x
	class(sol) = "list"
	if(sol$sign == "Moving"){
		if(newline){
			sol$sign = paste("Moving\naverage\nwith", sol$parameter["period"], "period")
		}else{
			sol$sign = paste("Moving average with", sol$parameter["period"], "period")
		}
	}else if(sol$sign == "Exponent"){
		if(newline){
			sol$sign = paste("Exponential\nsmoothing\nwith alpha =", sol$parameter["alpha"])
		}else{
			sol$sign = paste("Exponential smoothing with alpha =", sol$parameter["alpha"])
		}
	}else if(sol$sign == "Holt\'s"){
		if(newline){
			sol$sign = paste("Holt\'s\nexponential\nsmoothing\nwith alpha =", sol$parameter["alpha"], "\nbeta =", sol$parameter["beta"]
			, "\nET0 =", sol$parameter["ET0"], "\nF0 =", sol$parameter["F0"], "\nX0 =", sol$parameter["X0"])
		}else{
			sol$sign = paste("Holt\'s exponential smoothing with alpha =", sol$parameter["alpha"], ", beta =", sol$parameter["beta"]
			, ", ET0 =", sol$parameter["ET0"], ", F0 =", sol$parameter["F0"], ", X0 =", sol$parameter["X0"])
		}
	}else if(sol$sign == "Linear Regression"){
		if(newline){
			sol$sign = paste("Linear\nRegression\n", sol$parameter$lm.call)
		}else{
			sol$sign = paste("Linear Regression", sol$parameter$lm.call)
		}
	}else if(sol$sign == "User-defined package"){
		if(newline){
			sol$sign = "User-defined\npackage"
		}else{
			sol$sign = "User-defined package"
		}
	}else if(sol$sign == "Average"){
		if(newline){
			sol$sign = "Average\nValue"
		}else{
			sol$sign = "Average Value"
		}
	}else if(sol$sign == "Seasonal"){
		if(newline){
			sol$sign = paste("Seasonal\nforecasting\nwith", sol["index"], "\nand", sol["forecast.way"], "\nfor period =", sol$period)
		}else{
			sol$sign = paste("Seasonal forecasting with", sol["index"], "and", sol["forecast.way"], "for period =", sol$period)
		}
	}
	names(sol)[3] = "type_of_forecasting"
	names(sol)[4] = "forecasting period"
	class(sol) = "summary.FORECAST"
	return(sol[-2])
}

mean.FORECAST = function(x){
	return(mean(x$solution["Forecast", ], na.rm = TRUE))
}

data_filt = function(x){
	if(class(x) == "FORECAST"){
		if(sum(dimnames(x$solution)[[1]] == "x") == 0){
			return(NULL)
		}else{
			x = x$solution["x", ]
		}
	}else if(class(x) == "CMA"){
		x = x$x
	}
	x = subset(x, !is.na(x))
	return(x)
}

average_value = function(x){
	x = data_filt(x)
	sol = matrix(ncol = length(x) + 1, nrow = 2, dimnames = list(c("x", "Forecast"), c(1:length(x), "forecast")))
	sol["x", ] = c(x, NA)
	accu = x[1]
	for(i in 2:length(x)){
		sol["Forecast", i] = accu / (i - 1)
		accu = accu + x[i]
	}
	sol["Forecast", length(x) + 1] = accu / length(x)
	sol = list(solution = as.table(sol), parameter = NULL, sign = "Average", period = 1)
	class(sol) = "FORECAST"
	return(sol)
}

moving_average = function(x, period = 3){
	x = data_filt(x)
	if(length(x) <= period){
		return (NULL)
	}
	sol = matrix(ncol = length(x) + 1, nrow = 2, dimnames = list(c("x", "Forecast"), c(1:length(x), "forecast")))
	sol["x", ] = c(x, NA)
	for(i in 1:(length(x) - period + 1)){
		sol["Forecast", period + i] = sum(x[i:(i + period - 1)]) / period
	}
	sol = list(solution = as.table(sol), parameter = c(period = period), sign = "Moving", period = 1)
	class(sol) = "FORECAST"
	return(sol)
}

exponential_smoothing = function(x, alpha = 0.2){
	x = data_filt(x)
	sol = matrix(ncol = length(x) + 1, nrow = 2, dimnames = list(c("x", "Forecast"), c(1:length(x), "forecast")))
	sol["x", ] = c(x, NA)
	sol["Forecast", 1] = x[1]
	for(i in 2:(length(x) + 1)){
		sol["Forecast", i] = alpha * x[i - 1] + (1 - alpha) * sol["Forecast", i - 1]
	}
	sol = list(solution = as.table(sol), parameter = c(alpha = alpha), sign = "Exponent", period = 1)
	class(sol) = "FORECAST"
	return(sol)
}

holt_exponential_smoothing = function(x, alpha = 0.9, beta = 0.8, ET0 = 1, F0 = mean(x), X0 = F0){
	x = data_filt(x)
	n = length(x)
	sol = matrix(ncol = n + 2, nrow = 4, dimnames = list(c("x", "Forecast", "LT", "ET"), c(0:n, "forecast")))
	sol["x", ] = c(NA, x, NA)
	sol["ET", '0'] = ET0
	sol["Forecast", '0'] = F0
	x = c(X0, x)
	for(i in (1:n) + 1){
		sol["Forecast", i] = alpha * x[i - 1] + (1 - alpha) * sol["Forecast", i - 1] + sol["ET", i - 1]
		sol["LT", i] = alpha * (x[i] - x[i - 1]) + (1 - alpha) * (sol["Forecast", i] - sol["Forecast", i - 1])
		sol["ET", i] = beta * sol["LT", i] + (1 - beta) * sol["ET", i - 1]
	}
	sol["Forecast", "forecast"] = alpha * x[n + 1] + (1 - alpha) * sol["Forecast", n + 1] + sol["ET", n + 1]
	f = function(x, period){
		param = x$parameter
		x = x$solution
		n = ncol(x) - 2
		nowx = x["x", n + 1]
		nowf = x["Forecast", n + 1]
		nowET = x["ET", n + 1]
		alpha = param["alpha"]
		sol = alpha * nowx + (1 - alpha) * nowf + period * nowET
		names(sol) = NULL
		return(sol)
	}
	sol = list(solution = as.table(sol), parameter = c(alpha = alpha, beta = beta, ET0 = ET0, F0 = F0, X0 = X0), sign = "Holt\'s", period = 1, forecast_function = f)
	class(sol) = "FORECAST"
	return(sol)
}

forecast = function(x, period){
	if(class(x) != "FORECAST" || period <= 0){
		return(NULL)
	}
	if(!is.null(x$forecast_function)){
		if(is.null(x$forecast_function.expression)){
			return(x$forecast_function(x, period))
		}else{
			return(x$forecast_function(x, period, x$forecast_function.expression))
		}
	}else{
		f = x$solution["Forecast", ]
		f = subset(f, names(f) == "forecast")
		if(length(f) < period){
			return(NULL)
		}
		sol = f[period]
		names(sol) = NULL
		return(sol)
	}
}

find_ori = function(data){
	if(class(data) != "list"){
		data = list(data)
	}
	ori = vector()
	for(i in data){
		if(class(i) == "FORECAST"){
			if(sum(dimnames(i$solution)[[1]] == "x") == 0){
				next
			}
			ori = i$solution["x", ]
			ori = na.omit(ori)
			break
		}else if(class(i) == "CMA"){
			ori = i$x
			ori = na.omit(ori)
		}
	}
	if(length(ori) == 0){
		ori = data[[1]]
		ori = na.omit(ori)
		data[[1]] = NULL
	}
	return(list(data = data, ori = ori))
}

draw_time_series = function(..., color = c("black", "red", "blue", "green", "pink")){
	data = list(...)
	t = find_ori(data)
	data = t$data
	ori = t$ori
	n = length(ori)
	library(ggplot2)
	g = ggplot()
	g = g + geom_line(aes(x = 1:n, y = ori, colour = "Original data")) + geom_point(aes(x = 1:n, y = ori, colour = "Original data"))
	accu = 0
	breaks = c("Original data")
	for(i in data){
		if(!is.null(i)){
			accu = accu + 1
			sign = ""
			v = vector()
			if(class(i) == "FORECAST"){
				v = i$solution["Forecast", ]
				sign = summary(i, newline = TRUE)$type_of_forecasting
			}else if(class(i) == "CMA"){
				v = i$average
				sign = summary(i, newline = TRUE)
			}else{
				v = i
				sign = paste("No.", accu, "data")
			}
			breaks = append(breaks, sign)
			v = na.omit(v)
			start = n - length(v) + 1
			here = 1
			if(class(i) == "FORECAST"){
				here = i$period
			}else if(class(i) == "CMA"){
				here = sum(!is.na(x)) - length(ori)
			}
			a = aes_string(x = (start + here):(n + here), y = v, colour = shQuote(sign)) #!!!aes_string instead of aes!!!!
			g = g + geom_line(a) + geom_point(a)
		}
	}
	if(!is.null(color)){
		values = rep_len(color, accu + 1)
		names(values) = breaks
		g = g + scale_colour_manual("Colour:", breaks = breaks, values = values)
	}
	g = g + xlab("time") + ylab("value")
	return(g)
}

make_regression = function(lm, var, range, ..., period = 1){
	data = data.frame(...)
	if(length(data) > 0){
		l = length(range) - 1
		v = data[1, ]
		for(i in 1:l){
			data = rbind(data, v)
		}
		data[[var]] = range
	}else{
		eval(parse(text = paste("data = data.frame(", shQuote(var), " = range)")))
	}
	y = predict(lm, data)
	x = vector()
	for(i in 1:(length(range) - period)){
		index = which(lm$model[[var]] == range[i])
		if(length(index) == 0){
			break
		}
		x = append(x, lm$model[index[1], 1])
	}
	if(length(x) == length(range) - period){
		sol = matrix(ncol = length(y), nrow = 2, dimnames = list(c("x", "Forecast"), c(1:(length(y) - period), rep_len("forecast", period))))
		sol["x", ] = c(x, rep_len(NA, period))
		sol["Forecast", ] = y
	}else{
		sol = matrix(y, ncol = length(y), nrow = 1, dimnames = list(c("Forecast"), c(1:(length(y) - period), rep_len("forecast", period))))
	}
	sol = list(solution = as.table(sol), parameter = lm, sign = "Linear Regression", period = period)
	class(sol) = "FORECAST"
	return(sol)
}

make_package = function(y, x = NULL, period = 1){
	if(is.null(x)){
		sol = matrix(y, ncol = length(y), nrow = 1, dimnames = list(c("Forecast"), c(1:(length(y) - period), rep_len("forecast", period))))
	}else{
		y_cover = length(y) - 1
		n = length(x)
		ncol = n
		start = 1
		x = c(x, NA)
		if(y_cover > ncol){
			ncol = y_cover
			start = start - (y_cover - n)
			x = c(rep_len(NA, y_cover - n), x)
		}else{
			y = c(rep_len(NA, n - y_cover), y)
		}
		sol = matrix(c(x, y), byrow = TRUE, ncol = ncol + 1, nrow = 2, dimnames = list(c("x", "Forecast"), c(start:(start + ncol - period), rep_len("forecast", period))))
	}
	sol = list(solution = as.table(sol), parameter = NULL, sign = "User-defined package", period = period)
	class(sol) = "FORECAST"
	return(sol)
}

select_data = function(x, period = NULL, ori = NULL){
	data = NULL
	if(is.null(ori)){
		if(class(x) != "FORECAST" || sum(dimnames(x$solution)[[1]] == "x") == 0){
			return(NULL)
		}
		ori = x$solution["x", ]
	}
	if(class(x) != "FORECAST"){
		data = x
	}else{
		data = x$solution["Forecast", ]
	}
	data = subset(data, names(data) != "forecast")
	ori = na.omit(ori)
	data = na.omit(data)
	if(is.null(period) || period > length(ori)){
		period = length(ori)
	}
	if(period > length(data)){
		period = length(data)
	}
	ori = rev(ori)
	data = rev(data)
	ori = ori[1:period]
	data = data[1:period]
	return(list(ori = ori, data = data, period = period))
}

mean_error = function(x, period = NULL, ori = NULL){
	t = select_data(x, period, ori)
	if(is.null(t)){
		return(NULL)
	}
	return(mean(t$ori - t$data))
}

mean_percentage_error = function(x, period = NULL, ori = NULL){
	t = select_data(x, period, ori)
	if(is.null(t)){
		return(NULL)
	}
	return(mean((t$ori - t$data) / t$ori))
}

mean_absolute_deviation = function(x, period = NULL, ori = NULL){
	t = select_data(x, period, ori)
	if(is.null(t)){
		return(NULL)
	}
	return(mean(abs(t$ori - t$data)))
}

mean_square_error = function(x, period = NULL, ori = NULL){
	t = select_data(x, period, ori)
	if(is.null(t)){
		return(NULL)
	}
	return(mean((t$ori - t$data) ^ 2))
}

mean_absolute_percentage_error = function(x, period = NULL, ori = NULL){
	t = select_data(x, period, ori)
	if(is.null(t)){
		return(NULL)
	}
	return(mean(abs(t$ori - t$data) / t$ori))
}

compare_forecast = function(...){
	data = list(...)
	t = find_ori(data)
	data = t$data
	ori = t$ori
	accu = 0
	l = list()
	for(i in data){
		if(!is.null(i)){
			accu = accu + 1
			sign = ""
			if(class(i) == "FORECAST"){
				x = i$solution["Forecast", ]
				x = subset(x, names(x) != "forecast")
				l[[accu]] = x
				sign = summary(i)$type_of_forecasting
			}else{
				l[[accu]] = i
				l[[accu]] = l[[accu]][-length(l[[accu]])]
				sign = paste("No.", accu, "data")
			}
			names(l)[accu] = sign
		}
	}
	len = length(na.omit(l[[1]]))
	for(i in 2:accu){
		if(len > length(na.omit(l[[i]]))){
			len = length(na.omit(l[[i]]))
		}
	}
	if(len > length(ori)){
		len = length(ori)
	}
	sol = matrix(nrow = length(l), ncol = 5, dimnames = list(names(l), c("ME", "MPE", "MAD", "MSE", "MAPE")))
	sol[, 1] = sapply(l, mean_error,                     period = len, ori = ori)
	sol[, 2] = sapply(l, mean_percentage_error,          period = len, ori = ori)
	sol[, 3] = sapply(l, mean_absolute_deviation,        period = len, ori = ori)
	sol[, 4] = sapply(l, mean_square_error,              period = len, ori = ori)
	sol[, 5] = sapply(l, mean_absolute_percentage_error, period = len, ori = ori)
	return(as.table(sol))
}

centered_moving_average = function(x, period){
	x = data_filt(x)
	sol = list()
	x = na.omit(x)
	v = vector(length = length(x) + 1 - period)
	start = 0
	while(start + period <= length(x)){
		v[start + 1] = mean(x[(start + 1):(start + period)])
		start = start + 1
	}
	if(period %% 2 == 0){
		tmp = c(v, 0)
		tmp = tmp[-1]
		v = (v + tmp) / 2
		v = v[-length(v)]
	}
	sol[["x"]] = x
	blank = (length(x) - length(v)) / 2
	sol[["average"]] = c(rep_len(NA, blank), v, rep_len(NA, blank))
	names(sol[["average"]]) = 1:length(x)
	sol[["period"]] = period
	class(sol) = "CMA"
	return(sol)
}

print.CMA = function(x){
	print(as.table(x$average))
}

mean.CMA = function(x){
	return(mean(x$x, na.rm = TRUE))
}

summary.CMA = function(x, newline = FALSE, ...){
	sol = ""
	if(newline){
		sol = paste("CMA for\nperiod =", x$period)
	}else{
		sol = paste("CMA for period =", x$period)
	}
	class(sol) = "summary.CMA"
	return(sol)
}

seasonal_index = function(x, ...){
	UseMethod("seasonal_index")
}

seasonal_index.CMA = function(x, ...){
	effect = x$x / x$average
	tmp = 1:length(x$x)
	ind = vector(length = x$period)
	for(i in 1:x$period){
		y = effect[tmp %% x$period == i %% x$period]
		ind[i] = mean(y, na.rm = TRUE)
	}
	ind = ind / sum(ind) * x$period
	return(ind)
}

seasonal_index.lm = function(x, period, ...){
	if(class(x) != "lm"){
		x = lm(x~c(1:length(x)))
	}
	predict = predict(x)
	effect = x$model[[1]] / predict
	tmp = 1:length(predict)
	ind = vector(length = period)
	for(i in 1:period){
		y = effect[tmp %% period == i %% period]
		ind[i] = mean(y, na.rm = TRUE)
	}
	ind = ind / sum(ind) * period
	return(ind)
}

seasonal_effect = function(x, season.period, index = c("CMA", "lm"), forecast.way = c("Holt\'s", "lm", "average", "moving", "exponential"), ...){
	x = data_filt(x)
	index = match.arg(index)
	forecast.way = match.arg(forecast.way)
	y = 0
	if(index == "CMA"){
		y = centered_moving_average(x, season.period)
	}else if(index == "lm"){
		y = lm(x~c(1:length(x)))
	}
	season = seasonal_index(y, period = season.period)
	x.season_dropped = x / season
	o = 0
	f = NULL
	f.e = NULL
	if(forecast.way == "lm"){ #use of FORECASTed lm object is somehow BAD choice
		tmp = 1:length(x)
		r = lm(x.season_dropped ~ tmp)
		df = data.frame(tmp = 1 : (length(x) + season.period))
		o = predict(r, df)
		f.e = list(coe = r$coefficients, season = season)
		names(f.e$coe) = c("intercept", "x")
		f = function(x, period, expression){
			n = length(na.omit(x$solution["x", ]))
			return((expression$coe["x"] * (n + period) + expression$coe["intercept"]) * ifelse((n + period) %% x$period == 0, expression$season[x$period], expression$season[(n + period) %% x$period]))
		}
	}else{
		obj = NULL
		if(forecast.way == "Holt\'s"){
			obj = holt_exponential_smoothing(x, season.period, ...)
		}else if(forecast.way == "average"){
			obj = average_value(x)
		}else if(forecast.way == "moving"){
			obj = moving_average(x, ...)
		}else if(forecast.way == "exponential"){
			obj = exponential_smoothing(x, ...)
		}
		f.e = obj
		f = function(x, period, expression){
			return(forecast(expression, period) * ifelse((n + period) %% x$period == 0, expression$season[x$period], expression$season((n + period) %% x$period)))
		}
		o = obj$solution["Forecast", ]
		o = o[which(names(o) == "1"):length(o)]
		o = subset(o, names(o) != "forecast")
		o = append(o, mapply(forecast, x = obj, period = 1:season.period))
	}
	o_ori = o
	o = o * season
	sol = matrix(ncol = length(o), nrow = 5, dimnames = list(c("x", "Forecast", "deseasonal", "Forecast.deseasonal", "seasonal.index"), c(1:(length(o) - season.period), rep_len("forecast", season.period))))
	sol["x", ] = c(x, rep_len(NA, season.period))
	sol["Forecast", ] = o
	sol["deseasonal", ] = c(x.season_dropped, rep_len(NA, season.period))
	sol["Forecast.deseasonal", ] = o_ori
	sol["seasonal.index", ] = rep_len(season, length(o))
	sol = list(solution = as.table(sol), parameter = c(index = index, forecast.way = forecast.way), sign = "Seasonal", period = season.period, forecast_function = f, forecast_function.expression = f.e)
	class(sol) = "FORECAST"
	return(sol)
}
