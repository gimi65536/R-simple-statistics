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
	}
	names(sol)[3] = "type_of_forecasting"
	class(sol) = "summary.FORECAST"
	return(sol[-2])
}

average_value = function(x){
	sol = matrix(ncol = length(x) + 1, nrow = 2, dimnames = list(c("x", "Forecast"), c(1:length(x), "forecast")))
	sol["x", ] = c(x, NA)
	accu = x[1]
	for(i in 2:length(x)){
		sol["Forecast", i] = accu / (i - 1)
		accu = accu + x[i]
	}
	sol["Forecast", length(x) + 1] = accu / length(x)
	sol = list(solution = as.table(sol), parameter = NULL, sign = "Average")
	class(sol) = "FORECAST"
	return(sol)
}

moving_average = function(x, period = 3){
	if(length(x) <= period){
		return (NULL)
	}
	sol = matrix(ncol = length(x) + 1, nrow = 2, dimnames = list(c("x", "Forecast"), c(1:length(x), "forecast")))
	sol["x", ] = c(x, NA)
	for(i in 1:(length(x) - period + 1)){
		sol["Forecast", period + i] = sum(x[i:(i + period - 1)]) / period
	}
	sol = list(solution = as.table(sol), parameter = c(period = period), sign = "Moving")
	class(sol) = "FORECAST"
	return(sol)
}

exponential_smoothing = function(x, alpha = 0.2){
	sol = matrix(ncol = length(x) + 1, nrow = 2, dimnames = list(c("x", "Forecast"), c(1:length(x), "forecast")))
	sol["x", ] = c(x, NA)
	sol["Forecast", 1] = x[1]
	for(i in 2:(length(x) + 1)){
		sol["Forecast", i] = alpha * x[i - 1] + (1 - alpha) * sol["Forecast", i - 1]
	}
	sol = list(solution = as.table(sol), parameter = c(alpha = alpha), sign = "Exponent")
	class(sol) = "FORECAST"
	return(sol)
}

holt_exponential_smoothing = function(x, alpha = 0.9, beta = 0.8, ET0 = 1, F0 = mean(x), X0 = F0){
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
	sol = list(solution = as.table(sol), parameter = c(alpha = alpha, beta = beta, ET0 = ET0, F0 = F0, X0 = X0), sign = "Holt\'s")
	class(sol) = "FORECAST"
	return(sol)
}

holt_forecast = function(x, period){
	if(class(x) != "FORECAST" || x$sign != "Holt\'s" || period <= 0){
		return(NULL)
	}
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
		}
	}
	if(length(ori) == 0){
		ori = data[[1]]
		ori = na.omit(ori)
		data[[1]] = NULL
	}
	return(list(data = data, ori = ori))
}

draw_time_series = function(data, color = c("black", "red", "blue", "green", "pink"), ...){
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
			}else{
				v = i
				sign = paste("No.", accu, "data")
			}
			breaks = append(breaks, sign)
			v = na.omit(v)
			start = n - length(v) + 2
			a = aes_string(x = start:(n + 1), y = v, colour = shQuote(sign)) #!!!aes_string instead of aes!!!!
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

make_regression = function(lm, var, range, ...){
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
	for(i in 1:(length(range) - 1)){
		index = which(lm$model[[var]] == range[i])
		if(length(index) == 0){
			break
		}
		x = append(x, lm$model[index[1], 1])
	}
	if(length(x) == length(range) - 1){
		sol = matrix(ncol = length(y), nrow = 2, dimnames = list(c("x", "Forecast"), c(1:(length(y) - 1), "forecast")))
		sol["x", ] = c(x, NA)
		sol["Forecast", ] = y
	}else{
		sol = matrix(y, ncol = length(y), nrow = 1, dimnames = list(c("Forecast"), c(1:(length(y) - 1), "forecast")))
	}
	sol = list(solution = as.table(sol), parameter = c(lm.call = as.character(lm$call)[2]), sign = "Linear Regression")
	class(sol) = "FORECAST"
	return(sol)
}

make_package = function(y, x = NULL){
	if(is.null(x)){
		sol = matrix(y, ncol = length(y), nrow = 1, dimnames = list(c("Forecast"), c(1:(length(y) - 1), "forecast")))
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
		sol = matrix(c(x, y), byrow = TRUE, ncol = ncol + 1, nrow = 2, dimnames = list(c("x", "Forecast"), c(start:(start + ncol - 1), "forecast")))
	}
	sol = list(solution = as.table(sol), parameter = NULL, sign = "User-defined package")
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
	ori = na.omit(ori)
	data = na.omit(data)
	data = data[-length(data)]
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

compare_forecast = function(data){
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
				l[[accu]] = i$solution["Forecast", ]
				sign = summary(i)$type_of_forecasting
			}else{
				l[[accu]] = i
				sign = paste("No.", accu, "data")
			}
			names(l)[accu] = sign
		}
	}
	len = length(l[[1]])
	for(i in 2:accu){
		if(len > length(na.omit(l[[i]]))){
			len = length(na.omit(l[[i]]))
		}
	}
	if(len > length(x) + 1){
		len = length(x) + 1
	}
	sol = matrix(nrow = length(l), ncol = 5, dimnames = list(names(l), c("ME", "MPE", "MAD", "MSE", "MAPE")))
	sol[, 1] = sapply(l, mean_error, period = len - 1, ori = ori)
	sol[, 2] = sapply(l, mean_percentage_error, period = len - 1, ori = ori)
	sol[, 3] = sapply(l, mean_absolute_deviation, period = len - 1, ori = ori)
	sol[, 4] = sapply(l, mean_square_error, period = len - 1, ori = ori)
	sol[, 5] = sapply(l, mean_absolute_percentage_error, period = len - 1, ori = ori)
	return(as.table(sol))
}
