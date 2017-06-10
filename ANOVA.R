getName = function(df, group = "Group", space = " "){
	sol = names(df)
	if(is.null(sol)){
		sol = paste(group, 1:length(df), sep = space)
	}else{
		for(i in 1:length(df)){
			if(is.na(sol[i])){
				sol[i] = group + space + i
			}
		}
	}
	return (sol)
}

ANOVA.oneway = function(..., buildin = FALSE, table = FALSE, df = NULL, index = NULL){
	if(is.null(df)){
		args = list(...)
		if(length(index) != 0){
			args = args[index]
		}
		bartp = BartlettTest(..., index = index)[[2]][[2]]
	}else{
		args = as.list(df)
		if(length(index) != 0){
			args = args[index]
		}
		bartp = BartlettTest(df = df)[[2]][[2]]
	}
	name = getName(args)
	number.treatment = length(args)
	if(number.treatment == 0){
		return (NULL)
	}
	if(buildin){
		vec = vector()
		fac = vector()
		for(i in 1:number.treatment){
			x = args[[i]]
			x = x[!is.na(x)]
			vec = c(vec, x)
			fac = c(fac, rep(i, length(x)))
		}
		fac = as.factor(fac)
		df = data.frame(y = vec)
		return(aov(y ~ fac, data = df))
	}
	number.element = vector(mode = "integer", length = number.treatment)
	mean.treatment = vector(mode = "numeric", length = number.treatment)
	var.treatment = vector(mode = "numeric", length = number.treatment)
	shap = vector()
	for(i in 1:number.treatment){
		x = args[[i]]
		x = x[!is.na(x)]
		number.element[i] = length(x)
		mean.treatment[i] = mean(x)
		var.treatment[i] = var(x)
		shap = append(shap, ShapiroTest(args[[i]])[[2]])
	}
	number = sum(number.element)
	mean = sum(number.element * mean.treatment) / number
	SST = sum(number.element * (mean.treatment - mean) ^ 2)
	SSE = sum((number.element - 1) * var.treatment)
	MST = SST / (number.treatment - 1)
	MSE = SSE / (number - number.treatment)
	F = MST / MSE
	p = 1 - pf(F, number.treatment - 1, number - number.treatment)
	if(!table){
		t1 = list(degeree_of_freedom = number.treatment - 1, sum_of_sqares = SST, mean_squares = MST)
		t2 = list(degeree_of_freedom = number - number.treatment, sum_of_sqares = SSE, mean_squares = MSE)
		t3 = list(degeree_of_freedom = number - 1, sum_of_sqares = SST + SSE)
		t4 = list(F_statistic = F, p_value = p)
		t5 = list(number_of_treatment = number.treatment, number_of_all_elements = number, total_mean = mean, number_of_elements_in_treatments = number.element,
							mean_of_each_treatment = mean.treatment, var_of_each_treatment = var.treatment, Shapiro_pvalue = shap, Bartlett_pvalue = bartp)
		sol = list(Treatments = t1, Error = t2, Total = t3, F_statistic = t4, Additional = t5)
	}else{
		x1 = c(number.treatment - 1, SST, MST, F, number - number.treatment, SSE, MSE, NA, number - 1, SST + SSE, NA, NA)
		n1 = list(c("Treatment", "Error", "Total"), c("Degree of Freedom", "Sum of Squares", "Mean Squares", "F-statistic"))
		mat1 = matrix(data = x1, nrow = 3, ncol = 4, dimnames = n1, byrow = TRUE)
		t1 = as.table(mat1)
		x2 = c(number.element, number, mean.treatment, mean, var.treatment, NA)
		n2 = list(c("number", "mean", "variance"), c(name, "Total"))
		mat2 = matrix(data = x2, nrow = 3, ncol = number.treatment + 1, dimnames = n2, byrow = TRUE)
		t2 = as.table(mat2)
		mat3 = matrix(data = p, nrow = 1, dimnames = list("pvalue", ""))
		t3 = as.table(mat3)
		mat4 = matrix(data = shap, nrow = 1, dimnames = list("pvalue", name))
		t4 = as.table(mat4)
		mat5 = matrix(data = bartp, nrow = 1, dimnames = list("pvalue", ""))
		t5 = as.table(mat5)
		sol = list(ANOVA_table = t1, Additional_information = t2, p_value = t3, Shapiro_pvalue = t4, Bartlett_pvalue = t5)
	}
	return (sol)
}

ShapiroTest = function(x){ #check if the vector is normally distributed
	return(shapiro.test(x))
}

ShapiroTestALL = function(..., df = NULL, index = NULL){
	if(is.null(df)){
		args = list(...)
	}else{
		args = as.list(df)
	}
	if(length(index) != 0){
		args = args[index]
	}
	name = getName(args)
	n = length(args)
	if(n == 0){
		return (NULL)
	}
	w = vector(mode = "numeric", length = n)
	p = vector(mode = "numeric", length = n)
	for(i in 1:n){
		s = shapiro.test(args[[i]])
		w[i] = s[[1]]
		p[i] = s[[2]]
	}
	x = c(w, p)
	nama = list(c("W", "pvalue"), name)
	mat = matrix(data = x, nrow = 2, ncol = n, dimnames = nama, byrow = TRUE)
	return (as.table(mat))
}

BartlettTest = function(..., buildin = FALSE, table = FALSE, df = NULL, index = NULL){
	if(is.null(df)){
		args = list(...)
	}else{
		args = as.list(df)
	}
	if(length(index) != 0){
		args = args[index]
	}
	name = getName(args)
	number.treatment = length(args)
	if(number.treatment == 0){
		return (NULL)
	}
	if(buildin){
		vec = vector()
		fac = vector()
		for(i in 1:number.treatment){
			x = args[[i]]
			x = x[!is.na(x)]
			vec = c(vec, x)
			fac = c(fac, rep(i, length(x)))
		}
		fac = as.factor(fac)
		df = data.frame(y = vec)
		return(bartlett.test(y ~ fac, data = df))
	}
	number.element = vector(mode = "integer", length = number.treatment)
	var.treatment = vector(mode = "numeric", length = number.treatment)
	for(i in 1:number.treatment){
		x = args[[i]]
		x = x[!is.na(x)]
		number.element[i] = length(x)
		var.treatment[i] = var(x)
	}
	number = sum(number.element)
	SSE = sum((number.element - 1) * var.treatment)
	MSE = SSE / (number - number.treatment)
	C = 1 + 1 / (3 * (number.treatment - 1)) * (sum(1 / (number.element - 1)) - 1 / (number - number.treatment))
	B = ((number - number.treatment) * log(MSE) - sum((number.element - 1) * log(var.treatment))) / C
	p = 1 - pchisq(B, number.treatment - 1)
	if(!table){
		t1 = list(sum_of_sqares = SSE, mean_squares = MSE)
		t2 = list(B_statistic = B, p_value = p, degeree_of_freedom = number.treatment - 1)
		t3 = list(number_of_treatment = number.treatment, number_of_all_elements = number, number_of_elements_in_treatments = number.element,
							var_of_each_treatment = var.treatment, C = C)
		sol = list(Error = t1, B_statistic = t2, Additional = t3)
	}else{
		x1 = c(SSE, MSE, B, C, number.treatment - 1, p)
		n1 = list(c("SSE", "MSE", "B statistic", "C value", "Degree of Freedom", "pvalue"), "")
		mat1 = matrix(data = x1, ncol = 1, dimnames = n1)
		t1 = as.table(mat1)
		x2 = c(number.element, number, var.treatment, NA)
		n2 = list(c("number", "variance"), c(name, "Total"))
		mat2 = matrix(data = x2, nrow = 2, ncol = number.treatment + 1, dimnames = n2, byrow = TRUE)
		t2 = as.table(mat2)
		sol = list(Solution = t1, Additional_information = t2)
	}
	return (sol)
}

cross_differ = function(..., alpha = 0.05, way = "LSD", df = NULL, label = FALSE, index = NULL){
	if(is.null(df)){
		data = ANOVA.oneway(..., table = FALSE, index = index)
		name = getName(list(...))
	}else{
		if(label){
			data = ANOVA.oneway(df = df[-1,], table = FALSE, index = index)
			name = as.character(as.vector(df[1,]))
		}else{
			data = ANOVA.oneway(df = df, table = FALSE, index = index)
			name = getName(df)
		}
	}
	if(length(index) != 0){
		name = name[index]
	}
	number.treatment = data[[5]][[1]]
	number.element = data[[5]][[4]]
	number = data[[5]][[2]]
	mean.treatment = data[[5]][[5]]
	MSE = data[[2]][[3]]
	nama = list(name, name)
	sol = matrix(nrow = number.treatment, ncol = number.treatment, dimnames = nama)
	inf = matrix(nrow = number.treatment, ncol = number.treatment, dimnames = nama)
	cri = matrix(nrow = number.treatment, ncol = number.treatment, dimnames = nama)
	if(way == "Tukey"){
		n = 1 / mean(1 / number.element)
		q = qtukey(1 - alpha, number.treatment, number - number.treatment)
		critical = q * sqrt(MSE / n)
		for(i in 1 : (number.treatment - 1)){
			for(j in (i + 1) : number.treatment){
				sol[j, i] = (abs(mean.treatment[i] - mean.treatment[j]) > critical)
				inf[j, i] = abs(mean.treatment[i] - mean.treatment[j])
				cri[j, i] = critical
			}
		}
	}else{
		Alpha = alpha
		if(way == "Bonferroni"){
			Alpha = 2 * Alpha / number.treatment / (number.treatment - 1)
		}
		t = qt(1 - Alpha / 2, number - number.treatment)
		for(i in 1 : (number.treatment - 1)){
			for(j in (i + 1) : number.treatment){
				critical = t * sqrt(MSE * (1 / number.element[i] + 1 / number.element[j]))
				sol[j, i] = (abs(mean.treatment[i] - mean.treatment[j]) > critical)
				inf[j, i] = abs(mean.treatment[i] - mean.treatment[j])
				cri[j, i] = critical
			}
		}
	}
	return (list(solution = as.table(sol), difference = as.table(inf), critical_value = as.table(cri)))
}

ANOVA.twoway = function(..., buildin = FALSE, table = FALSE, df = NULL, label = FALSE, index = NULL){
	blockName = NULL
	if(is.null(df)){
		args = list(...)
	}else{
		if(label){
			args = as.list(subset(df, select = -1))
			blockName = as.character(df[[1]])
		}else{
			args = as.list(df)
		}
	}
	if(length(index) != 0){
		args = args[index]
	}
	treatmentnName = getName(args, group = "Treatment")
	number.treatment = length(args)
	if(number.treatment == 0){
		return (NULL)
	}
	number.element = -1
	for(i in 1:number.treatment){
		x = args[[i]]
		x = x[!is.na(x)]
		if(number.element == -1){
			number.element = length(x)
		}else if(length(x) != number.element){
			return (NULL)
		}
	}
	number.block = number.element
	if(is.null(blockName)){
		blockName = paste("Block", 1:number.block);
	}
	if(buildin){
		vec = vector()
		fac = vector()
		for(i in 1:number.treatment){
			x = args[[i]]
			x = x[!is.na(x)]
			vec = c(vec, x)
			fac = c(fac, rep(i, length(x)))
		}
		fac = as.factor(fac)
		fac2 = rep(1:number.element, number.treatment)
		fac2 = as.factor(fac2)
		df = data.frame(y = vec)
		return(aov(y ~ fac + fac2, data = df))
	}
	#number.element = vector(mode = "integer", length = number.treatment)
	mean.treatment = vector(mode = "numeric", length = number.treatment)
	mean.block = vector(mode = "numeric", length = number.block)
	var.treatment = vector(mode = "numeric", length = number.treatment)
	var.block = vector(mode = "numeric", length = number.block)
	data = matrix(nrow = number.block, ncol = number.treatment)
	for(i in 1:number.treatment){
		x = args[[i]]
		x = x[!is.na(x)]
		data[, i] = x
		#number.element[i] = length(x)
		mean.treatment[i] = mean(x)
		var.treatment[i] = var(x)
	}
	for(i in 1:number.block){
		mean.block[i] = mean(data[i, ])
		var.block[i] = var(data[i, ])
	}
	number = number.treatment * number.block
	mean = sum(number.element * mean.treatment) / number
	SST = sum(number.element * (mean.treatment - mean) ^ 2)
	SSB = sum(number.treatment * (mean.block - mean) ^ 2)
	SSE = 0
	for(i in 1:number.treatment){
		for(j in 1:number.block){
			SSE = SSE + (data[j, i] - mean.treatment[i] - mean.block[j] + mean) ^ 2
		}
	}
	MST = SST / (number.treatment - 1)
	MSB = SSB / (number.block - 1)
	MSE = SSE / (number - number.treatment - number.block + 1)
	F = MST / MSE
	FB = MSB / MSE
	p = 1 - pf(F, number.treatment - 1, number - number.treatment - number.block + 1)
	pB = 1 - pf(FB, number.block - 1, number - number.treatment - number.block + 1)
	if(!table){
		t1 = list(degeree_of_freedom = number.treatment - 1, sum_of_sqares = SST, mean_squares = MST)
		t15 = list(degeree_of_freedom = number.block - 1, sum_of_sqares = SSB, mean_squares = MSB)
		t2 = list(degeree_of_freedom = number - number.treatment - number.block + 1, sum_of_sqares = SSE, mean_squares = MSE)
		t3 = list(degeree_of_freedom = number - 1, sum_of_sqares = SST + SSB + SSE)
		t4 = list(F_statistic = F, p_value = p)
		t45 = list(F_statistic = FB, p_value = pB)
		t5 = list(number_of_treatment = number.treatment, number_of_all_elements = number, total_mean = mean, number_of_elements_in_treatments = number.element,
							mean_of_each_treatment = mean.treatment, var_of_each_treatment = var.treatment, mean_of_each_block = mean.block, var_of_each_block = var.block)
		sol = list(Treatments = t1, Blocks = t15, Error = t2, Total = t3, F_statistic_for_Treatments = t4, F_statistic_for_Blocks = t45, Additional = t5)
	}else{
		x1 = c(number.treatment - 1, SST, MST, F, number.block - 1, SSB, MSB, FB, number - number.treatment - number.block + 1, SSE, MSE, NA, number - 1, SST + SSB + SSE, NA, NA)
		n1 = list(c("Treatment", "Block", "Error", "Total"), c("Degree of Freedom", "Sum of Squares", "Mean Squares", "F-statistic"))
		mat1 = matrix(data = x1, nrow = 4, ncol = 4, dimnames = n1, byrow = TRUE)
		t1 = as.table(mat1)
		x2 = c(rep(number.element, number.treatment), number, mean.treatment, mean, var.treatment, NA)
		n2 = list(c("number", "mean", "variance"), c(treatmentnName, "Total"))
		mat2 = matrix(data = x2, nrow = 3, ncol = number.treatment + 1, dimnames = n2, byrow = TRUE)
		t2 = as.table(mat2)
		mat3 = matrix(data = p, nrow = 1, dimnames = list("pvalue", ""))
		t3 = as.table(mat3)
		x4 = c(rep(number.treatment, number.block), number, mean.block, mean, var.block, NA)
		n4 = list(c(blockName, "Total"), c("number", "mean", "variance"))
		mat4 = matrix(data = x4, nrow = number.block + 1, ncol = 3, dimnames = n4, byrow = FALSE)
		t4 = as.table(mat4)
		mat5 = matrix(data = pB, nrow = 1, dimnames = list("pvalue", ""))
		t5 = as.table(mat5)
		sol = list(ANOVA_table = t1, Additional_information_for_Treatments = t2, p_value_for_Treatments = t3, Additional_information_for_Blocks = t4, p_value_for_Blocks = t5)
	}
	return (sol)
}

cross_differ_twoway = function(..., alpha = 0.05, way = "LSD", df = NULL, label = FALSE, index = NULL){
	blockName = NULL
	if(is.null(df)){
		data = ANOVA.twoway(..., table = FALSE, index = index)
		treatmentnName = getName(list(...))
	}else{
		if(label){
			data = ANOVA.twoway(df = subset(df, select = -1), table = FALSE, index = index)
			treatmentnName = getName(subset(df, select = -1))
			blockName = as.character(df[[1]])
		}else{
			data = ANOVA.twoway(df = df, table = FALSE, index = index)
			treatmentnName = getName(df)
		}
	}
	if(length(index) != 0){
		treatmentnName = treatmentnName[index]
	}
	number.treatment = data[[7]][[1]]
	number.element = data[[7]][[4]]
	number.block = number.element
	number = data[[7]][[2]]
	mean.treatment = data[[7]][[5]]
	mean.block = data[[7]][[7]]
	MSE = data[[3]][[3]]
	if(is.null(blockName)){
		blockName = paste("Block", 1:number.block)
	}
	nama = list(treatmentnName, treatmentnName)
	sol = matrix(nrow = number.treatment, ncol = number.treatment, dimnames = nama)
	inf = matrix(nrow = number.treatment, ncol = number.treatment, dimnames = nama)
	cri = matrix(nrow = number.treatment, ncol = number.treatment, dimnames = nama)
	nama2 = list(blockName, blockName)
	sol2 = matrix(nrow = number.block, ncol = number.block, dimnames = nama2)
	inf2 = matrix(nrow = number.block, ncol = number.block, dimnames = nama2)
	cri2 = matrix(nrow = number.block, ncol = number.block, dimnames = nama2)
	if(way == "Tukey"){
		n = number.element
		q = qtukey(1 - alpha, number.treatment, number - number.treatment - number.block + 1)
		critical = q * sqrt(MSE / n)
		for(i in 1 : (number.treatment - 1)){
			for(j in (i + 1) : number.treatment){
				sol[j, i] = (abs(mean.treatment[i] - mean.treatment[j]) > critical)
				inf[j, i] = abs(mean.treatment[i] - mean.treatment[j])
				cri[j, i] = critical
			}
		}
		n2 = number.treatment
		q2 = qtukey(1 - alpha, number.block, number - number.treatment - number.block + 1)
		critical2 = q2 * sqrt(MSE / n2)
		for(i in 1 : (number.block - 1)){
			for(j in (i + 1) : number.block){
				sol2[j, i] = (abs(mean.block[i] - mean.block[j]) > critical)
				inf2[j, i] = abs(mean.block[i] - mean.block[j])
				cri2[j, i] = critical2
			}
		}
	}else{
		Alpha = alpha
		if(way == "Bonferroni"){
			Alpha = 2 * Alpha / number.treatment / (number.treatment - 1)
		}
		t = qt(1 - Alpha / 2, number - number.treatment - number.block + 1)
		for(i in 1 : (number.treatment - 1)){
			for(j in (i + 1) : number.treatment){
				critical = t * sqrt(MSE * (1 / number.element + 1 / number.element))
				sol[j, i] = (abs(mean.treatment[i] - mean.treatment[j]) > critical)
				inf[j, i] = abs(mean.treatment[i] - mean.treatment[j])
				cri[j, i] = critical
			}
		}
		Alpha = alpha
		if(way == "Bonferroni"){
			Alpha = 2 * Alpha / number.block / (number.block - 1)
		}
		t2 = qt(1 - Alpha / 2, number - number.treatment - number.block + 1)
		for(i in 1 : (number.block - 1)){
			for(j in (i + 1) : number.block){
				critical = t2 * sqrt(MSE * (1 / number.treatment + 1 / number.treatment))
				sol2[j, i] = (abs(mean.block[i] - mean.block[j]) > critical)
				inf2[j, i] = abs(mean.block[i] - mean.block[j])
				cri2[j, i] = critical
			}
		}
	}
	return (list(for_treatment = list(solution = as.table(sol), difference = as.table(inf), critical_value = as.table(cri)),
		for_block = list(solution = as.table(sol2), difference = as.table(inf2), critical_value = as.table(cri2))))
}

ANOVA.2factor = function(df, buildin = FALSE, table = FALSE, balance = TRUE, index = NULL){
	if(length(index) != 0){
		df = df[index]
	}
	f = as.factor(df[[1]])
	level = levels(f)
	for(i in 2 : length(f)){
		if(is.na(f[i])){
			f[i] = f[i - 1]
		}
	}
	df = subset(df, select = -1)
	name = getName(df)
	number.factor = c(length(df), length(level))
	data = vector()
	f1 = vector()
	f2 = vector()
	for(i in 1 : number.factor[1]){
		data = append(data, df[[i]])
		f1 = append(f1, rep_len(name[i], length(f)))
		f2 = append(f2, as.character(f))
	}
	be_deleted = union(which(is.na(data)), which(is.na(f2)))
	if(length(be_deleted) > 0){
		data = data[-be_deleted]
		f1 = f1[-be_deleted]
		f2 = f2[-be_deleted]
	}
	if(buildin){
		df = data.frame(y = data)
		return(aov(y ~ f1 + f1 * f2 + f2, data = df))
	}
	pnor = ShapiroTest(data)[[2]]
	nama.intersect = list(level, name)
	number.intersect = number.factor[1] * number.factor[2]
	number.element = matrix(nrow = 2, ncol = max(number.factor))
	mean.treatment = matrix(nrow = 2, ncol = max(number.factor))
	var.treatment = matrix(nrow = 2, ncol = max(number.factor))
	number.element.intersect = matrix(ncol = number.factor[1], nrow = number.factor[2], dimnames = nama.intersect)
	mean.intersect = matrix(ncol = number.factor[1], nrow = number.factor[2], dimnames = nama.intersect)
	var.intersect = matrix(ncol = number.factor[1], nrow = number.factor[2], dimnames = nama.intersect)
	rA = rep_len(0, number.factor[1])
	rB = rep_len(0, number.factor[2])
	for(n in 1:2){
		for(i in 1 : number.factor[n]){
			if(n == 1){
				w = which(f1 == name[i])
			}else{
				w = which(f2 == level[i])
			}
			x = data[w]
			number.element[n, i] = length(w)
			if(length(w) > 0){
				mean.treatment[n, i] = mean(x)
				var.treatment[n, i] = var(x)
			}else{
				mean.treatment[n, i] = 0
				var.treatment[n, i] = 0
			}
		}
	}
	datav = list()
	for(i in 1 : number.factor[2]){
		for(j in 1 : number.factor[1]){
			w = intersect(which(f1 == name[j]), which(f2 == level[i]))
			x = data[w]
			number.element.intersect[i, j] = length(w)
			if(length(w) > 0){
				mean.intersect[i, j] = mean(x)
				var.intersect[i, j] = var(x)
			}else{
				mean.intersect[i, j] = 0
				var.intersect[i, j] = 0
			}
			rA[j] = rA[j] + number.element.intersect[i, j]
			rB[i] = rB[i] + number.element.intersect[i, j]
			datav[[length(datav) + 1]] = x
		}
	}
	rA = rA / number.factor[2]
	rB = rB / number.factor[1]
	number = length(data)
	mean = mean(data)
	SST = sum((data - mean) ^ 2)
	SSA = 0
	SSB = 0
	x = mean.treatment[1, ]
	x = x[!is.na(x)]
	for(i in 1 : number.factor[1]){
		f1_ind = which(name == f1[i])
		SSA = SSA + (mean.treatment[1, i] - mean) ^ 2 * rA[f1_ind]
	}
	SSA = SSA * number.factor[2]
	x = mean.treatment[2, ]
	x = x[!is.na(x)]
	for(i in 1 : number.factor[2]){
		f2_ind = which(level == f2[i])
		SSB = SSB + (mean.treatment[2, i] - mean) ^ 2 * rB[f2_ind]
	}
	SSB = SSB * number.factor[1]
	SSAB = 0
	SSE = 0
	for(i in 1 : number.factor[2]){
		for(j in 1 : number.factor[1]){
			SSAB = SSAB + (mean.intersect[i, j] - mean.treatment[1, j] - mean.treatment[2, i] + mean) ^ 2 * number.element.intersect[i, j]
		}
	}
	pvar = BartlettTest(df = as.data.frame(datav))[[2]][[2]]
	for(i in 1 : number){
		f1_ind = which(name == f1[i])
		f2_ind = which(level == f2[i])
		SSE = SSE + (data[i] - mean.intersect[f2_ind, f1_ind]) ^ 2
	}
	MSA = SSA / (number.factor[1] - 1)
	MSB = SSB / (number.factor[2] - 1)
	MSAB = SSAB / ((number.factor[1] - 1) * (number.factor[2] - 1))
	MSE = SSE / (number - number.factor[1] * number.factor[2])
	FA = MSA / MSE
	FB = MSB / MSE
	FAB = MSAB / MSE
	pA = 1 - pf(FA, number.factor[1] - 1, number - number.factor[1] * number.factor[2])
	pB = 1 - pf(FB, number.factor[2] - 1, number - number.factor[1] * number.factor[2])
	pAB = 1 - pf(FAB, (number.factor[1] - 1) * (number.factor[2] - 1), number - number.factor[1] * number.factor[2])
	if(!table){
		t1 = list(degeree_of_freedom = number.factor[1] - 1, sum_of_sqares = SSA, mean_squares = MSA)
		t2 = list(degeree_of_freedom = number.factor[2] - 1, sum_of_sqares = SSB, mean_squares = MSB)
		t3 = list(degeree_of_freedom = (number.factor[1] - 1) * (number.factor[2] - 1), sum_of_sqares = SSAB, mean_squares = MSAB)
		t4 = list(degeree_of_freedom = number - number.factor[1] * number.factor[2], sum_of_sqares = SSE, mean_squares = MSE)
		t5 = list(degeree_of_freedom = number - 1, sum_of_sqares = SSA + SSB + SSAB + SSE)
		t6 = list(F_statistic = FA, p_value = pA)
		t7 = list(F_statistic = FB, p_value = pB)
		t8 = list(F_statistic = FAB, p_value = pAB)
		t9 = list(number_of_factor1 = number.factor[1], number_of_factor2 = number.factor[2], number_of_intersect = number.intersect,
					number_of_elements_in_factor1 = number.element[1, ][!is.na(number.element[1, ])],
					number_of_elements_in_factor2 = number.element[2, ][!is.na(number.element[2, ])],
					number_of_elements_in_intersect = number.element.intersect,
					mean_of_each_in_factor1 = mean.treatment[1, ][!is.na(mean.treatment[1, ])],
					mean_of_each_in_factor2 = mean.treatment[2, ][!is.na(mean.treatment[2, ])],
					mean_of_each_intersect = mean.intersect,
					var_of_each_in_factor1 = var.treatment[1, ][!is.na(var.treatment[1, ])],
					var_of_each_in_factor2 = var.treatment[2, ][!is.na(var.treatment[2, ])],
					var_of_each_intersect = var.intersect, number_of_all_elements = number, total_mean = mean)
		sol = list(Factor1 = t1, Factor2 = t2, Intersect = t3, Error = t4, Total = t5,
					F_statistic_for_Factor1 = t6, F_statistic_for_Factor2 = t7, F_statistic_for_Intersect = t8, Additional = t9,
					p_value_for_Normality = pnor, p_value_for_Equal_Variance = pvar)
	}else{
		x1 = c(number.factor[1] - 1, SSA, MSA, FA, number.factor[2] - 1, SSB, MSB, FB, (number.factor[1] - 1) * (number.factor[2] - 1), SSAB, MSAB, FAB,
					number - number.factor[1] * number.factor[2], SSE, MSE, NA, number - 1, SSA + SSB + SSAB + SSE, NA, NA)
		n1 = list(c("Factor1", "Factor2", "Intersect", "Error", "Total"), c("Degree of Freedom", "Sum of Squares", "Mean Squares", "F-statistic"))
		mat1 = matrix(data = x1, nrow = 5, ncol = 4, dimnames = n1, byrow = TRUE)
		t1 = as.table(mat1)
		x2 = c(number.element[1, ][!is.na(number.element[1, ])], number, mean.treatment[1, ][!is.na(mean.treatment[1, ])], mean, var.treatment[1, ][!is.na(var.treatment[1, ])], NA)
		n2 = list(c("number", "mean", "variance"), c(name, "Total"))
		mat2 = matrix(data = x2, nrow = 3, ncol = number.factor[1] + 1, dimnames = n2, byrow = TRUE)
		t2 = as.table(mat2)
		mat3 = matrix(data = pA, nrow = 1, dimnames = list("pvalue", ""))
		t3 = as.table(mat3)
		x4 = c(number.element[2, ][!is.na(number.element[2, ])], number, mean.treatment[2, ][!is.na(mean.treatment[2, ])], mean, var.treatment[2, ][!is.na(var.treatment[2, ])], NA)
		n4 = list(c(level, "Total"), c("number", "mean", "variance"))
		mat4 = matrix(data = x4, nrow = number.factor[2] + 1, ncol = 3, dimnames = n4, byrow = FALSE)
		t4 = as.table(mat4)
		mat5 = matrix(data = pB, nrow = 1, dimnames = list("pvalue", ""))
		t5 = as.table(mat5)
		t6 = as.table(number.element.intersect)
		t7 = as.table(mean.intersect)
		t8 = as.table(var.intersect)
		mat9 = matrix(data = pAB, nrow = 1, dimnames = list("pvalue", ""))
		t9 = as.table(mat9)
		mat10 = matrix(data = pnor, nrow = 1, dimnames = list("pvalue", ""))
		t10 = as.table(mat10)
		mat11 = matrix(data = pvar, nrow = 1, dimnames = list("pvalue", ""))
		t11 = as.table(mat11)
		sol = list(ANOVA_table = t1, Additional_information_for_Factor1 = t2, p_value_for_Factor1 = t3, Additional_information_for_Factor2 = t4, p_value_for_Factor2 = t5,
					Number_of_Elements_in_Intersect = t6, Mean_of_Each_Intersect = t7, Var_of_Each_Intersect = t8, p_value_for_Intersect = t9,
					p_value_for_Normality = t10, p_value_for_Equal_Variance = t11)
	}
	return (sol)
}
