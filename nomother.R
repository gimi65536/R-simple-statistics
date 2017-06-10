Rank = function(x1, x2){
	data = c(x1, x2)
	n = c(length(x1), length(x2))
	s = vector()
	r = rank(data)
	level1 = levels(as.factor(data))
	level2 = levels(as.factor(r))
	sol = list()
	x = list(x1, x2)
	for(i in 1:2){
		sol[[i]] = vector()
		for(j in 1 : n[i]){
			sol[[i]] = append(sol[[i]], as.numeric(level2[which(level1 == x[[i]][j])]))
		}
		s = append(s, sum(sol[[i]]))
		length(sol[[i]]) = max(n)
	}
	sol[[3]] = s
	return (sol)
}

Rank2 = function(df){
	data = vector()
	n = vector()
	if(length(df) != 0){
		for(i in 1 : length(df)){
			x = df[[i]]
			x = x[!is.na(x)]
			data = append(data, x)
			n = append(n, length(x))
		}
	}else{
		return (list())
	}
	s = vector()
	k = length(df)
	r = rank(data)
	accu = 0
	sol = list()
	for(i in 1 : k){
		sol[[i]] = vector()
		sol[[i]] = append(sol[[i]], r[(accu + 1) : (accu + n[i])])
		s = append(s, sum(sol[[i]]))
		length(sol[[i]]) = max(n)
		accu = accu + n[i]
	}
	sol[[k + 1]] = s
	return (sol)
}

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

WilcoxonTest = function(..., df = NULL, index = c(1, 2), alternative = "two.sided"){
	if(alternative != "greater" && alternative != "less" && alternative != "two.sided"){
		alternative = "two.sided"
	}
	if(is.null(df)){
		df = list(...)
		df = as.data.frame(df)
	}
	x1 = df[[index[1]]]
	x2 = df[[index[2]]]
	x1 = x1[!is.na(x1)]
	x2 = x2[!is.na(x2)]
	x = list(x1, x2)
	n = c(length(x1), length(x2))
	if(is.null(names(df))){
		name = c("Treatment 1", "Treatment 2")
	}else{
		name = c(names(df)[index[1]], names(df)[index[2]])
	}
	l = Rank(x1, x2)
	s = l[[3]]
	x[[3]] = l[[1]]
	x[[4]] = l[[2]]
	m = matrix(c(x[[3]], s[1], x[[4]], s[2]), ncol = 2, dimnames = list(c(rep_len("", max(n)), "Total"), name))
	E1 = n[1] * (sum(n) + 1) / 2
	sigma = sqrt(n[1] * n[2] * (sum(n) + 1) / 12)
	z = (m[max(n) + 1, 1] - E1) / sigma
	p = pnorm(z)
	if(alternative == "greater"){
		p = 1 - p
	}else if(alternative == "two.sided"){
		if(p > 0.5){
			p = 1 - p
		}
		p = 2 * p
	}
	m2 = matrix(c(n, n * (sum(n) + 1) / 2, sigma, z, p), ncol = 1, dimnames = list(c(paste("Number", 1:2), paste("Expected", 1:2), "Sigma", "zvalue", "pvalue"), ""))
	return (list(Rank = as.table(m), Solution = as.table(m2)))
}

signed_test = function(x, more = "equal"){
	if(more != "negative" && more != "positive" && more != "equal"){
		more = "equal"
	}
	n = sum(x != 0)
	neg = sum(x < 0)
	pos = sum(x > 0)
	if(more == "negative"){
		z = (neg - n / 2) / sqrt(n / 4)
	}else{
		z = (pos - n / 2) / sqrt(n / 4)
	}
	p = pnorm(z)
	if(more == "equal"){
		if(p > 0.5){
			p = 1 - p
		}
		p = 2 * p
	}else{
		p = 1 - p
	}
	m = matrix(c(n, pos, neg, z, p), ncol = 1, dimnames = list(c("All nonzero", "Positive", "Negative", paste("zvalue for", more, "mode"), "pvalue"), ""))
	return(as.table(m))
}

WilcoxonSignTest = function(..., df = NULL, index = c(1, 2), more = "equal"){
	if(more != "negative" && more != "positive" && more != "equal"){
		more = "equal"
	}
	if(is.null(df)){
		df = list(...)
		df = as.data.frame(df)
	}
	x1 = df[[index[1]]]
	x2 = df[[index[2]]]
	remove = is.na(x1) | is.na(x2)
	x1 = x1[!remove]
	x2 = x2[!remove]
	if(is.null(names(df))){
		name = c("Treatment 1", "Treatment 2")
	}else{
		name = c(names(df)[index[1]], names(df)[index[2]])
	}
	name = append(name, c("Difference", "Positive Rank", "Negative Rank"))
	array = x1 - x2
	n = length(array)
	n_ = length(array[array != 0])
	r = rank(abs(array[array != 0]))
	accu = 0
	pos = rep_len(NA, n)
	neg = rep_len(NA, n)
	for(i in 1 : n){
		if(array[i] > 0){
			pos[i] = r[i - accu]
		}else if(array[i] < 0){
			neg[i] = r[i - accu]
		}else{
			accu = accu + 1
		}
	}
	rs_p = sum(pos[!is.na(pos)])
	rs_n = sum(neg[!is.na(neg)])
	m1 = matrix(c(x1, NA, x2, NA, array, NA, pos, rs_p, neg, rs_n), ncol = 5, dimnames = list(c(rep_len("", n), "Total"), name))
	m1 = as.table(m1)
	exp = n_ * (n_ + 1) / 4
	sigma = sqrt(n_ * (n_ + 1) * (2 * n_ + 1) / 24)
	if(more == "negative"){
		z = (rs_n - exp) / sigma
	}else{
		z = (rs_p - exp) / sigma
	}
	p = pnorm(z)
	if(more == "equal"){
		if(p > 0.5){
			p = 1 - p
		}
		p = 2 * p
	}else{
		p = 1 - p
	}
	m2 = matrix(c(n_, exp, sigma, z, p), ncol = 1, dimnames = list(c("Number", "Expected", "Sigma", paste("zvalue for", more, "mode"), "pvalue"), ""))
	return (list(Rank = as.table(m1), Solution = as.table(m2)))
}

KWTest = function(..., df = NULL, index = NULL){
	if(is.null(df)){
		args = list(...)
		if(length(index) != 0){
			args = args[index]
		}
	}else{
		args = as.list(df)
		if(length(index) != 0){
			args = args[index]
		}
	}
	name = getName(args)
	l = Rank2(args)
	k = length(args)
	n = vector()
	data = vector()
	H = 0
	for(i in 1 : k){
		x = l[[i]]
		x = x[!is.na(x)]
		data = append(data, append(l[[i]], l[[k + 1]][i]))
		n = append(n, length(x))
		H = H + sum(x) ^ 2 / n[i]
	}
	m1 = matrix(data, ncol = k, dimnames = list(c(rep_len("", max(n)), "Total"), name))
	N = sum(n)
	H = (12 / N / (N + 1) * H) - 3 * (N + 1)
	p = 1 - pchisq(H, k - 1)
	m2 = matrix(c(k, k - 1, H, p), ncol = 1, dimnames = list(c("Treatment", "Degree of Freedom", "H", "pvalue"), ""))
	return (list(Rank = as.table(m1), Solution = as.table(m2)))
}

FriedmanTest = function(..., df = NULL, index = NULL){
	if(is.null(df)){
		args = list(...)
		if(length(index) != 0){
			args = args[index]
		}
	}else{
		args = as.list(df)
		if(length(index) != 0){
			args = args[index]
		}
	}
	name = getName(args)
	b = vector()
	k = length(args)
	for(i in 1 : k){
		x = args[[i]]
		x = x[!is.na(x)]
		b = append(b, length(x))
		args[[i]] = x
	}
	b = min(b)
	m0 = matrix(ncol = k, nrow = b, dimnames = list(paste("Block", 1 : b), name))
	m = matrix(ncol = k, nrow = b + 1, dimnames = list(c(paste("Block", 1 : b), "Total"), name))
	for(i in 1 : k){
		x = args[[i]]
		x = x[!is.na(x)]
		m0[, i] = x[1 : b]
	}
	for(i in 1 : b){
		m[i, ] = rank(m0[i, ])
	}
	m[b + 1, ] = apply(m, 2, sum, na.rm = TRUE)
	T = m[b + 1, ]
	Fr = (12 * sum(T ^ 2) / b / k / (k + 1)) - 3 * b * (k + 1)
	p = 1 - pchisq(Fr, k - 1)
	m2 = matrix(c(k, b, k - 1, Fr, p), ncol = 1, dimnames = list(c("Treatment", "Block", "Degree of Freedom", "Fr", "pvalue"), ""))
	return (list(Data = as.table(m0), Rank = as.table(m), Solution = as.table(m2)))
}

PearsonRankCorrelation = function(..., df = NULL, index = c(1, 2), alternative = "equal"){
	if(alternative != "greater" && alternative != "less" && alternative != "two.sided"){
		alternative = "two.sided"
	}
	if(is.null(df)){
		df = list(...)
		df = as.data.frame(df)
		names(df) = c("Treatment 1", "Treatment 2")
	}
	x1 = df[[index[1]]]
	x2 = df[[index[2]]]
	x1 = x1[!is.na(x1)]
	x2 = x2[!is.na(x2)]
	x = list(x1, x2)
	n = c(length(x1), length(x2))
	n = min(n)
	if(is.null(names(df))){
		name = c("Treatment 1", "Treatment 2")
	}else{
		name = c(names(df)[index[1]], names(df)[index[2]])
	}
	m0 = matrix(nrow = n, ncol = 2, dimnames = list(rep_len("", n), name))
	m0[, 1] = x1[1 : n]
	m0[, 2] = x2[1 : n]
	m = matrix(nrow = n, ncol = 2, dimnames = list(rep_len("", n), name))
	m[, 1] = rank(x1)
	m[, 2] = rank(x2)
	print(name)
	print(sum(m[, 1]))
	print(sum(m[, 2]))
	print(sum(m[, 1] * m[, 2]))
	print(sum(m[, 1] ^ 2))
	print(sum(m[, 2] ^ 2))
	Sa = sqrt((sum(m[, 1] ^ 2) - sum(m[, 1]) ^ 2 / n) / (n - 1))
	Sb = sqrt((sum(m[, 2] ^ 2) - sum(m[, 2]) ^ 2 / n) / (n - 1))
	Sab = (sum(m[, 1] * m[, 2]) - sum(m[, 1]) * sum(m[, 2]) / n) / (n - 1)
	cor = Sab / Sa / Sb
	z = NULL
	if(n <= 30){
		p = cor.test(x1, x2, method = "spearman", alternative = alternative)[[3]]
	}else{
		z = cor / sqrt(1 / (n - 1))
		p = pnorm(z)
		if(alternative == "greater"){
			p = 1 - p
		}else if(alternative == "two.sided"){
			if(p > 0.5){
				p = 1 - p
			}
			p = 2 * p
		}
	}
	m2 = matrix(c(n, Sa, Sb, Sab, cor, z, p), ncol = 1, dimnames = list(c("Number", "Sa", "Sb", "Sab", "Cor", "zvalue", paste("pvalue for", alternative)), ""))
	return (list(Data = as.table(m0), Rank = as.table(m), Solution = as.table(m2)))
}
