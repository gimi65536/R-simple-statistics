ftest = function(x1, x2, ratio = 1, alternative = c("two.sided", "greater", "less")){
	alternative = match.arg(alternative)
	x1 = x1[!is.na(x1)]
	x2 = x2[!is.na(x2)]
	n = c(length(x1), length(x2))
	v = c(var(x1), var(x2))
	f = v[1] / v[2] / ratio
	p = pf(f, n[1] - 1, n[2] - 1)
	if(alternative == "greater"){
		p = 1 - p
	}else if(alternative == "two.sided"){
		if(p > 0.5){
			p = 1 - p
		}
		p = 2 * p
	}
	sol = matrix(data = c(n[1], v[1], f, p, n[2], v[2], NA, NA), byrow = TRUE, nrow = 2, dimnames = list(1:2, c("number", "variance", "f", "pvalue")))
	return(as.table(sol))
}

mean_differ = function(x1, x2, difference = 0, alternative = c("two.sided", "greater", "less"), alpha = 0.05){
	alternative = match.arg(alternative)
	x1 = x1[!is.na(x1)]
	x2 = x2[!is.na(x2)]
	n = c(length(x1), length(x2))
	var_equal = ftest(x1, x2)[1, 4] > alpha
	if(var_equal){
		t = mean(x1) - mean(x2)
		sp2 = ((n[1] - 1) * var(x1) + (n[2] - 1) * var(x2)) / (n[1] + n[2] - 2)
		t = t / sqrt(sp2 / n[1] + sp2 / n[2])
		d = n[1] + n[2] - 2
	}else{
		sp2 = NA
		s = c(var(x1) / n[1], var(x2) / n[2])
		d = sum(s) ^ 2 / (s[1] ^ 2 / (n[1] - 1) + s[2] ^ 2 / (n[2] - 1))
		t = (mean(x1) - mean(x2)) / sqrt(sum(s))
	}
	p = pt(t, d)
	if(alternative == "greater"){
		p = 1 - p
	}else if(alternative == "two.sided"){
		if(p > 0.5){
			p = 1 - p
		}
		p = 2 * p
	}
	sol = matrix(data = c(n[1], sp2, d, t, p, n[2], NA, NA, NA, NA), byrow = TRUE, nrow = 2, dimnames = list(1:2, c("number", "sp2", "freedom", "t", "pvalue")))
	return(as.table(sol))
}