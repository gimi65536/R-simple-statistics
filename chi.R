godness_fit = function(observe, exp, ratio = FALSE, test_for_normality = FALSE){
	n = max(length(observe), length(exp))
	sol = matrix(nrow = n + 1, ncol = 4, dimnames = list(c(paste("Category", 1:n), "Total"), c("Observed Frequency", "Expected Frequency", "f - e", "Chi-sq")))
	ob = rep_len(observe, n)
	ex = rep_len(exp, n)
	if(ratio){
		ex = sum(ob) * ex
	}
	sol[, 1] = c(ob, sum(ob))
	sol[, 2] = c(ex, sum(ex))
	sol[, 3] = c(ob - ex, NA)
	chi = (ob - ex) ^ 2 / ex
	sol[, 4] = c(chi, sum(chi))
	if(test_for_normality){
		p = 1 - pchisq(sum(chi), n - 3)
	}else{
		p = 1 - pchisq(sum(chi), n - 1)
	}
	return (list(table = as.table(sol), pvalue = p))
}

test_normality = function(x, n){
	if(n <= 3){
		return (NULL)
	}
	x = sort(x)
	kind_table = table(data.frame(x = x)$x)
	kind = length(kind_table)
	tmp = n
	five = 5
	for(i in 1:kind){
		five = five - kind_table[i]
		if(five <= 0){
			tmp = tmp - 1
			if(tmp == 0){
				break
			}
			five = 5
		}
	}
	if(tmp > 0){
		return (NULL)
	}
	exp = length(x) / n
	accu = 0
	cut = vector()
	for(i in 1:length(x)){
		accu = accu + kind_table[i]
		if(accu == exp){
			cut = append(cut, i)
			accu = 0
		}else if(accu < exp && i < length(x) && accu + kind_table[i + 1] > exp){
			if(exp - accu <= accu + kind_table[i + 1] - exp){
				cut = append(cut, i)
				accu = accu - exp
			}else{
				cut = append(cut, i + 1)
				accu = accu + kind_table[i + 1] - exp
				i = i + 1
			}
		}
		if(length(cut) == n - 1){
			break
		}
	}
	if(length(cut) != n - 1){
		return (NULL)
	}
	z = (x - mean(x)) / sd(x)
	critical = vector(length = n + 1)
	ob = vector(length = n); ex = vector(length = n)
	ob[1] = sum(kind_table[1:cut[1]])
	critical[1] = -Inf; critical[n + 1] = Inf
	for(i in 2:n){
		critical[i] = mean(z[cut[i - 1]], z[cut[i - 1] + 1])
		if(i == n){
			ob[n] = sum(kind_table[(cut[n - 1] + 1):kind])
		}else{
			ob[i] = sum(kind_table[(cut[i - 1] + 1):cut[i]])
		}
	}
	if(sum(ob < 5) > 0){
		return (NULL)
	}
	criticalp = pnorm(critical)
	for(i in 1:n){
		ex[i] = criticalp[i + 1] - criticalp[i]
	}
	sol = godness_fit(ob, ex, ratio = TRUE, test_for_normality = TRUE)
	for(i in 1:n){
		dimnames(sol$table)[[1]][i] = paste(critical[i], "<= z <=", critical[i + 1])
	}
	return (sol)
}

contingency_table_nominal = function(..., df = NULL){
	if(is.null(df)){
		df = list(...)
		df = as.data.frame(df)
	}
	f = as.factor(t(df)) #by row
	#f= as.factor(as.matrix(df)) by column
	col = length(df)
	fas = levels(f)
	row = length(fas)
	colname = names(df)
	conj = matrix(nrow = row, ncol = col, dimnames = list(fas,colname))
	for(i in 1 : col){
		x = df[[i]]
		for(j in fas){
			conj[j, i] = sum(x == j)
		}
	}
	return (contingency_table(conj))
}

contingency_table_2factor = function(..., df = NULL, index = c(1, 2)){
	if(is.null(df)){
		df = list(...)
		df = as.data.frame(df)
	}
	x1 = df[[index[1]]]
	x2 = df[[index[2]]]
	x1 = x1[!is.na(x1)]
	x2 = x2[!is.na(x2)]
	if(length(x1) != length(x2)){
		return (NULL)
	}
	if(is.null(names(df))){
		name = c("Factor 1", "Factor 2")
	}else{
		name = c(names(df)[index[1]], names(df)[index[2]])
	}
	f1 = as.factor(x1)
	f2 = as.factor(x2)
	level1 = levels(f1)
	level2 = levels(f2)
	nama = list(paste(name[2], level2, sep = '.'), paste(name[1], level1, sep = '.'))
	m = matrix(nrow = length(level2), ncol = length(level1), dimnames = nama)
	for(i in 1 : length(level2)){
		for(j in 1 : length(level1)){
			m[i, j] = sum(x1 == level1[j] & x2 == level2[i])
		}
	}
	return (contingency_table(m))
}

contingency_table = function(m){
	if(typeof(m) != "matrix"){
		m = as.matrix(m)
	}
	row = nrow(m)
	col = ncol(m)
	conj = matrix(nrow = row + 1, ncol = col + 1)
	if(is.null(rownames(m))){
		rownames(conj) = c(rep_len(NA, row), "Total")
	}else{
		rownames(conj) = c(rownames(m), "Total")
	}
	if(is.null(colnames(m))){
		colnames(conj) = c(rep_len(NA, col), "Total")
	}else{
		colnames(conj) = c(colnames(m), "Total")
	}
	nama = list(rownames(conj), colnames(conj))
	for(i in 1 : row){
		for(j in 1 : col){
			conj[i, j] = m[i, j]
		}
		conj[i, col + 1] = sum(conj[i, -(col + 1)])
	}
	for(i in 1 : col){
		conj[row + 1, i] = sum(conj[-(row + 1), i])
	}
	conj[row + 1, col + 1] = sum(conj[-(row + 1), col + 1])
	fre = conj / conj[row + 1, col + 1]
	exp = matrix(nrow = row, ncol = col, dimnames = dimnames(m))
	expfre = matrix(nrow = row, ncol = col, dimnames = dimnames(m))
	chi = matrix(nrow = row + 1, ncol = col + 1, dimnames = nama)
	for(i in 1 : row){
		for(j in 1 : col){
			exp[i, j] = conj[row + 1, j] * conj[i, col + 1] / conj[row + 1, col + 1]
			expfre[i, j] = fre[row + 1, j] * fre[i, col + 1]
			chi[i, j] = (conj[i, j] - exp[i, j]) ^ 2 / exp[i, j]
		}
		chi[i, col + 1] = sum(chi[i, ][-(col + 1)])
	}
	for(i in 1 : (col + 1)){
		chi[row + 1, i] = sum(chi[, i][-(row + 1)])
	}
	p = 1 - pchisq(chi[row + 1, col + 1], (row - 1) * (col - 1))
	return (list(conjunction_matrix = conj, frequency = fre, expected_time = exp, expected_frequency = expfre, chisq = chi, pvalue = p))
}