# R-simple-statistics
應修習台大資管統計學一下的同學們殷殷期盼，將本人所撰寫的弱弱的R code貼上來。

## ANOVA.R  
**這個檔案中所有的函數都需要在輸出的地方進行優化**    
**這個檔案充斥許多名詞錯誤與麻煩的寫法，有空的話必須大幅改寫**  

### ANOVA.oneway
One way ANOVA  
(...)為資料，若沒給df東西則使用這裡的資料  
buildin = false 為是否使用R內建的ANOVA  
table = false 為是否使用table形式輸出（否則就是一串list）  
df = NULL 為使用的data frame，如果有東西就**只**使用這邊的資料，否則使用(...)的資料  
index = NULL為資料中要檢定的部分，若是NULL則為全部資料  

### cross_differ
資料輸入方法與ANOVA.oneway相同，以one way ANOVA的結果來檢測兩筆資料的差異是否顯著。  
(...)為資料，若沒給df東西則使用這裡的資料  
alpha = 0.05 顯著水準
way = "LSD" 使用哪種方法檢測差異。接受的輸入有"LSD" "Bonferroni" "Tukey"。*函數的寫法需要修正，使用match.arg才是建議方法*  
df = NULL 為使用的data frame，如果有東西就**只**使用這邊的資料，否則使用(...)的資料  
label = false data frame的第一列是否要列入。*這個功能挺廢的，或許當初在寫的時候有什麼特殊需求吧？*  
index = NULL為資料中要檢定的部分，若是NULL則為全部資料  

### ANOVA.twoway
用於Experimental design的ANOVA
(...)為資料，若沒給df東西則使用這裡的資料  
buildin = false 為是否使用R內建的ANOVA  
table = false 為是否使用table形式輸出（否則就是一串list）  
df = NULL 為使用的data frame，如果有東西就**只**使用這邊的資料，否則使用(...)的資料  
index = NULL為資料中要檢定的部分，若是NULL則為全部資料  

### ANOVA.2factor
雙因子ANOVA
df 為使用的data frame  
buildin = false 為是否使用R內建的ANOVA  
table = false 為是否使用table形式輸出（否則就是一串list）  
balance = true 每個treatment內的元素數量是否相等。*目前是閒置變數，元素數量永遠視為相等*  
index = NULL為資料中要檢定的部分，若是NULL則為全部資料  

## chi.R

### godness_fit
檢測觀測結果與期望值是否相等  
observe 觀測結果  
exp 期望值  
ratio = false 期望值是否為比例  
test_for_normality = false 是否為常態分布檢驗  

### test_normality
給定一筆資料，檢測資料是否有常態分布  
x 資料  
n 預計要分成幾組數值  

### contingency_table_nominal
給定一組data frame，每個column當成第一因子的各個level，而每列資料的數值本身代表第二因子的level  
例如column有"驅逐艦" "戰艦"，而表格元素 1:"無破" 2:"大破"，則在"驅逐艦"column中出現的"2"代表一筆「"驅逐艦"且"大破"」的資料，寫入contingency table  
(...)為資料，若沒給df東西則使用這裡的資料  
df = NULL 為使用的data frame，如果有東西就**只**使用這邊的資料，否則使用(...)的資料  

### contingency_table_2factor
給定一組data frame，第一個column中出現的值為第一因子的各個level，第二個column中出現的值為第二因子的各個level  
例如第一個column代表"艦種"，1:"驅逐艦" 2:"戰艦"，第二個column代表"損傷"，1:"無破" 2:"大破"，則在某一row中，有第一個column=1 第二個column=2，代表一筆「"驅逐艦"且"大破"」的資料，寫入contingency table  
(...)為資料，若沒給df東西則使用這裡的資料  
df = NULL 為使用的data frame，如果有東西就**只**使用這邊的資料，否則使用(...)的資料  
index = c(1, 2) 把哪兩column視為因子  

### contingency_table
給定一個二維矩陣，算出column與row的獨立性
m 給定矩陣

## difference.R
**proportion相關的檢定並沒有寫上去**  

### ftest
檢測兩筆資料的變異數之比例是否大/小/等於ratio  
x1 第一筆資料  
x2 第二筆資料  
ratio = 1 虛無假設中 x1的變異數/x2的變異數 等於 ratio  
alternative = "two.sided" H1的描述是變異數比例大於(greater)/小於(less)/不等於(two.sided) ratio *函數的寫法需要修正，使用match.arg才是建議方法*  

### mean_differ
檢測兩筆資料的平均數是否為difference
x1 第一筆資料  
x2 第二筆資料  
difference = 0 虛無假設中 x1的平均數/x2的平均數 等於 difference  
alternative = "two.sided" H1的描述是變異數比例大於(greater)/小於(less)/不等於(two.sided) ratio *函數的寫法需要修正，使用match.arg才是建議方法*  
alpha = 0.05 在檢定兩資料變異數是否相等時要用多少的顯著水準  

## forecast.R
**繪圖時需要ggplot2**  
**使用自定義的"FORECAST"類型物件**  
**尚未寫出季節性趨勢**  

### average_value
算出資料使用average value方法算出來的預測值  
本期的值為前面所有期數的值之平均  
x 資料  

### moving_average
算出資料使用移動平均法算出來的預測值  
本期的值為前n期的值之平均  
x 資料  
period = 3 往前平均的期數  

### exponential_smoothing
算出資料使用exponential smoothing方法算出來的預測值  
x 資料  
alpha = 0.2  

### holt_exponential_smoothing
算出資料使用Holt's exponential smoothing方法算出來的預測值  
x 資料  
alpha = 0.9  
beta = 0.8  
ET0 = 1  
F0 = mean(x) 假定第零期的預測值  
X0 = F0 假定第零期的資料  

### holt_forecast
算出使用Holt's exponential smoothing方法算出來的n期後預測值  
x holt_exponential_smoothing函數輸出的物件  
period 期數  

### draw_time_series
畫出每個預測的趨勢圖  
data list，裝著所有要畫的FORECAST物件  
color = c("black", "red", "blue", "green", "pink") 要使用的顏色，自動循環  
(...) 多餘參數，目前無用  

### make_regression
把迴歸物件轉換成趨勢物件  
lm 迴歸物件，也就是內建函數lm輸出的物件，"lm"型態  
var 會更動的自變數「名稱」  
range 自變數的值  
(...) 輸入其他固定的自變數的值  

### make_package
把資料轉成趨勢物件  
y 預測值向量，最後一個值為未來一期的預測值  
x = NULL 原值向量，代表原始資料  

### mean_error
計算趨勢的mean error  
x FORECAST物件  
period = NULL 拿尾n期進行計算，若為NULL則拿預測值與原始值都有的期間來計算  
ori = NULL 原始資料，若為NULL但x有原始資料則會拿x的原始資料  

### mean_percentage_error
計算趨勢的mean percentage error  
x FORECAST物件  
period = NULL 拿尾n期進行計算，若為NULL則拿預測值與原始值都有的期間來計算  
ori = NULL 原始資料，若為NULL但x有原始資料則會拿x的原始資料  

### mean_absolute_deviation
計算趨勢的mean absolute deviation  
x FORECAST物件  
period = NULL 拿尾n期進行計算，若為NULL則拿預測值與原始值都有的期間來計算  
ori = NULL 原始資料，若為NULL但x有原始資料則會拿x的原始資料  

### mean_square_error
計算趨勢的mean square error  
x FORECAST物件  
period = NULL 拿尾n期進行計算，若為NULL則拿預測值與原始值都有的期間來計算  
ori = NULL 原始資料，若為NULL但x有原始資料則會拿x的原始資料  

### mean_absolute_percentage_error
計算趨勢的mean absolute percentage error  
x FORECAST物件  
period = NULL 拿尾n期進行計算，若為NULL則拿預測值與原始值都有的期間來計算  
ori = NULL 原始資料，若為NULL但x有原始資料則會拿x的原始資料  

### compare_forecast
比較一串趨勢物件個別的誤差  
data list，裝著所有要比較的FORECAST物件  

# 待補充...
