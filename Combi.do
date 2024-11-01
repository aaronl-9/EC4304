use data.dta, clear
preserve

drop if _n > 480

global train_n = 280
global train_start = 1
global train_end = 380
// no validation set here (all validation set obs go into train set)
global test_start = $train_end + 1 // 381 
global test_end = $test_start + 99 // 480 
global test_num_forecasts = $test_end - $test_start // (num of forecasts to make is actually 100 but leaving it as 99 for the purpose of the loop later)

// HAC for Granger causality tests

global hac_lag = round(0.75 * (280-12)^(1/3))

local time = tm(1981m1)
newey indpro_1step L(1/12).indpro_1step L(1/12).news_sentiment L(1/12).oecd_cli ///
	if time >= `time' & _n <= $train_end , lag($hac_lag)
	
testparm L(1/12).news_sentiment // Cannot reject null of no Gr. causality
testparm L(1/12).oecd_cli

// AR(1)

/* DO NOT USE 
local step 1 3 6 12
foreach j of local step {
	gen y`j'_forecast = .
	gen y`j'_forecast_err = .
	
	forvalues i = 0/$test_num_forecasts {
		local start = $train_start + `i'
		local end = $train_end + `i'
		
		qui reg F`j'.indpro_1step indpro_1step ///
			in `start'/`end', r
		qui predict temp_fct, xb
		qui replace y`j'_forecast = temp_fct[_n+`j'] if _n == `i'
		drop temp_fct
		
		qui replace y`j'_forecast_err = (y`j'_forecast - F`j'.indpro_1step) if _n == `i'
	}
	egen mse = mean(y`j'_forecast_err^2) 
	scalar rmse = sqrt(mse)
	display "`j'-step-ahead RMSE: " rmse
	drop mse
}
*/

* 3 step ahead/AR(3)
local step = 3
gen y3_forecast = .
forvalues i = 0/$test_num_forecasts { //99
	local start = $train_start + `i'
	local end = $train_end + `i'
	
	qui reg F`step'.indpro_1step L(1/3).indpro_1step in `start'/`end', r
	qui predict temp_fct, xb
	
	local fct_index = `end' + `step'
	if `fct_index' <= _N {
		qui replace y`step'_forecast = temp_fct[`fct_index'] if _n == `fct_index'
	}
	drop temp_fct
}

gen sqe_ar1_3 = ((y3_forecast-indpro_1step)^2)
summarize sqe_ar1_3
local mse_mean = r(mean)
scalar rmsfe_ar1_3 = sqrt(`mse_mean')
di rmsfe_ar1_3
drop sqe_ar1_3

* 6 step ahead/AR(3)
local step = 6
gen y`step'_forecast = .
forvalues i = 0/$test_num_forecasts { //99
	local start = $train_start + `i'
	local end = $train_end + `i'
	
	qui reg F`step'.indpro_1step L(1/3)indpro_1step in `start'/`end', r
	qui predict temp_fct, xb
	
	local fct_index = `end' + `step'
	if `fct_index' <= _N {
		qui replace y`step'_forecast = temp_fct[`fct_index'] if _n == `fct_index'
	}
	drop temp_fct
}

gen sqe_ar1_6 = ((y6_forecast-indpro_1step)^2)
summarize sqe_ar1_6
local mse_mean = r(mean)
scalar rmsfe_ar1_6 = sqrt(`mse_mean')
di rmsfe_ar1_6
drop sqe_ar1_6

* 12 step ahead/AR(3)
local step = 12
gen y`step'_forecast = .
forvalues i = 0/$test_num_forecasts { //99
	local start = $train_start + `i'
	local end = $train_end + `i'
	
	qui reg F`step'.indpro_1step L(1/3).indpro_1step in `start'/`end', r
	qui predict temp_fct, xb
	
	local fct_index = `end' + `step'
	if `fct_index' <= _N {
		qui replace y`step'_forecast = temp_fct[`fct_index'] if _n == `fct_index'
	}
	drop temp_fct
}

gen sqe_ar1_12 = ((y12_forecast-indpro_1step)^2)
summarize sqe_ar1_12
local mse_mean = r(mean)
scalar rmsfe_ar1_12 = sqrt(`mse_mean')
di rmsfe_ar1_12
drop sqe_ar1_12

drop *forecast*

* 3 step ahead/AR(3)
local step = 3
gen y3_forecast = .
forvalues i = 0/$test_num_forecasts { //99
	local start = $train_start + `i'
	local end = $train_end + `i'
	
	qui reg F`step'.indpro_1step L(1/2).indpro_1step in `start'/`end', r
	qui predict temp_fct, xb
	
	local fct_index = `end' + `step'
	if `fct_index' <= _N {
		qui replace y`step'_forecast = temp_fct[`fct_index'] if _n == `fct_index'
	}
	drop temp_fct
}

gen sqe_ar1_3 = ((y3_forecast-indpro_1step)^2)
summarize sqe_ar1_3
local mse_mean = r(mean)
scalar rmsfe_ar1_3 = sqrt(`mse_mean')
di rmsfe_ar1_3
drop sqe_ar1_3

* 6 step ahead/AR(2)
local step = 6
gen y`step'_forecast = .
forvalues i = 0/$test_num_forecasts { //99
	local start = $train_start + `i'
	local end = $train_end + `i'
	
	qui reg F`step'.indpro_1step L(1/2).indpro_1step in `start'/`end', r
	qui predict temp_fct, xb
	
	local fct_index = `end' + `step'
	if `fct_index' <= _N {
		qui replace y`step'_forecast = temp_fct[`fct_index'] if _n == `fct_index'
	}
	drop temp_fct
}

gen sqe_ar1_6 = ((y6_forecast-indpro_1step)^2)
summarize sqe_ar1_6
local mse_mean = r(mean)
scalar rmsfe_ar1_6 = sqrt(`mse_mean')
di rmsfe_ar1_6
drop sqe_ar1_6

* 12 step ahead/AR(2)
local step = 12
gen y`step'_forecast = .
forvalues i = 0/$test_num_forecasts { //99
	local start = $train_start + `i'
	local end = $train_end + `i'
	
	qui reg F`step'.indpro_1step L(1/2).indpro_1step in `start'/`end', r
	qui predict temp_fct, xb
	
	local fct_index = `end' + `step'
	if `fct_index' <= _N {
		qui replace y`step'_forecast = temp_fct[`fct_index'] if _n == `fct_index'
	}
	drop temp_fct
}

gen sqe_ar1_12 = ((y12_forecast-indpro_1step)^2)
summarize sqe_ar1_12
local mse_mean = r(mean)
scalar rmsfe_ar1_12 = sqrt(`mse_mean')
di rmsfe_ar1_12
drop sqe_ar1_12

restore
