use "Master.dta", clear
preserve

global train_n = 380
global train_start = 1
global train_end = $train_n

// val_start = $train_end + 1 // 281
// global val_end = $val_start + 99 // 380

global test_start = $train_end + 1 // 381
global test_end = $test_start + 99 // 480 

global hac_lag = floor(0.75 * (380-12)^(1/3))

global test_num_forecasts = $test_end - $test_start // 0-99 loop

drop if _n > $test_end
tsset time
	
***** 3 Step ahead forecast *****

/* Notations:

ar_1sa_rmse: AR 1 step ahead RMSE
adl_3sa_rmse: ADL with CLI 3 step ahead RMSE
adln_12sa_rmse: ADL with CLI + news 12 step ahead RMSE

*/

* AR(12) forecasts
local horizon = 6
gen ar_6sa_fit = .
forval i = 0/$test_num_forecasts {
	local current_train_start = $train_start + `i'
	local current_train_end = $train_end + `i' - `horizon' + 1
	// hac standard errors 
	qui newey indpro_6step L(6/17).indpro_1step ///
		if _n >= `current_train_start' & _n <= `current_train_end', lag($hac_lag)
	predict temp_fit, xb 
	replace ar_6sa_fit = temp_fit if _n == ($test_start + `i')
	drop temp_fit 
}

gen ar_6sa_sqerror = (indpro_6step - ar_6sa_fit)^2
egen ar_6sa_mse = mean(ar_6sa_sqerror)
gen ar_6sa_rmse = sqrt(ar_6sa_mse)
di ar_6sa_rmse // Verify

* ADL with CLI VS benchmark AR(AIC) [ADL(9,4)]

local horizon = 6
gen adl_6sa_fit = .
forval i = 0/$test_num_forecasts {
	local current_train_start = $train_start + `i'
	local current_train_end = $train_end + `i' - `horizon' + 1
	// hac standard errors 
	qui newey indpro_6step L(6/14).indpro_1step L(6/9).cli ///
		if _n >= `current_train_start' & _n <= `current_train_end', lag($hac_lag)
	predict temp_fit, xb 
	replace adl_6sa_fit = temp_fit if _n == ($test_start + `i')
	drop temp_fit 
}

gen adl_6sa_sqerror = (indpro_6step - adl_6sa_fit)^2
egen adl_6sa_mse = mean(adl_6sa_sqerror)
gen adl_6sa_rmse = sqrt(adl_6sa_mse)
di adl_6sa_rmse // Verify

* ADL with CLI+news [ADL(9,4,1)]

local horizon = 6
gen adln_6sa_fit = .
forval i = 0/$test_num_forecasts {
	local current_train_start = $train_start + `i' 
	local current_train_end = $train_end + `i' - `horizon' + 1
	// hac standard errors 
	qui newey indpro_6step L(6/14).indpro_1step L(6/9).cli L(6).news ///
		if _n >= `current_train_start' & _n <= `current_train_end', lag($hac_lag)
	predict temp_fit, xb 
	replace adln_6sa_fit = temp_fit if _n == ($test_start + `i')
	drop temp_fit 
}

gen adln_6sa_sqerror = (indpro_6step - adln_6sa_fit)^2
egen adln_6sa_mse = mean(adln_6sa_sqerror)
gen adln_6sa_rmse = sqrt(adln_6sa_mse)
di adln_6sa_rmse // Verify

/* Loss differentials */
* AR(12) & ADL(9,4) & ADL(9,4,1)
gen ld1_6sa = ar_6sa_sqerror - adl_6sa_sqerror
gen ld2_6sa = adl_6sa_sqerror - adln_6sa_sqerror

dfuller ld1_6sa // 
dfuller ld2_6sa // 

ac ld1_6sa
ac ld2_6sa
corrgram ld1_6sa

// DM
dmariano indpro_6step ar_6sa_fit adl_6sa_fit if _n >= $test_start , crit(MSE) maxlag(4) kernel(bartlett)
dmariano indpro_6step adl_6sa_fit adln_6sa_fit if _n >= $test_start , crit(MSE) maxlag(4) kernel(bartlett)