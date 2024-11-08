use "Master.dta", clear // can replace with "data.dta"

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
	
***** 12 Step ahead forecast *****

/* Notations:

ar_1sa_rmse: AR 1 step ahead RMSE
adl_3sa_rmse: ADL with CLI 3 step ahead RMSE
adln_12sa_rmse: ADL with CLI + news 12 step ahead RMSE

*/

* AR(2) forecasts
local horizon = 12
gen ar_12sa_fit = .
forval i = 0/$test_num_forecasts {
	local current_train_start = $train_start + `i'
	local current_train_end = $train_end + `i' - `horizon' + 1
	// hac standard errors 
	qui newey indpro_12step L(12/13).indpro_1step ///
		if _n >= `current_train_start' & _n <= `current_train_end', lag($hac_lag)
	predict temp_fit, xb 
	replace ar_12sa_fit = temp_fit if _n == ($test_start + `i')
	drop temp_fit 
}

gen ar_12sa_sqerror = (indpro_12step - ar_12sa_fit)^2
egen ar_12sa_mse = mean(ar_12sa_sqerror)
gen ar_12sa_rmse = sqrt(ar_12sa_mse)
di ar_12sa_rmse // Verify RMSE

* ADL with CLI VS benchmark AR(AIC) [ADL(11,5)]

local horizon = 12
gen adl_12sa_fit = .
forval i = 0/$test_num_forecasts {
	local current_train_start = $train_start + `i'
	local current_train_end = $train_end + `i' - `horizon' + 1
	// hac standard errors 
	qui newey indpro_12step L(12/22).indpro_1step L(12/16).cli ///
		if _n >= `current_train_start' & _n <= `current_train_end', lag($hac_lag)
	predict temp_fit, xb 
	replace adl_12sa_fit = temp_fit if _n == ($test_start + `i')
	drop temp_fit 
}

gen adl_12sa_sqerror = (indpro_12step - adl_12sa_fit)^2
egen adl_12sa_mse = mean(adl_12sa_sqerror)
gen adl_12sa_rmse = sqrt(adl_12sa_mse)
di adl_12sa_rmse // Verify

* ADL with CLI+news [ADL(10,4,12)]

local horizon = 12
gen adln_12sa_fit = .
forval i = 0/$test_num_forecasts {
	local current_train_start = $train_start + `i' 
	local current_train_end = $train_end + `i' - `horizon' + 1
	// hac standard errors 
	qui newey indpro_12step L(12/21).indpro_1step L(12/15).cli L(12/23).news ///
		if _n >= `current_train_start' & _n <= `current_train_end', lag($hac_lag)
	predict temp_fit, xb 
	replace adln_12sa_fit = temp_fit if _n == ($test_start + `i')
	drop temp_fit 
}

gen adln_12sa_sqerror = (indpro_12step - adln_12sa_fit)^2
egen adln_12sa_mse = mean(adln_12sa_sqerror)
gen adln_12sa_rmse = sqrt(adln_12sa_mse)
di adln_12sa_rmse // Verify RMSE

/* Loss differentials */
* AR(2) & ADL(11,5) & ADL(10,4,12)
gen ld1_12sa = ar_12sa_sqerror - adl_12sa_sqerror
gen ld2_12sa = adl_12sa_sqerror - adln_12sa_sqerror

dfuller ld1_12sa // stationary at p=0.0003
dfuller ld2_12sa // stationary at p=0.0006

ac ld1_12sa
ac ld2_12sa
// corrgram ld1_12sa 

// DM
dmariano indpro_12step ar_12sa_fit adl_12sa_fit if _n >= $test_start , crit(MSE) maxlag(4) kernel(bartlett)
dmariano indpro_12step adl_12sa_fit adln_12sa_fit if _n >= $test_start , crit(MSE) maxlag(4) kernel(bartlett)