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
local horizon = 3
gen ar_3sa_fit = .
forval i = 0/$test_num_forecasts {
	local current_train_start = $train_start + `i'
	local current_train_end = $train_end + `i' - `horizon' + 1
	// hac standard errors 
	qui newey indpro_3step L(3/14).indpro_1step ///
		if _n >= `current_train_start' & _n <= `current_train_end', lag($hac_lag)
	predict temp_fit, xb 
	replace ar_3sa_fit = temp_fit if _n == ($test_start + `i')
	drop temp_fit 
}

gen ar_3sa_sqerror = (indpro_3step - ar_3sa_fit)^2
egen ar_3sa_mse = mean(ar_3sa_sqerror)
gen ar_3sa_rmse = sqrt(ar_3sa_mse)
di ar_3sa_rmse // Verify

* ADL with CLI VS benchmark AR(AIC) [ADL(10,5)]

local horizon = 3
gen adl_3sa_fit = .
forval i = 0/$test_num_forecasts {
	local current_train_start = $train_start + `i'
	local current_train_end = $train_end + `i' - `horizon' + 1
	// hac standard errors 
	qui newey indpro_3step L(3/12).indpro_1step L(3/7).cli ///
		if _n >= `current_train_start' & _n <= `current_train_end', lag($hac_lag)
	predict temp_fit, xb 
	replace adl_3sa_fit = temp_fit if _n == ($test_start + `i')
	drop temp_fit 
}

gen adl_3sa_sqerror = (indpro_3step - adl_3sa_fit)^2
egen adl_3sa_mse = mean(adl_3sa_sqerror)
gen adl_3sa_rmse = sqrt(adl_3sa_mse)
di adl_3sa_rmse // Verify

* ADL with CLI+news [ADL(10,5,2)]

local horizon = 3
gen adln_3sa_fit = .
forval i = 0/$test_num_forecasts {
	local current_train_start = $train_start + `i' 
	local current_train_end = $train_end + `i' - `horizon' + 1
	// hac standard errors 
	qui newey indpro_3step L(3/12).indpro_1step L(3/7).cli L(3/4).news ///
		if _n >= `current_train_start' & _n <= `current_train_end', lag($hac_lag)
	predict temp_fit, xb 
	replace adln_3sa_fit = temp_fit if _n == ($test_start + `i')
	drop temp_fit 
}

gen adln_3sa_sqerror = (indpro_3step - adln_3sa_fit)^2
egen adln_3sa_mse = mean(adln_3sa_sqerror)
gen adln_3sa_rmse = sqrt(adln_3sa_mse)
di adln_3sa_rmse // Verify

/* Loss differentials */
* AR(12) & ADL(10,5) & ADL(10,5,2)
gen ld1_3sa = ar_3sa_sqerror - adl_3sa_sqerror
gen ld2_3sa = adl_3sa_sqerror - adln_3sa_sqerror

dfuller ld1_3sa, lag(12) // stationary at p=0.0329
dfuller ld2_3sa, lag(12) // not stationary at p=0.4733

ac ld1_3sa
ac ld2_3sa
corrgram ld1_3sa

// DM
dmariano indpro_3step ar_3sa_fit adl_3sa_fit if _n >= $test_start , crit(MSE) maxlag(4) kernel(bartlett)
dmariano indpro_3step adl_3sa_fit adln_3sa_fit if _n >= $test_start , crit(MSE) maxlag(4) kernel(bartlett)