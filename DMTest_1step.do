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
	
***** 1 Step ahead forecast *****

/* Notations:

ar_1sa_rmse: AR 1 step ahead RMSE
adl_3sa_rmse: ADL with CLI 3 step ahead RMSE
adln_12sa_rmse: ADL with CLI + news 12 step ahead RMSE

*/

* AR(3) forecasts
gen ar_1sa_fit = .
forval i = 0/$test_num_forecasts {
	local current_train_start = $train_start + `i'
	local current_train_end = $train_end + `i'
	// hac standard errors 
	qui newey indpro_1step L(1/3).indpro_1step ///
		if _n >= `current_train_start' & _n <= `current_train_end', lag($hac_lag)
	predict temp_fit, xb 
	replace ar_1sa_fit = temp_fit if _n == ($test_start + `i')
	drop temp_fit 
}

gen ar_1sa_sqerror = (indpro_1step - ar_1sa_fit)^2
egen ar_1sa_mse = mean(ar_1sa_sqerror)
gen ar_1sa_rmse = sqrt(ar_1sa_mse)
di ar_1sa_rmse

* ADL with CLI VS benchmark AR(AIC) [ADL(10,8)]

// get point forecasts (rolling window)
gen adl_1sa_fit = .
forval i = 0/$test_num_forecasts {
	local current_train_start = $train_start + `i'
	local current_train_end = $train_end + `i'
	// hac standard errors 
	qui newey indpro_1step L(1/10).indpro_1step L(1/8).cli ///
		if _n >= `current_train_start' & _n <= `current_train_end', lag($hac_lag)
	predict temp_fit, xb 
	replace adl_1sa_fit = temp_fit if _n == ($test_start + `i')
	drop temp_fit 
}

gen adl_1sa_sqerror = (indpro_1step - adl_1sa_fit)^2
egen adl_1sa_mse = mean(adl_1sa_sqerror)
gen adl_1sa_rmse = sqrt(adl_1sa_mse)
di adl_1sa_rmse

* ADL with CLI+news [ADL(2,10,1)]

gen adln_1sa_fit = .
forval i = 0/$test_num_forecasts {
	local current_train_start = $train_start + `i'
	local current_train_end = $train_end + `i'
	// hac standard errors 
	qui newey indpro_1step L(1/2).indpro_1step L(1/10).cli L.news ///
		if _n >= `current_train_start' & _n <= `current_train_end', lag($hac_lag)
	predict temp_fit, xb 
	replace adln_1sa_fit = temp_fit if _n == ($test_start + `i')
	drop temp_fit 
}

gen adln_1sa_sqerror = (indpro_1step - adln_1sa_fit)^2
egen adln_1sa_mse = mean(adln_1sa_sqerror)
gen adln_1sa_rmse = sqrt(adln_1sa_mse)
di adln_1sa_rmse // Verify correct RMSE

/* Loss differentials */
* AR(3) & ADL(10,8) & ADL(2,10,1)
gen ld1_1sa = ar_1sa_sqerror - adl_1sa_sqerror
gen ld2_1sa = adl_1sa_sqerror - adln_1sa_sqerror

dfuller ld1_1sa // stationary at p=0.0000
dfuller ld2_1sa

ac ld1_1sa
ac ld2_1sa
// corrgram ld1_1sa 

// DM
dmariano indpro_1step ar3_1sa_fit adl_1sa_fit if _n >= $test_start , crit(MSE) maxlag(4) kernel(bartlett)
dmariano indpro_1step adl_1sa_fit adln_1sa_fit if _n >= $test_start , crit(MSE) maxlag(4) kernel(bartlett)
























