	*--------TWFE SANDBOX--------*
	*The goal of this code is to create simulated panel data and analyze it using the various TWFE estimators available.


	*--------Packages--------*
	*ssc install eventstudyinteract, replace	// Sun & Abraham (2021)
	*ssc install csdid, replace					// Callaway & Sant'Anna (2020)
	*ssc install did_multiplegt, replace		// Chaisemartin & d'Haultfoeuille (2020)
	*ssc install event_plot, replace			// Plots
	
	
	*--------Making the Data--------*
	clear

	local units = 30	// Can add more units or longer timeframe as desired
	local start = 1
	local end 	= 60

	local time = `end' - `start' + 1
	local obsv = `units' * `time'
	set obs `obsv'

	egen id	   = seq(), b(`time')  
	egen t 	   = seq(), f(`start') t(`end') 	

	sort  id t
	xtset id t

	set seed 23484			// Trying different seeds can produce different effect sizes and cohort timings

	gen Y 	   		= 0		// Outcome 	
	gen D 	   		= 0		// Treatment 
	gen cohort      = .  	// Cohort
	gen effect      = .		// Treatment effect size
	gen first_treat = .		// When the treatment happens for each cohort
	gen rel_time	= .     // time - first_treat (e.g., rel_time of treatment period = 0)

	levelsof id, local(lvls)
	foreach x of local lvls {			// Loop to assign cohorts
		local chrt = runiformint(0,5)	// Default is 5 treatment cohorts (plus never-treated)
		replace cohort = `chrt' if id==`x'
	}


	levelsof cohort , local(lvls)  
	foreach x of local lvls {			// Loop that assigns cohort treatment effects
		
		local eff = runiformint(2,10)
			replace effect = `eff' if cohort==`x'
				
		local timing = runiformint(`start',`end' + 20)	// Assigns cohort treatment timing
		replace first_treat = `timing' if cohort==`x'
		replace first_treat = . if first_treat > `end'
			replace D = 1 if cohort==`x' & t>= `timing' 
	}

	replace rel_time = t - first_treat
	replace first_treat = 0 if first_treat == . 
	
	replace Y = id + t + cond(D==1, effect * rel_time, 0) + rnormal() // DGP
	
	gen time_since_treat = rel_time if rel_time >-1
	replace time_since_treat = 0 if time_since_treat == .

	xtline Y, overlay legend(off)		// Visualize

	summ rel_time						// Creating relative time indicators
	local relmin = abs(r(min))
	local relmax = abs(r(max))

		// leads
		cap drop F_*
		forval x = 2/`relmin' {  // drop the first lead
			gen F_`x' = rel_time == -`x'
		}

		
		//lags
		cap drop L_*
		forval x = 0/`relmax' {
			gen L_`x' = rel_time ==  `x'
		}
		
	gen never_treat = first_treat == 0

	gen last_cohort = first_treat ==r(max) // Dummy for the latest- or never-treated cohort
	
	
	*--------SA--------*
	eventstudyinteract Y L_* F_*, vce(cluster id) absorb(id t) cohort(first_treat) control_cohort(never_treat)
	matrix sa_b = e(b_iw) // Store the estimates for later
	matrix sa_v = e(V_iw)
	
	*Plot the results:
	event_plot e(b_iw)#e(V_iw), default_look graph_opt(xtitle("Periods since the event") ytitle("Average effect") xlabel(-10(1)10) ///
			title("eventstudyinteract")) stub_lag(L_#) stub_lead(F_#) trimlag(10) trimlead(10) together
			
	*We can use the eventstudyweights package to diagnose how "bad" the naive TWFE estimation is. 	
	*Each column of this resulting table is a coefficient estimate corresponding to a lead / lag, and every row is a cohort-period. Higher number of negative weights is indicative of worse TWFE issues.
	eventstudyweights L_* F_*, absorb(i.t i.id) cohort(first_treat) rel_time(t)
	mat list e(weights)


	*--------CS--------*
	csdid Y, ivar(id) time(t) gvar(first_treat) method(drimp) notyet
	estat event, estore(cs) //Store the estimates for later
	
	*We can use the event option of the csdid_estat aggregator to get the event-study estimates
	csdid_estat event, window(-10 30) estore(cs) 
	
	event_plot cs, default_look graph_opt(xtitle("Periods since the event") ytitle("Average effect") ///
	title("csdid") xlabel(-10(1)30)) stub_lag(Tp#) stub_lead(Tm#) together	
	
	
	*--------CdH--------*
	did_multiplegt Y id t D, robust_dynamic dynamic(10) placebo(10) breps(20) cluster(id)
	matrix dcdh_b = e(didmgt_estimates) // Store the estimates for later
	matrix dcdh_v = e(didmgt_variances)

	
	*--------OLS--------*
	xtreg Y i.t F_* L_*, cluster(id) fe
	estimates store ols //Store estimates for later
	event_plot, default_look stub_lag(L_#) stub_lead(F_#) together graph_opt(xtitle("Days since the event") ytitle("OLS coefficients") xlabel(-14(1)5) ///
	title("OLS"))
	
	
	*-------Visualization-------*
	event_plot dcdh_b#dcdh_v cs sa_b#sa_v ols, ///
	stub_lag(Effect_# Tp# L_# L_#) stub_lead(Placebo_# Tm# F_# F_#) plottype(scatter) ciplottype(rcap) ///
	together perturb(-0.325(0.13)0.325) trimlead(5) trimlag(20) noautolegend ///
	graph_opt(title("Event study estimators in a simulated panel (30 units, 45 periods)", size(medlarge)) ///
		xtitle("Periods since the event") ytitle("Average causal effect") xlabel(-5(1)20) ylabel(-50(50)200) ///
		legend(order(1 "de Chaisemartin-D'Haultfoeuille" 3 "Callaway-Sant'Anna" 5 "Sun-Abraham" 7 "OLS") rows(2) region(style(none))) ///
		xline(-0.5, lcolor(gs8) lpattern(dash)) yline(0, lcolor(gs8)) graphregion(color(white)) bgcolor(white) ylabel(, angle(horizontal)) ///
	) ///
		lag_opt1(msymbol(+) color(cranberry)) lag_ci_opt1(color(cranberry)) ///
		lag_opt2(msymbol(Dh) color(navy)) lag_ci_opt2(color(navy)) ///
		lag_opt3(msymbol(Th) color(forest_green)) lag_ci_opt3(color(forest_green)) ///
		lag_opt4(msymbol(Sh) color(dkorange)) lag_ci_opt4(color(dkorange)) 
graph export "_estimators_example.png", replace
	