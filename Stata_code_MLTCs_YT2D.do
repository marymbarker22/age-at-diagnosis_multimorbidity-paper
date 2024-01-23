

*****************
**Table 1********
*****************
cd "..."
cls
use "...", clear
mdesc exp age_index agegr sex eth imd bmi

baselinetable   							/*
*/	age_index(cts tab("p50 (p25, p75)")) 	/*
*/	agegr(cat countformat(%15.0fc))			/*
*/  sex(cat countformat(%15.0fc))  			/*
*/  eth(cat countformat(%15.0fc)) 			/*
*/  imd(cat countformat(%15.0fc)) 			/*
*/	bmi(cts tab("p50 (p25, p75)"))			/*
*/	anxiety_co(cat countformat(%15.0fc))  		/*
*/	depression_co(cat countformat(%15.0fc)) 	/* 
*/	hypertension_3(cat countformat(%15.0fc)) 	/*
*/	obesity(cat countformat(%15.0fc)) 			/*  
*/	, by(exp, total) reportmissing notable maxcatvalues(20)	/*
*/	exportexcel("...Table1", replace)

preserve
keep if bmi !=.
sum bmi, d
egen float bmig = cut(bmi), at(0 30 200) icodes label
label define bmig 0 `"<30"', modify
label define bmig 1 `"≥30"', modify
groups exp bmig, percent(exp) 
restore

******************
**Table S1********
******************
gen idn = 1

preserve
collapse (sum) anxiety_co depression_co obesity hypertension_3 idn, by(exp agegr)
foreach var of varlist anxiety_co depression_co obesity hypertension_3 {
	gen p`var' = `var'*100/idn
}
tempfile parts
save `parts', replace

collapse (sum) anxiety_co depression_co obesity hypertension_3 idn, by(exp)
foreach var of varlist anxiety_co depression_co obesity hypertension_3 {
	gen p`var' = `var'*100/idn
}
tempfile tots
save `tots', replace

use `parts', clear
append using `tots'
replace agegr = 7 if agegr == .
label define agegr_l 7 `"All"', add
gsort -exp agegr
renames anxiety_co depression_co obesity hypertension_3 panxiety_co pdepression_co pobesity phypertension_3 \ a d o h pa pd po ph
foreach var of varlist p* {
	tostring `var', replace format(%7.2f) force
}
foreach var of varlist a-h {
	tostring `var', replace format(%7.0fc) force
	replace `var' = `var' + " (" + p`var' + ")"
}
drop p*
order exp agegr idn, first
tostring idn, replace format(%7.0fc) force
reshape wide idn-h, i(agegr) j(exp)
order agegr idn1 a1 d1 h1 o1 idn0 a0 d0 h0 o0
gen cohort = 1
tempfile cohort1
save `cohort1', replace
restore

keep if bmi != .
collapse (sum) anxiety_co depression_co hypertension_3 idn, by(exp agegr)
foreach var of varlist anxiety_co depression_co hypertension_3 {
	gen p`var' = `var'*100/idn
}
tempfile parts
save `parts', replace

collapse (sum) anxiety_co depression_co hypertension_3 idn, by(exp)
foreach var of varlist anxiety_co depression_co hypertension_3 {
	gen p`var' = `var'*100/idn
}
tempfile tots
save `tots', replace

use `parts', clear
append using `tots'
replace agegr = 7 if agegr == .
label define agegr_l 7 `"All"', add
gsort -exp agegr
renames anxiety_co depression_co hypertension_3 panxiety_co pdepression_co phypertension_3 \ a d h pa pd ph
foreach var of varlist p* {
	tostring `var', replace format(%7.2f) force
}
foreach var of varlist a-h {
	tostring `var', replace format(%7.0fc) force
	replace `var' = `var' + " (" + p`var' + ")"
}
drop p*
order exp agegr idn, first
tostring idn, replace format(%7.0fc) force
reshape wide idn-h, i(agegr) j(exp)
order agegr idn1 a1 d1 h1 idn0 a0 d0 h0 
gen cohort = 2
tempfile cohort2
save `cohort2', replace

use `cohort1', clear
append using `cohort2'
save "...TableS1"


***********************************************
**Fig 1 [Risk Ratio for each condition]********
***********************************************
cd "..."
use "...", clear

local adj1 = "i.exp#i.agegr i.agegr i.sex i.eth i.imd"
local adj2 = "i.exp#i.agegr i.agegr i.sex i.eth i.imd bmi"

foreach did in anxiety_co depression_co obesity hypertension_3 {
	
	if "`did'" != "obesity" {
		
		forval k = 1/2 {
			qui glm `did' `adj`k'', link(log) family(binomial) 
			preserve
			qui parmest, fast eform
			qui split parm, p(#)
			qui keep if parm1 == "1.exp"
			qui sencode parm2, replace
			qui keep estimate parm2 min95 max95
			qui gen did = "`did'"
			qui gen model = `k'
			qui gen vars = "`adj`k''"
			qui tempfile rel_`did'_`k'
			qui save `rel_`did'_`k'', replace
			di "out = `did' | model = `k'"
			restore
			}
		}
	
		else {
			qui glm `did' `adj1', link(log) family(binomial)
			preserve
			qui parmest, fast eform
			qui split parm, p(#)
			qui keep if parm1 == "1.exp"
			qui sencode parm2, replace
			qui keep estimate parm2 min95 max95
			qui gen did = "`did'"
			qui gen model = 1
			qui gen vars = "`adj1'"
			qui tempfile rel_`did'_1
			qui save `rel_`did'_1', replace
			di "out = `did' | model = 1"
			restore
		}
}
	
clear
foreach did in anxiety_co depression_co obesity hypertension_3 {
	forval k = 1/2 {
		cap append using `rel_`did'_`k''
	}
}
qui order did model vars parm2 estimate min95 max95
renames did parm2 estimate model \ outcome ageg rr Model
label define parm2 1 `"16-27"', modify
label define parm2 2 `"28-31"', modify
label define parm2 3 `"32-35"', modify
label define parm2 4 `"36-39"', modify
label define parm2 5 `"40-43"', modify
label define parm2 6 `"44-47"', modify
label define parm2 7 `"48-50"', modify
split outcome, p("_")
replace outcome = outcome1
drop outcome1 outcome2
replace outcome = proper(outcome)
drop vars
sort Model outcome ageg
bys Model outcome (ageg): gen agegn = _n

twoway (scatter rr agegn if Model == 1, sort mcolor(red) msize(small)  msymbol(square)) 			///
       (scatter rr agegn if Model == 2, sort mcolor(blue) msize(small) msymbol(square)) 			///
       (rspike min95 max95 agegn if Model == 1, sort lcolor(red)  lwidth(thin)) 					///
       (rspike min95 max95 agegn if Model == 2, sort lcolor(blue) lwidth(thin))						///
       , yscale(log) xlabel(1 "16-27" 2 "28-31" 3 "32-35" 4 "36-39" 5 "40-43" 6 "44-47" 7 "48-50", grid)	///
       by(, legend(position(6))) legend(order(1 "Model 1" 2 "Model 2") rows(1) colgap(small)) 				///
       by(outcome, yrescale note("") scale(0.8)) subtitle(, size(medsmall)) 								///
       xtitle("Age group at diagnosis (years)") ytitle("Risk Ratio (T2D vs no T2D)")						///
       name("Fig1", replace) nodraw 
graph save "Fig1" "...Fig1.gph", replace


************************************************************
**Fig 2 [Predicted probabilitites of each condition]********
************************************************************
cd "..."
use "...", clear

local adj1 = "c.age_index##i.exp i.sex i.eth i.imd"
local adj2 = "c.age_index##i.exp i.sex i.eth i.imd bmi"

foreach did in anxiety_co depression_co obesity hypertension_3 {
	
	if "`did'" != "obesity" {
		
		forval k = 1/2 {
			qui glm `did' `adj`k'', link(log) family(binomial) 
			preserve
			qui tempfile tmp
			qui margins exp, at(age_index=(16(1)50)) saving ("`tmp'")
			qui use "`tmp'" , clear
			qui gen age_index = _at + 15
			qui gen per    = _margin*100
			qui gen per_lb = _ci_lb*100
			qui gen per_ub = _ci_ub*100
			qui keep per per_lb per_ub age_index _m1
			qui rename _m1 exp
			qui gen did = "`did'"
			qui gen model = `k'
			qui gen vars = "`adj`k''"
			qui tempfile abs_`did'_`k'
			qui save `abs_`did'_`k'', replace
			di "out = `did' | model = `k'"
			restore
			}
		}
		
		else {
			qui glm `did' `adj`k'', link(log) family(binomial) 
			preserve
			qui tempfile tmp
			qui margins exp, at(age_index=(16(1)50)) saving ("`tmp'")
			qui use "`tmp'" , clear
			qui gen age_index = _at + 15
			qui gen per    = _margin*100
			qui gen per_lb = _ci_lb*100
			qui gen per_ub = _ci_ub*100
			qui keep per per_lb per_ub age_index _m1
			qui rename _m1 exp
			qui gen did = "`did'"
			qui gen model = 1
			qui gen vars = "`adj1'"
			qui tempfile abs_`did'_1
			qui save `abs_`did'_1', replace
			di "out = `did' | model = 1"
			restore
		}
}		

clear
foreach did in anxiety_co depression_co obesity hypertension_3 {
	forval k = 1/2 {
		cap append using `abs_`did'_`k''
	}
}
qui order did model vars, first
renames did model exp \ outcome Model t2d
split outcome, p("_")
replace outcome = outcome1
drop outcome1 outcome2
replace outcome = proper(outcome)
drop vars

twoway (line per age_index if t2d == 0 & Model == 1, sort lcolor(red)  lpattern(shortdash)) 	///
       (line per age_index if t2d == 0 & Model == 2, sort lcolor(blue) lpattern(shortdash)) 	///
       (line per age_index if t2d == 1 & Model == 1, sort lcolor(red)  lpattern(solid)) 		///
       (line per age_index if t2d == 1 & Model == 2, sort lcolor(blue) lpattern(solid)) 		///
       (rarea per_lb per_ub age_index if t2d == 0 & Model == 1, color(red%30)  lwidth(none))	///
       (rarea per_lb per_ub age_index if t2d == 0 & Model == 2, color(blue%30) lwidth(none))	///
       (rarea per_lb per_ub age_index if t2d == 1 & Model == 1, color(red%30)  lwidth(none))	///
       (rarea per_lb per_ub age_index if t2d == 1 & Model == 2, color(blue%30) lwidth(none))	///
       , by(, legend(position(6))) xlabel(16(2)50, grid) ymtick(##2) 							///
       legend(order(												///
       1 "No type 2 diabetes, Model 1" 								///
       2 "No type 2 diabetes, Model 2" 								///
       3 "Type 2 diabetes, Model 1" 								///
       4 "Type 2 diabetes, Model 2") 								///
       rows(2) colgap(medsmall) rowgap(small) keygap(1))						///
       by(outcome, yrescale note("") scale(0.75)) subtitle(, size(medsmall))	///
       xtitle("Age at diagnosis (years)") ytitle("Probability (%)") name("Fig2", replace)
graph save "Fig2" "...Fig2.gph", replace


***************************************************************
**Fig 3/4 [Ratio and predicted number of comorbidities]********
***************************************************************

//RATIO COMORBIDITIES//
cd "..."
use "...", clear
gen n_comorb = anxiety_co + depression_co + obesity + hypertension_3
poisson n_comorb i.exp#i.agegr i.agegr i.sex i.eth i.imd, vce(robust)

preserve
parmest, fast eform
split parm, p(#)
keep if parm1 == "1.exp"
sencode parm2, replace
qui order parm2 estimate min95 max95
drop eq-parm1
rename parm2 ageg
label define parm2 1 `"16-27"', modify
label define parm2 2 `"28-31"', modify
label define parm2 3 `"32-35"', modify
label define parm2 4 `"36-39"', modify
label define parm2 5 `"40-43"', modify
label define parm2 6 `"44-47"', modify
label define parm2 7 `"48-50"', modify
gen Model = 1, before(ageg)

twoway (scatter estimate ageg, sort mcolor(red) msize(medsmall) msymbol(square)) 	///
       (rspike min95 max95 ageg, sort lcolor(red) lwidth(thin)) 					///
       , xlabel(1 "16-27" 2 "28-31" 3 "32-35" 4 "36-39" 5 "40-43" 6 "44-47" 7 "48-50", grid)							///
       scale(0.8) xtitle("Age group at diagnosis (years)") ytitle("Ratio of number of morbidities (T2D vs no T2D)")		///
       ymtick(##2) legend(off) name("Fig3", replace) nodraw
graph save "Fig3" "...Fig3.gph", replace
restore


//PREDICTED NUMBER OF COMORBIDITIES//
poisson n_comorb c.age_index##i.exp i.sex i.eth i.imd, vce(robust)
tempfile tmp
margins exp, at(age_index=(16(1)50)) saving ("`tmp'")
use "`tmp'" , clear
gen age_index = _at + 15
keep _margin _ci_lb _ci_ub _m1 age_index
renames _margin _ci_lb _ci_ub _m1 \ ncomb min95 max95 t2d 
gen Model = 1
order Model t2d age_index, first
foreach var of varlist ncomb min95 max95 {
	replace `var' = `var'*100
}  

twoway (line ncomb age_index if t2d == 0, sort lcolor(red)  lpattern(shortdash)) 	///
       (line ncomb age_index if t2d == 1, sort lcolor(red)  lpattern(solid)) 		///
       (rarea min95 max95 age_index if t2d == 0, color(red%30)  lwidth(none))		///
       (rarea min95 max95 age_index if t2d == 1, color(red%30)  lwidth(none))		///
       , legend(position(6) order(1 "No type 2 diabetes" 2 "Type 2 diabetes")		///
       rows(1) colgap(medsmall) rowgap(small) keygap(1)) scale(0.75)				///
       xlabel(16(2)50, grid) xtitle("Age at diagnosis (years)") ytitle("Number of comorbiditites (per 100 individuals)") ///
       ymtick(##2) ylabel(, format(%7.1f)) name("Fig4", replace) nodraw
graph save "Fig4" "...Fig4.gph", replace
graph close _all


***************************************************************
**Fig S3 [Distibution of number of diseases]*******************
***************************************************************
clear all
cd "..."
use "...", clear
gen n_comorb = anxiety_co + depression_co + obesity + hypertension_3

plotshares n_comorb if exp == 0, over(age_index)					///
	xtitle("Age at diagnosis (years)") ytitle("Proportion")			///
	legend(on rows(1) pos(6) order(1 "0" 2 "1" 3 "2" 4 "3" 5 "4")) 	///
	xlabel(16(2)50) ylabel(0(0.1)1, format(%3.1f)) scheme(rainbow)	///
	xsize(7) ysize(6) name("shr0", replace) clear title("No type 2 diabetes") output(frequency) nostack nodraw frame("exp0")
	
plotshares n_comorb if exp == 1, over(age_index)					///
	xtitle("Age at diagnosis (years)") ytitle("Proportion")			///
	legend(on rows(1) pos(6) order(1 "0" 2 "1" 3 "2" 4 "3" 5 "4")) 	///
	xlabel(16(2)50) ylabel(0(0.1)1, format(%3.1f)) scheme(rainbow)	///
	xsize(7) ysize(6) name("shr1", replace) clear title("Type 2 diabetes") output(frequency) nostack nodraw frame("exp1")
	
graph combine shr1 shr0, rows(1) ycommon scale(1.3) nocopies name("FigS3", replace) nodraw
graph save "FigS3" "...FigS3.gph", replace

forval k = 0/1 {
	frame change exp`k'
	drop x_val2-x_val5
	gen t2d = `k'
	renames x_val1 y_val1 y_val2 y_val3 y_val4 y_val5 \ age_index ncom0 ncom1 ncom2 ncom3 ncom4
	order t2d age_index, first
	tempfile exp`k'
	save `exp`k''
}
use `exp0', clear
append using `exp1'
gen tot = ncom0 + ncom1 + ncom2 + ncom3 + ncom4