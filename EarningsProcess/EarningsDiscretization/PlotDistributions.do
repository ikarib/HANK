# delim;
set more 1;
gr drop _all;

set scheme s1color;

*global OutputDir /Volumes/FILES/Projects/MollContinuousTime/Fortran/EarningsEstimation4/WithRhoAllMoments_v3;
global OutputDir /Volumes/FILES/Projects/MollContinuousTime/Fortran/DiscretizeEarnings2/2point_3_10;
local nbin1 = 50;
local nbin2 = 50;

clear;
insheet using $OutputDir/yannsim.txt, delim(",");

renvars v*, subst(v lyann);

gen yann = exp(lyann1);
gen d1lyann = lyann2-lyann1;
gen d5lyann = lyann5-lyann1;

sum yann, meanonly;
gen yann_relmean = yann/r(mean);
gen lyann_relmean = log(yann_relmean);

hist yann_relmean, bin(`nbin1') frac name(hist_yann)
	xtitle("Annual Earnings", size(large))
	ytitle("Fraction", size(large))
	color(red) lcolor(black) scheme(s1color);
gr export $OutputDir/hist_yann.pdf, replace;
gr export $OutputDir/hist_yann.eps, replace;

hist lyann_relmean, bin(`nbin1') frac name(hist_lyann)
	xtitle("Log Annual Earnings", size(large))
	ytitle("Fraction", size(large))
	color(red) lcolor(black) scheme(s1color);
gr export $OutputDir/hist_lyann.pdf, replace;
gr export $OutputDir/hist_lyann.eps, replace;

hist d1lyann, bin(`nbin2') frac name(hist_d1lyann)
	xtitle("1 Year Log Earnings Changes", size(large))
	ytitle("Fraction", size(large))
	color(red) lcolor(black) scheme(s1color);
gr export $OutputDir/hist_d1lyann.pdf, replace;
gr export $OutputDir/hist_d1lyann.eps, replace;

hist d5lyann, bin(`nbin2') frac name(hist_d5lyann)
	xtitle("5 Year Log Earnings Changes", size(large))
	ytitle("Fraction", size(large))
	color(red) lcolor(black) scheme(s1color);
gr export $OutputDir/hist_d5lyann.pdf, replace;
gr export $OutputDir/hist_d5lyann.eps, replace;

