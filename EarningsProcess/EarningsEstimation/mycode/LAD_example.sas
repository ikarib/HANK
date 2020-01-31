* Sample SAS program using the LAD;

libname source1 '/LADdata/data1;          * first 10% sample ;
libname source2 '/LADdata/data2;          * second 10% sample ;
libname Out '/LADuser/kais/data';          * user's directory ;

* This sample program's objective is to use the 20% LAD to find the summary statistics such as the median, skewness, kurtosis, and quantiles of earnings of Canadians that appear on their T4 slips, according to year (in this case, 1982 to 2016). Data for earnings are from the yearly LAD files.

* The first step is to create a datafile containing all the information that we need to produce our tables. This datafile will be called and will be saved in the 'out' directory. The Longitudinal Identifier Number (LIN__I) is used to merge the annual LAD datasets. ;

data out.SAOnt;
merge 
source1.lad2000(where=(prco_i2000 = 5) keep=lin__i  prco_i2000 saspyi2000 t4e__i2000)
source2.lad2000(where=(prco_i2000 = 5) keep=lin__i  prco_i2000 saspyi2000 t4e__i2000)
source1.lad2001(where=(prco_i2001 = 5) keep=lin__i prco_i2001 saspyi2001 t4e__i2001)
source2.lad2001(where=(prco_i2001 = 5) keep=lin__i  prco_i2001 saspyi2001 t4e__i2001)
source1.lad2002(where=(prco_i2002 = 5) keep=lin__i  prco_i2002 saspyi2002 t4e__i2002)
source2.lad2002(where=(prco_i2002 = 5) keep=lin__i  prco_i2002 saspyi2002 t4e__i2002)
source1.reg2002(keep=lin__i sxco_i flag_i2000-flag_i2002 wgt2_i) 
source2. reg2002(keep=lin__i sxco_i flag_i2000-flag_i2002 wgt2_i);

by lin__i ;

If flag_i2000=1 and flag_i2001=1 and flag_i2002=1; *person must be taxfiler in all 3 years;

* We create a flag variable that identifies the SA recipients for each year. The result is three variables, flag_sa2000, flag_sa2001 and flag_sa2002, taking a value of either 1 or 0.

If (t4e__i2000=0 and saspyi2000>0) then flag_sa2000 = 1 ; 
          else flag_sa2000 = 0 ;
if (t4e__i2001=0 and saspyi2001>0)  then flag_sa2001 = 1 ;
          else flag_sa2001 = 0 ; 
if (t4e__i2002=0 and saspyi2002>0) then flag_sa2002 = 1 ; 
          else flag_sa2002 = 0 ;

run ;

* The SAS 'freq' procedure is used to produce our tables. We would also need to make sure that confidentiality guidelines standards are respected. ;

proc freq data = out.SAOnt;

          tables sxco_i*flag_sa2000*flag_sa2001*flag_sa2002 /missing; 
          weight wgt2_i ;

run ;

* End of the sample program ;