//////////////////////////////////////////////////////////////////////////////
// STATA code to generate Tables 5 and C.1
//////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// 1. Housekeeping
// Install packages
clear
ssc install inequal7
ssc install pshare

////////////////////////////////////////////////////////////////////////////////
// 1. Global macros for file paths
// Set file paths (RAWDIR MUST BE CHANGED)

global rawdir "C:\Users\glv2\Dropbox\HANK\ReplicationMaterial\STATA"
global datadir "${rawdir}\data"
global progdir "${rawdir}\prog"
global outputdir "${rawdir}\output"

//////////////////////////////////////////////////////////////////////////////
// 2. Read data

cd $datadir
use SCF_04_small.dta, clear

//////////////////////////////////////////////////////////////////////////////
// 3. Stats for Balance sheet of US households
//Table C.1

// housepos: real estate, durables=consumer durables,  liqpos=deposits, gb=government bonds, cb=corporate bonds,
// indirect_eq+direct_eq=corporate equity, netbus=equity in noncorp. business, houseneg=mortgage debt
// nrevccredit=nonrevolving consumer credit, revccredit=revolving consumer credit
// Entries in Table C.1 are obtained as "sum of weight x mean" from the output of the sum command below 

sum housepos durables liqpos gb cb indirect_eq direct_eq netbus houseneg nrevccredit revccredit [fw=fwgt], d

//////////////////////////////////////////////////////////////////////////////
// 5. Compute GINI coefficients and shares
//Table 5b

inequal7 netbrliq netbrilliq networth [fw=fwgt]

pshare estimate netbrliq netbrilliq networth [fw=fwgt], percentiles(25 50 90 99 99.9) 



