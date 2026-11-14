**Cleaning and coding: GBC GWAS

*recode to bidchip
encode Infinium, gen(Infinium_code)
recode Infinium_code (3=1 "Humanomniexpress-24 v1.1")(2=2 "GSA-24 v3.0")(1=3 "Confluence"),gen(bidchip_SP)


*rs60541816_CA_A_AA_G
gen nrs60541816_CA_A_AA_G_SP=.
replace nrs60541816_CA_A_AA_G_SP=1 if rs60541816_CA_A_AA_G==1
replace nrs60541816_CA_A_AA_G_SP=2 if rs60541816_CA_A_AA_G==0
replace nrs60541816_CA_A_AA_G_SP=0 if rs60541816_CA_A_AA_G==2
tab nrs60541816_CA_A_AA_G_SP


*rs2131242countedallele_aalternat
gen nrs2131242_SP =.
replace nrs2131242_SP=1 if rs2131242countedallele_aalternat==1
replace nrs2131242_SP=2 if rs2131242countedallele_aalternat==0
replace nrs2131242_SP=0 if rs2131242countedallele_aalternat==2
tab nrs2131242_SP


*Calculate GRS
gen grs_meta_SP = rs6733452_CA_A_AA_G + rs15734_CA_A_AA_G + rs67274221_CA_A_AA_G + rs35604230_CA_A_AA_G + rs77155094_CA_T_AA_C + rs2312959_CA_A_AA_C + nrs60541816_CA_A_AA_G_SP + nrs2131242_SP


*Categories of GRS based on controls
pctile grs_meta_SP_cont = grs_meta_SP if gbc==0, nq(3)
xtile grs_meta_SP_cat = grs_meta_SP, cutpoints(grs_meta_SP_cont)
tabstat grs_meta_SP ,by(grs_meta_SP_cat) statistics (min max mean median)

