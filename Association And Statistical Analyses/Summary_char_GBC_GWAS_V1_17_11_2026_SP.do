*Summary statistics
*Gender
bysort ForManuscript_Phase: ta gender_s gbc ,col m chi
ta gender_s gbc ,col m chi


*Age
bysort ForManuscript_Phase: ta Age_cat_RM_09112025 gbc,col m chi
ta Age_cat_RM_09112025 gbc ,col m chi

*Mean and SD controls
bysort ForManuscript_Phase: sum Age if gbc==0
*Mean and SD cases
bysort ForManuscript_Phase: sum Age if gbc==1
bysort ForManuscript_Phase : ttest Age,by(gbc)

bysort gbc: sum Age
ttest Age,by(gbc)


*Current residense
bysort ForManuscript_Phase: ta Cur_res_Zone gbc,col m chi
ta Cur_res_Zone gbc,col m chi


*birth risk region
bysort ForManuscript_Phase: ta Birth_res_risk_zone gbc,col m chi
ta Birth_res_risk_zone gbc,col m chi


*Education self
bysort ForManuscript_Phase: ta education_self gbc,col m chi
ta education_self gbc,col m chi


*gallstone
bysort ForManuscript_Phase: ta stone_in_GB_YN_original_s gbc,col m chi
ta stone_in_GB_YN_original_s gbc,col m chi

