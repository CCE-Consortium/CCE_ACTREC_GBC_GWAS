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


*BMI
bysort ForManuscript_Phase: ta BMI_cat_international gbc,col m chi
ta BMI_cat_india gbc,col m chi

*Mean and SD controls
bysort ForManuscript_Phase: sum BMI if gbc==0
*Mean and SD cases
bysort ForManuscript_Phase: sum BMI if gbc==1
bysort ForManuscript_Phase : ttest BMI ,by(gbc)

bysort gbc: sum BMI 
ttest BMI ,by(gbc)


*tobacco chewing
bysort ForManuscript_Phase: ta Tobacco_chewing_YN gbc,col m chi
ta Tobacco_chewing_YN gbc,col m chi

*Mean and SD controls
bysort ForManuscript_Phase: sum Any_tobacco_chewing_duration if gbc==0 & Tobacco_chewing_YN==1
*Mean and SD cases
bysort ForManuscript_Phase: sum Any_tobacco_chewing_duration if gbc==1 & Tobacco_chewing_YN==1
bysort ForManuscript_Phase : ttest Any_tobacco_chewing_duration ,by(gbc)

bysort gbc: sum Any_tobacco_chewing_duration if Tobacco_chewing_YN==1
ttest Any_tobacco_chewing_duration ,by(gbc)


*tobacco smoking
bysort ForManuscript_Phase: ta Smoking_YN gbc,col m chi
ta Smoking_YN gbc,col m chi

*Mean and SD controls
bysort ForManuscript_Phase: sum Any_smoke_duration if gbc==0 & Smoking_YN==1
*Mean and SD cases
bysort ForManuscript_Phase: sum Any_smoke_duration if gbc==1 & Smoking_YN==1
bysort ForManuscript_Phase : ttest Any_smoke_duration ,by(gbc)

bysort gbc: sum Any_smoke_duration if Smoking_YN==1
ttest Any_smoke_duration ,by(gbc)


*alcohol consumption
bysort ForManuscript_Phase: ta Alcohol_YN gbc,col m chi
ta Alcohol_YN gbc,col m chi

*Mean and SD controls
bysort ForManuscript_Phase: sum Any_alco_GramPD if gbc==0 & Alcohol_YN==1
*Mean and SD cases
bysort ForManuscript_Phase: sum Any_alco_GramPD if gbc==1 & Alcohol_YN==1
bysort ForManuscript_Phase : ttest Any_alco_GramPD ,by(gbc)

bysort gbc: sum Any_alco_GramPD if Alcohol_YN==1
ttest Any_alco_GramPD ,by(gbc)


*musturd oil consumption
bysort ForManuscript_Phase: ta Oil_YN gbc,col m chi
ta Oil_YN gbc,col m chi

*Mean and SD controls
bysort ForManuscript_Phase: sum con_per_person_merged if gbc==0 & Oil_YN==1
*Mean and SD cases
bysort ForManuscript_Phase: sum con_per_person_merged if gbc==1 & Oil_YN==1
bysort ForManuscript_Phase : ttest con_per_person_merged ,by(gbc)

bysort gbc: sum con_per_person_merged if Oil_YN==1
ttest con_per_person_merged ,by(gbc)


*gallstone
bysort ForManuscript_Phase: ta stone_in_GB_YN_original_s gbc,col m chi
ta stone_in_GB_YN_original_s gbc,col m chi


*gallstone:before 1 year
bysort ForManuscript_Phase: ta stone_in_GB_YN_RM_09122025 gbc,col m chi
ta stone_in_GB_YN_RM_09122025 gbc,col m chi


*family hist of cancer First degree
bysort ForManuscript_Phase: ta CA_GB_YN gbc,col m chi
ta CA_GB_YN gbc,col m chi
