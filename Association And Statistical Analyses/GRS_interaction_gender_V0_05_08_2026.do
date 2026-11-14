
**interaction analysis: GRS and gender
xi:logistic gbc grs_meta_SP_cat1#gender_s Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 
estimates store grs_1
xi:logistic gbc grs_meta_SP_cat1 gender_s Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 
estimates store grs_2
lrtest grs_1 grs_2
est tab,p(%12.10g)


***calculating excess relative risk 

*egen joint_exp_SP = group(grs_meta_SP_cat1 gender_s), label (this code is already exist that's why it has written as comment)
logistic gbc i.joint_exp_SP Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
nlcom exp(_b[4.joint_exp_SP]) - exp(_b[2.joint_exp_SP]) - exp(_b[3.joint_exp_SP]) + 1
nlcom exp(_b[6.joint_exp_SP]) - exp(_b[5.joint_exp_SP]) - exp(_b[3.joint_exp_SP]) + 1

