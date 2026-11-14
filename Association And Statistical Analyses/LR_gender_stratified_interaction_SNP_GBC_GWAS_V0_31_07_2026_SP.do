*logistic regression: stratified by gender

bysort gender_s: ta gbc

*rs6733452_CA_A_AA_G
xi:logistic gbc rs6733452_CA_A_AA_G Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 if gender_s==0
est tab,p(%12.10g)
xi:logistic gbc rs6733452_CA_A_AA_G Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 if gender_s==1
est tab,p(%12.10g)

*rs15734_CA_A_AA_G
xi:logistic gbc rs15734_CA_A_AA_G Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 if gender_s==0
est tab,p(%12.10g)
xi:logistic gbc rs15734_CA_A_AA_G Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 if gender_s==1
est tab,p(%12.10g)

*rs67274221_CA_A_AA_G
xi:logistic gbc rs67274221_CA_A_AA_G Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 if gender_s==0
est tab,p(%12.10g)
xi:logistic gbc rs67274221_CA_A_AA_G Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 if gender_s==1
est tab,p(%12.10g)

*rs35604230_CA_A_AA_G
xi:logistic gbc rs35604230_CA_A_AA_G Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 if gender_s==0
est tab,p(%12.10g)
xi:logistic gbc rs35604230_CA_A_AA_G Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 if gender_s==1
est tab,p(%12.10g)

*rs60541816_CA_A_AA_G
xi:logistic gbc rs60541816_CA_A_AA_G Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 if gender_s==0
est tab,p(%12.10g)
xi:logistic gbc rs60541816_CA_A_AA_G Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 if gender_s==1
est tab,p(%12.10g)

*rs77155094_CA_T_AA_C
xi:logistic gbc rs77155094_CA_T_AA_C Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 if gender_s==0
est tab,p(%12.10g)
xi:logistic gbc rs77155094_CA_T_AA_C Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 if gender_s==1
est tab,p(%12.10g)

*rs2312959_CA_A_AA_C
xi:logistic gbc rs2312959_CA_A_AA_C Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 if gender_s==0
est tab,p(%12.10g)
xi:logistic gbc rs2312959_CA_A_AA_C Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 if gender_s==1
est tab,p(%12.10g)

*rs2131242countedallele_aalternat
xi:logistic gbc rs2131242countedallele_aalternat Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 if gender_s==0
est tab,p(%12.10g)
xi:logistic gbc rs2131242countedallele_aalternat Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 if gender_s==1
est tab,p(%12.10g)




*Interaction analysis
* wald test used for only this 'rs' because observations are differ
xi:logistic gbc rs6733452_CA_A_AA_G##gender_s Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
testparm rs6733452_CA_A_AA_G#gender_s

xi:logistic gbc rs15734_CA_A_AA_G gender_s Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
estimates store a
xi:logistic gbc rs15734_CA_A_AA_G##gender_s Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
estimates store b
lrtest a b

xi:logistic gbc rs67274221_CA_A_AA_G gender_s Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
estimates store c
xi:logistic gbc rs67274221_CA_A_AA_G##gender_s Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
estimates store d
lrtest c d

xi:logistic gbc rs35604230_CA_A_AA_G gender_s Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
estimates store e
xi:logistic gbc rs35604230_CA_A_AA_G##gender_s Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
estimates store f
lrtest e f

xi:logistic gbc rs60541816_CA_A_AA_G gender_s Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
estimates store g
xi:logistic gbc rs60541816_CA_A_AA_G##gender_s Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
estimates store h
lrtest g h

xi:logistic gbc rs77155094_CA_T_AA_C gender_s Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
estimates store i
xi:logistic gbc rs77155094_CA_T_AA_C##gender_s Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
estimates store j
lrtest i j

xi:logistic gbc rs2312959_CA_A_AA_C gender_s Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
estimates store k
xi:logistic gbc rs2312959_CA_A_AA_C##gender_s Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
estimates store l
lrtest k l

xi:logistic gbc rs2131242countedallele_aalternat gender_s Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
estimates store m
xi:logistic gbc rs2131242countedallele_aalternat##gender_s Age PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
estimates store n
lrtest m n












