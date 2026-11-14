*****Gender stratified and whole datset: Manin table 2
local snps_all rs6733452_CA_A_AA_G rs15734_CA_A_AA_G rs67274221_CA_A_AA_G rs35604230_CA_A_AA_G rs60541816_CA_A_AA_G rs77155094_CA_T_AA_C rs2312959_CA_A_AA_C rs2131242countedallele_aalternat

levelsof gender_s, local(gender_levels)

foreach snp of local snps_all {

    display "================================="
    display "SNP: `snp'"
    display "================================="

    foreach s of local gender_levels {

        quietly count if `snp' != . & gender_s == `s'
        local N = r(N)
        local total_alleles = 2 * `N'

        quietly count if `snp' == 0 & gender_s == `s'
        local count_0 = r(N)

        quietly count if `snp' == 1 & gender_s == `s'
        local count_1 = r(N)

        local major_freq = (2*`count_0' + `count_1')/`total_alleles'
        local minor_freq = 1 - `major_freq'

        display "Gender = `s'"
        display "N = `N'"
        display "Major allele frequency = " round(`major_freq',0.0001)
        display "Minor allele frequency = " round(`minor_freq',0.0001)
        display ""
    }
}




**gender stratified separately for cases and controls

local snps_all rs6733452_CA_A_AA_G rs15734_CA_A_AA_G rs67274221_CA_A_AA_G rs35604230_CA_A_AA_G rs60541816_CA_A_AA_G rs77155094_CA_T_AA_C rs2312959_CA_A_AA_C rs2131242countedallele_aalternat


levelsof gender_s, local(gender_levels)
levelsof gbc, local(case_levels)

foreach snp of local snps_all {

    display "================================="
    display "SNP: `snp'"
    display "================================="

    foreach s of local gender_levels {

        foreach c of local case_levels {

            quietly count if `snp' != . & gender_s == `s' & gbc == `c'
            local N = r(N)

            if `N' > 0 {

                local total_alleles = 2 * `N'

                quietly count if `snp' == 0 & gender_s == `s' & gbc == `c'
                local count_0 = r(N)

                quietly count if `snp' == 1 & gender_s == `s' & gbc == `c'
                local count_1 = r(N)

                local major_freq = (2*`count_0' + `count_1')/`total_alleles'
                local minor_freq = 1 - `major_freq'

                display "Gender = `s'   Case/Control = `c'"
                display "N = `N'"
                display "Major allele frequency = " round(`major_freq',0.0001)
                display "Minor allele frequency = " round(`minor_freq',0.0001)
                display ""
            }
        }
    }
}
