import pandas as pd
import os
from pyliftover import LiftOver as LO
import time
import math

def liftover_method(Convert_From : int, Convert_To : int, ChainPath : str, InputBIM_DF : pd.DataFrame, OutputPath : str) -> None :
    NewBim = pd.DataFrame()
    OldPositions, OldChromosomes = list(), list()
    
    Chromosomes, BasePairs, SnpIDs = list(), list(), list()
    
    AdnanSami = LO(ChainPath)
    
    # for every chromosome, base pair position and SNP ID entry in the BIM File
    for chr, bp, snp in zip(list(InputBIM_DF["CHR"]), list(InputBIM_DF["BP"]), list(InputBIM_DF["SNP"])):
        OldChr = chr    # both the Old variables contain the chr and base pair as per the original BIM file build...whether it's 37 or 38
        OldBp = bp
        bp = int(bp)    # converting to int to prevent simple errors in the convert_coordinate method
        LiftOverReturn = []
        # As per Build GRCh38, chromosome codes are encoded as chr1, or chr13 etc along with different base pair positions.
        # Sometimes, the base pair positions are in GRCh38, but the chromosome codes are just 1 or 13 etc instead of chr1 or chr13 for example.
        # Alternatively, you may have a BIM file of GRCh37 build, but it'll have chr1 or chr13 instead of just 1 or 13 etc in the chromosome column.
        # this part of the code fixes that.
        if Convert_From == 38:          # the file should have chr1, chr2, chr19 etc...not just the number, i.e, not just 1, 2, 19 etc
            if "chr" not in str(chr):   # if there isn't a chr prefix
                chr = f"chr{chr}"       # then put the prefix
        elif Convert_From == 37:        # the file should NOT have chr1, chr2, chr19 etc...It should just have the number, i.e, 1, 2, 19 etc
            if "chr" in str(chr):       # if there is a chr prefix
                chr = str(chr).replace("chr", "")   # remove it
                chr = int(chr)          # and convert the variable into an int type
        # Now, to create the new Chr column term for the converted BIM File
        ChrTerm = ""
        if Convert_To == 38:                        # if converted to GRCh38, then the chr prefix must be attached in front of the chromosome number.
            ChrTerm = "chr" + str(chr)
        elif Convert_To == 37:                      # if converted to GRCh37, then the chr prefix must be removed from the front of the chromosome number.
            ChrTerm = str(chr).replace("chr", "")
            ChrTerm = int(ChrTerm)
        SnpIDs.append(snp)
        Chromosomes.append(ChrTerm)
        OldChromosomes.append(OldChr)
        OldPositions.append(OldBp)
        # by this point, we have saved the SNP id, the old chromosome number, old base pair position, new chromosome number. Only thing that's remaining is the new converted base pair position
        LiftOverReturn = AdnanSami.convert_coordinate(chr, bp)  # if converted successfully, this will return a list, with the second element [1] being the new base pair position
        if LiftOverReturn == [] or LiftOverReturn == None:      # else, this will return an empty list, denoting that the conversion failed.
            BasePairs.append(OldBp)                             # in such a case, we will simply use the pre-existing base pair position & pretend that it's the new one. It can be removed during matching in later stages
        else:
            TempBp = int(str(LiftOverReturn[0]).replace("[(", "").replace(")]", "").split(",")[1])  # extracting the new bp position
            BasePairs.append(TempBp)
    # the below DF contains the mapping of old / preexisting CHR:BP codes to new ones..Where there were errors, you'll find that the old and the new Base pair positions are the same
    NewBim = pd.DataFrame({"New_CHR" : Chromosomes, "New_BP" : BasePairs, "CHR" : OldChromosomes, "BP" : OldPositions, "SNP" : SnpIDs})
    # this one matches the input BIM file, along with the allele information, morgans column etc, with the above mapping dataframe. Thus we can now see the allele info as well for the new base pair positions
    # it essentially has the CHR and base pair positions for each SNP, both the old and converted ones (or old ones in case of errors in conversion)
    Merged = pd.merge(NewBim, InputBIM_DF, on=["CHR", "BP", "SNP"])
    ConvertedBIM_Chunk = pd.DataFrame()
    # this DF will now contains the new converted BIM file. Mind you, there will still be SNPs on the old Base pair positions, if they faced errors while switching the builds
    ConvertedBIM_Chunk.insert(0, "CHR", Merged["New_CHR"])          # new chromosome codes with the same old CHR column header
    ConvertedBIM_Chunk.insert(1, "SNP", Merged["SNP"])
    ConvertedBIM_Chunk.insert(2, "Morgans", Merged["Morgans"])
    ConvertedBIM_Chunk.insert(3, "BP", Merged["New_BP"])            # new base pair positions with the same old BP column header
    ConvertedBIM_Chunk.insert(4, "A1", Merged["A1"])
    ConvertedBIM_Chunk.insert(5, "A2", Merged["A2"])
    # the below file contains this so called converted BIM file
    ConvertedBIM_Chunk.to_csv(f"{OutputPath}/FULL_BIMfile_Lifted_from_{Convert_From}_to_{Convert_To}.csv", index=False)
    # the below file will contain the mapping information...as in
    NewBim.to_csv(f"{OutputPath}/Conversion_{Convert_From}-to-{Convert_To}_JUST_THE_SNPs_MAPPING_Including-Problematic-Sites.csv", index=False)
    # Now...we can make a removal list for problematic SNP locations.
    # Also, if we were to compare with another file, based on the same build AFTER conversion, we can probably get the locations of some of these SNPs that were not converted by the program
    # Firstly, getting a list of SNP IDs that were not converted.
    ProblemSNPs = NewBim[NewBim["New_BP"] == NewBim["BP"]]
    ProblemSNPs.to_csv(f"{OutputPath}/ALL_ProblematicSNPs_Not-Converted_{Convert_From}-to-{Convert_To}_JUST_THE_SNPs_MAPPING.csv", index=False)
    ProblemSNPs_number = len(ProblemSNPs)
    print(f"\n------------------\nThe original BIM file had\t{len(InputBIM_DF)}\tbase pair positions.\nUpon conversion, the BIM file has the following number of SNPs, that were not converted from\tGRCh{Convert_From}\tto\tGRCh{Convert_To}\t:\t{ProblemSNPs_number}\n------------------\n")
    ProblemSNPs["SNP"].to_csv(f"{OutputPath}/ALL_ProblematicSNPs_Not-Converted_{Convert_From}-to-{Convert_To}_JUST_THE_REMOVAL-LIST.txt", index=False, sep="\t", header=None)
    ProblemSNPs_List = list(ProblemSNPs["SNP"]) 
    # Secondly, intaking a second BIM file, if you want to...& trying to extract the positions of the above problematic SNPs from that BIM file, assuming the SNP IDs match
    if ProblemSNPs_number > 0:
        ModBIM = str(input(f"Enter\tY\tif you want to enter another BIM file, with the\tGRCh{Convert_To}\tbuild number,\nand preferably with the same type of SNP Ids, or rsIDs, to extract the base pair position of the above problematic SNPs\n\t\t:\t\t"))
        if ModBIM == "Y":
            ReferenceBIMPath = str(input(f"Enter the full path to such a BIM file\tBuild\t--\tGRCh{Convert_To}\t:\t"))
            RefBImDF = pd.read_csv(ReferenceBIMPath, sep="\s+", header=None, names = ["CHR", "SNP", "Morgans_Ref", "BP_Ref", "A1_Ref", "A2_Ref"])
            SNPs_In_RefBImDF = list(RefBImDF["SNP"])
            Modded_Chr, Modded_Bp, Modded_Snp, Modded_Morgans, Modded_A1, Modded_A2, StillProblemSNPs = list(), list(), list(), list(), list(), list(), list()
            Count, ProblematicSNPs_Found_Count, ProblematicSNPs_Fixed_Count, ProblematicSNPs_Remain_Count, Total_ProblemSites_Found, Total_ProblemSites_Fixed, Total_ProblemSites_Remain = 0, 0, 0, 0, 0, 0, 0
            Progress_PercentagePoints = 10
            PercentVerified = int(math.ceil(len(ConvertedBIM_Chunk) * (Progress_PercentagePoints / 100)))
            TotalPercentVerified = 0
            print(f"\n----------------------------\nThere are a total of\t{len(ConvertedBIM_Chunk)}\tbase pairs in the lifted BIM file\nThere are a total of\t{len(RefBImDF)}\tbase pairs in the reference BIM file\n-------------------------\n")
            # creating the selected df earlier itself. This will contains all the SNPs that are in the problematic list, that is, they did not get lifted over & they are also present in the reference dataset
            ReferenceBIM_Copy = RefBImDF.copy(deep=True)
            ReferenceBIM_Copy.columns = ["CHR", "SNP", "Morgans", "BP_Ref", "A1", "A2"]
            SelectedProblemSNPs = ReferenceBIM_Copy[ReferenceBIM_Copy["SNP"].isin(ProblemSNPs_List)]
            for row in ConvertedBIM_Chunk.itertuples():
                Count+=1
                if Count >= PercentVerified:
                    TotalPercentVerified += Progress_PercentagePoints
                    ProblematicSNPs_Remain_Count = ProblematicSNPs_Found_Count - ProblematicSNPs_Fixed_Count
                    print(f"{Count} Base pairs - \t{TotalPercentVerified}%\t - have been checked. {ProblematicSNPs_Found_Count} problematic SNPs found, out of which {ProblematicSNPs_Fixed_Count} have been fixed using the reference & {ProblematicSNPs_Remain_Count} problematic sites still remain")
                    Count = 0
                    Total_ProblemSites_Found += ProblematicSNPs_Found_Count
                    Total_ProblemSites_Fixed += ProblematicSNPs_Fixed_Count
                    Total_ProblemSites_Remain += ProblematicSNPs_Remain_Count
                    ProblematicSNPs_Found_Count = 0
                    ProblematicSNPs_Fixed_Count = 0
                    ProblematicSNPs_Remain_Count = 0
                Temp_Chr = row.CHR
                Temp_Snp = row.SNP
                Temp_Bp = row.BP
                Temp_Morgans = row.Morgans
                Temp_A1 = row.A1
                Temp_A2 = row.A2
                Modded_Chr.append(Temp_Chr)
                Modded_Snp.append(Temp_Snp)
                Modded_A1.append(Temp_A1)
                Modded_A2.append(Temp_A2)
                Modded_Morgans.append(Temp_Morgans)
                NewBP = Temp_Bp
                if Temp_Snp not in ProblemSNPs_List:    # The SNP that is being read is not problematic at all
                    NewBP = Temp_Bp
                    # Modded_Bp.append(Temp_Bp)
                    # Since the SNPs are not problematic, we skip it & simply add the SNP to the list.
                elif Temp_Snp in ProblemSNPs_List and Temp_Snp not in SNPs_In_RefBImDF:
                    ProblematicSNPs_Found_Count += 1
                    NewBP = Temp_Bp
                    # ProblematicSNPs_Remain_Count += 1
                    # Modded_Bp.append(NewBP)
                    StillProblemSNPs.append(Temp_Snp)
                elif Temp_Snp in ProblemSNPs_List and Temp_Snp in SNPs_In_RefBImDF: 
                    # these SNPs are in the problematic list, that is, they did not get lifted over & they are also present in the reference dataset
                    ProblematicSNPs_Found_Count += 1    
                    sel = SelectedProblemSNPs[(SelectedProblemSNPs["SNP"] == Temp_Snp)]
                    sel_bp = list(sel["BP_Ref"])[0]
                    sel_a1 = list(sel["A1"])[0]
                    sel_chr = list(sel["CHR"])[0]
                    if len(sel)==1:
                        if not (sel_bp == Temp_Bp) and ( sel_a1 == Temp_A1) and (sel_chr == Temp_Chr):
                            # sel = RefBImDF[(RefBImDF["CHR"] == Temp_Chr) & (RefBImDF["SNP"] == Temp_Snp) & ~(RefBImDF["BP_Ref"] == Temp_Bp) & (RefBImDF["A1_Ref"] == Temp_A1) & (RefBImDF["A2_Ref"] == Temp_A2)]
                            ProblematicSNPs_Fixed_Count += 1
                            NewBP = int(sel_bp)
                            # Modded_Bp.append(NewBP)
                    else:
                        # ProblematicSNPs_Remain_Count += 1
                        NewBP = Temp_Bp
                        # Modded_Bp.append(NewBP)
                        StillProblemSNPs.append(Temp_Snp)
                else:
                    print("\n\n\t\tFATAL ERROR in Conditional Segment, Inner Loop\n\n\n")
                    exit()
                Modded_Bp.append(NewBP)
                # if Temp_Snp in ProblemSNPs_List and Temp_Snp in SNPs_In_RefBImDF:
                #     ProblematicSNPs_Found_Count += 1    
                #     sel = RefBImDF[(RefBImDF["CHR"] == Temp_Chr) & (RefBImDF["SNP"] == Temp_Snp) & ~(RefBImDF["BP_Ref"] == Temp_Bp) & (RefBImDF["A1_Ref"] == Temp_A1) & (RefBImDF["A2_Ref"] == Temp_A2)]
                #     if len(sel)==1:
                #         ProblematicSNPs_Fixed_Count += 1
                #         NewBP = int(sel["BP_Ref"])
                #         Modded_Bp.append(NewBP)
                #     else:
                #         ProblematicSNPs_Remain_Count += 1
                #         NewBP = Temp_Bp
                #         Modded_Bp.append(NewBP)
                #         StillProblemSNPs.append(row.SNP)
                # elif Temp_Snp in ProblemSNPs_List and Temp_Snp not in SNPs_In_RefBImDF:
                #     ProblematicSNPs_Remain_Count += 1
                #     NewBP = Temp_Bp
                #     Modded_Bp.append(NewBP)
                #     StillProblemSNPs.append(row.SNP)
                # else:
                #     Modded_Bp.append(Temp_Bp)
                #     continue
            print(f"A total of {Total_ProblemSites_Found} probematic sites were found, with Base pair positions still as per GRCh{Convert_From}. {Total_ProblemSites_Fixed} sites were fixed using the reference BIM file & {Total_ProblemSites_Remain} sites still remain problematic\n")
            print(f"Length of Modded Chr, BP, SNP, A1, A2\t:\t{len(Modded_Chr)}\t{len(Modded_Bp)}\t{len(Modded_Snp)}\t{len(Modded_A1)}\t{len(Modded_A2)}")
            ModdedBIM = pd.DataFrame({"CHR":Modded_Chr, "SNP":Modded_Snp, "Morgans" : Modded_Morgans,"BP":Modded_Bp, "A1":Modded_A1, "A2":Modded_A2})
            # ModdedBIM.insert(2, "Morgans", ModdedBIM["SNP"])
            # ModdedBIM["Morgans"] = [0 for i in range(0, len(ModdedBIM))]
            # ModdedBIM = pd.merge(ConvertedBIM_Chunk, RefBImDF[["CHR", "SNP", "BP_Ref"]], on=["CHR", "SNP"], how="left")
            # ModdedBIM.loc[(~(ModdedBIM["BP"] == ModdedBIM["BP_Ref"])), "BP"] = ModdedBIM["BP_Ref"]
            print(ModdedBIM.head(5))
            print(ModdedBIM.tail(5))
            ModdedBIM.to_csv(f"{OutputPath}/FULL_BIMfile_Lifted_from_{Convert_From}_to_{Convert_To}_WithProblemSNPs-Modded.csv", index=False)
            ModdedBIM.to_csv(f"{OutputPath}/FULL_BIMfile_Lifted_from_{Convert_From}_to_{Convert_To}_WithProblemSNPs-Modded.bim", index=False, header=None, sep="\t")
            # New_ProblemSNPs = ModdedBIM[(~(ModdedBIM["BP"] == ModdedBIM["BP_Ref"]) & (ModdedBIM["CHR"] == ModdedBIM["CHR"]) & (ModdedBIM["SNP"] == ModdedBIM["SNP"]))]
            New_ProblemSNPs = ModdedBIM[ModdedBIM["SNP"].isin(StillProblemSNPs)]
            New_ProblemSNPs.to_csv(f"{OutputPath}/Remaining_ProblematicSNPs_Not-Converted_{Convert_From}-to-{Convert_To}_AfterModding.csv", index=False)
            New_ProblemSNPs_number = len(New_ProblemSNPs)
            print(f"\n----------------------------\nThere are a total of\t{len(ModdedBIM)}\tbase pair positions.\nUpon modding, the BIM file has the following number of SNPs, that were still not converted from\tGRCh{Convert_From}\tto\tGRCh{Convert_To}\t:\t{New_ProblemSNPs_number}")
            pd.DataFrame({"Still_ProblematicSNPs" : StillProblemSNPs}).to_csv(f"{OutputPath}/Remaining_ProblematicSNPs_Not-Converted_{Convert_From}-to-{Convert_To}_AfterModding_REMOVAL-LIST.txt", header=None, index=False, sep="\t")
        else:
            pass
    

bimFile = str(input("Enter the FULL path of the BIM file, with the \t.bim\t extension\nThe BED & FAM files should be in the same location, with the same name\n\n\t\tEnter:\t"))
c_from = int(input("Enter\t37\tif converting from GRCh37\t\tOR\nEnter\t38\tif converting from GRCh38\n\n\t\tEnter\t:\t"))
c_to = 37 if c_from==38 else 38
ChainFilePath = str(input("Enter the full path to the CORRESPONDING liftover chain path file with the extensions\tover.chain.gz\tfor this liftover\t:\t"))
OutputPath = str(input("Enter the path to the output folder where the lifted BIM file etc will be saved\t:\t"))

BimDF = pd.read_csv(bimFile, sep="\s+", header=None, names = ["CHR", "SNP", "Morgans", "BP", "A1", "A2"])

liftover_method(Convert_From=c_from, Convert_To=c_to, ChainPath=ChainFilePath, InputBIM_DF=BimDF, OutputPath=OutputPath)
