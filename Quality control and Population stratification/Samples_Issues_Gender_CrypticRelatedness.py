import threading
import pandas as pd
import time
import os

# This program will deal with aspects of samples being somehow (cryptically / mysteriously) related to each other, despite not being from the same family at all

#################################################
    # First part    :   Creating the gender check, call rates & cryptic relatedness file
    
InputBEDFilePath = str(input("Enter the full path to the BED fileset in which QC needs to be done\t:\t"))
PlinkPath = str(input("Enter the full path to PLINK 1.9 executable\t:\t"))
CrypticRelatedness_MinThreshold = str(input("Enter the float value of a percentage of PI_HAT, beyond which two samples will be considered related to each other\t:\t"))
# Extracting the file path and the file name from the above given path
FileName, FilePath = "", ""
PathComponents = InputBEDFilePath.split("/")[1:]
for pc in PathComponents[0:-1]:
    FilePath += f"/{pc}"
FileName = str(PathComponents[len(PathComponents)-1])

# now to create commands for checking for gender issues and cryptic relatedness, which will also require the call rate of samples (using the missing command)
CmdsCommonString = f"{PlinkPath} --bfile {InputBEDFilePath}"
CmdsCommonString_2 = f"--out {FilePath}/{FileName}"
SexCheckCommand, CryptRelCommand, CallRateCommand = f"{CmdsCommonString} --check-sex {CmdsCommonString_2}_Gender", \
    f"{CmdsCommonString} --genome --max 1 -min {CrypticRelatedness_MinThreshold} {CmdsCommonString_2}_CryptRel", \
        f"{CmdsCommonString} --missing {CmdsCommonString_2}_CallRate"

# now to execute all three commands
def command_execute(CommandString : str) -> None:
    os.system(CommandString)
    
CommandList = [SexCheckCommand, CryptRelCommand, CallRateCommand]
CommandThreadsList = list()
for cmd, threadname in zip(CommandList, ["SexCheck_Thread", "CryptRel_Thread", "CallRate_Thread"]):
    CommandThread = threading.Thread(target=command_execute, name=threadname, kwargs={"CommandString": cmd})
    CommandThreadsList.append(CommandThread)
    
for t in CommandThreadsList:
    t.start()
    
for t in CommandThreadsList:
    t.join()

#################################################
#   By this point, we have created all the files needed to proceed with the analysis of which samples to remove and why

# getting the IDs of all those samples, with gender issues
GenderDF = pd.read_csv(f"{FilePath}/{FileName}_Gender.sexcheck", sep="\s+")
GenderDF["FID"] = GenderDF["FID"].astype("str")
GenderDF["IID"] = GenderDF["IID"].astype("str")
# ["FID"] = ["FID"].astype("str")
GenderProblems = GenderDF[GenderDF["STATUS"] == "PROBLEM"]      # this are all the samples with some sort of gender issues
GenderProblems.columns = ["FID", "IID", "Sex_GivenData", "Sex_GeneticData", "STATUS", "F"]
GenderProblems_Full_Len = len(GenderProblems)
GenderProblems.to_csv(f"{FilePath}/{FileName}_Gender_Problems_Full_{GenderProblems_Full_Len}-Samples.csv", index=False) # this is the full file of the gender issues check
# now we need to separate the gender issues into two types : One where plink was unable to guess the gender of the sample & another where plink was able to guess

# selecting all those gender issues samples where plink was unable to identify the gender based on DNA data.
# these surely need to be removed.
GenderProblems_Undetermined = GenderProblems[GenderProblems["Sex_GeneticData"] == 0]
GenderProblems_Undetermined_Len = len(GenderProblems_Undetermined)
GenderProblems_Undetermined.to_csv(f"{FilePath}/{FileName}_Gender_Problems_Undetermined_{GenderProblems_Undetermined_Len}-Samples.csv", index=False)

# making the removal list...of samples that surely need to be removed from the data
GenderRemoval = GenderProblems_Undetermined[["FID", "IID"]] # contains IDs that need to be verified and then removed
GenderRemoval_Len = len(GenderRemoval)
GenderRemoval.to_csv(f"{FilePath}/{FileName}_Gender_Problems_RemovalList_{GenderRemoval_Len}-Samples.txt", index=False, sep="\t")

# selecting all those samples where plink was able to guess the gender of the sample.
# The questionairre data of these samples will be analyzed to determine their gender and thus, avoid, as much as possible, to remove these additional samples from the data
GenderProblems_Guessed = GenderProblems[~ (GenderProblems["Sex_GeneticData"] == 0)]
GenderProblems_Guessed_Len = len(GenderProblems_Guessed)
GenderProblems_Guessed.to_csv(f"{FilePath}/{FileName}_Gender_Problems_Guessed_{GenderProblems_Guessed_Len}-Samples.csv", index=False)

######################
# Now for cryptic relatedness, we need to keep ready the call rates, age, gender and case control data of the samples in hand, for comparison purposes

SamplesCallRate = pd.read_csv(f"{FilePath}/{FileName}_CallRate.imiss", sep="\s+")
SamplesCallRate["FID"] = SamplesCallRate["FID"].astype("str")
SamplesCallRate["IID"] = SamplesCallRate["IID"].astype("str")
SamplesCallRate["Number_Of_Variants_Present"] = SamplesCallRate["N_GENO"] - SamplesCallRate["N_MISS"]
SamplesCallRate["Call_Rate"] = SamplesCallRate["Number_Of_Variants_Present"] / SamplesCallRate["N_GENO"]
SamplesCallRate.to_csv(f"{FilePath}/{FileName}_CallRate.imiss.csv", index=False)

# here we have all the call rates information of each sample in the data.
# next, we need their age, gender and phenotype info. The FAM file will contain the pheno and gender, but age needs to be given via some other means.
    # THEREFORE, instead of suffering through this, the decision to provide the pedigree file, such that it ALSO CONTAINS THE AGE has been taken and been made mandatory

PedigreeFilePath = str(input("Enter the full path of the CSV pedigree file containing the columns\t:\tFID, IID, Affected_Status, Gender, Age\t:\t"))
PedigreeFile = pd.read_csv(PedigreeFilePath)
PedigreeFile["FID"] = PedigreeFile["FID"].astype("str")
PedigreeFile["IID"] = PedigreeFile["IID"].astype("str")

######################
# now for cryptic relatedness, AND also to filter the call rates and the pedigree files and only to keep the IDs found in the .genome file below.

CryptRelDF = pd.read_csv(f"{FilePath}/{FileName}_CryptRel.genome", sep="\s+")
CryptRelDF["FID1"] = CryptRelDF["FID1"].astype("str")
CryptRelDF["FID2"] = CryptRelDF["FID2"].astype("str")
CryptRelDF["IID1"] = CryptRelDF["IID1"].astype("str")
CryptRelDF["IID2"] = CryptRelDF["IID2"].astype("str")

# now we get the unique IDs of all the samples that are present in .genome crypt rel file
AllContaminatedIDs_Temp = pd.Series()
AllContaminatedIDs = list()
AllContaminatedIDs.extend(CryptRelDF["FID1"].unique())
AllContaminatedIDs.extend(CryptRelDF["FID2"].unique())
AllContaminatedIDs_Temp = pd.Series(AllContaminatedIDs)
AllContaminatedIDs = pd.Series(AllContaminatedIDs_Temp.unique())
AllContaminatedIDs.name = "FID"

# now using the above list, we can filter the call rates and the pedigree file

SamplesCallRate = pd.merge(SamplesCallRate, AllContaminatedIDs, on="FID")
PedigreeFile = pd.merge(PedigreeFile, AllContaminatedIDs, on="FID")
# Now we start with the actual filtering and selection of contaminated samples


###################################################

# Firstly, the most simple one  :   Removing planned duplicates in a proper manner

# Getting hold of all the planned duplicates
PlannedDuplicates = pd.DataFrame()

                    # select all rows where either FID1 or FID2 contain _R
                    # AND
                    # where FID1 = FID2 with their _R removed
PlannedDuplicates = CryptRelDF[ ((CryptRelDF["FID1"].str.contains("_R")) | (CryptRelDF["FID2"].str.contains("_R")))\
            & \
           (CryptRelDF["FID1"].str.replace("_R","") == CryptRelDF["FID2"].str.replace("_R",""))\
           ]

# Now for the samples that have been contaminated with another sample (no planned duplicates here)
CryptRelSamples = pd.DataFrame()

                    # select all rows where neither FID1 nor FID2 contain _R
                    # AND
                    # where FID1 != FID2
CryptRelSamples = CryptRelDF[ ~ ((CryptRelDF["FID1"].str.contains("_R")) | (CryptRelDF["FID2"].str.contains("_R")))\
            &
            ~ (CryptRelDF["FID1"] == CryptRelDF["FID2"])\
    ]

# Now for the samples that are contaminated with another sample's _R ID. These kind of samples all need to be removed.
Unplanned_Duplicates = pd.DataFrame()

                    # select all rows where either FID1 or FID2 contains _R
                    # AND
                    # where FID1 != FID2, when their _R is removed
Unplanned_Duplicates = CryptRelDF[ ((CryptRelDF["FID1"].str.contains("_R")) | (CryptRelDF["FID2"].str.contains("_R")))\
            & \
           (~(CryptRelDF["FID1"].str.replace("_R","") == CryptRelDF["FID2"].str.replace("_R","")))\
           ]

####################################################

# Now to decide which samples to remove from each of the above 3 dataframes

# For planned duplicates, we simply compare their call rates and remove the sample with the lower call rate.
PlannedDuplicates["ToBeRemoved"] = ""       # the sample IDs that need to be removed will be addded to this column
PlannedDuplicates["RemovalReason"] = ""     # the reason for removal will be mentioned
PlannedDuplicates["CallRate_FID1"] = 0.0    # this is the placeholder column for mentioning the call rate of the sample..
PlannedDuplicates["CallRate_FID2"] = 0.0    # same as above

# for each row, we get the call rates of the FID1 and FID2 samples from the SamplesCallRates DF and compare them
# if somehow they end up being equal, then we remove the _R sample on the principle that it is the duplicate one
for index, row in PlannedDuplicates.iterrows():
    FID1 = row["FID1"]  # getting the IDs of the samples
    FID2 = row["FID2"]
    FID1_CallRate = float(SamplesCallRate[SamplesCallRate["FID"] == FID1]["Call_Rate"]) # gettin the call rates of the samples
    FID2_CallRate = float(SamplesCallRate[SamplesCallRate["FID"] == FID2]["Call_Rate"])
    PlannedDuplicates["CallRate_FID1"][index] = FID1_CallRate   # adding the call rates of the samples
    PlannedDuplicates["CallRate_FID2"][index] = FID2_CallRate
    # determining which sample to keep and which one to remove on the basis of their call rates
    if FID1_CallRate > FID2_CallRate:
        PlannedDuplicates["ToBeRemoved"][index] = FID2
        PlannedDuplicates["RemovalReason"][index] = f"{FID2} has a lower call rate than {FID1}"
    elif FID1_CallRate < FID2_CallRate:
        PlannedDuplicates["ToBeRemoved"][index] = FID1
        PlannedDuplicates["RemovalReason"][index] = f"{FID1} has a lower call rate than {FID2}"
    else:
        # the unlikely case that both the original sample and the planned duplicate sample has EXACTLY the same call rate, then we remove the one with _R
        PlannedDuplicateFID = ""
        if "_R" in str(FID1):
            PlannedDuplicateFID = FID1
        else:
            PlannedDuplicateFID = FID2
        PlannedDuplicates["ToBeRemoved"][index] = PlannedDuplicateFID
        PlannedDuplicates["RemovalReason"][index] = f"Since both samples have the same call rate, {PlannedDuplicateFID} is removed since it's a duplicate"
    
# for samples that are contaminated or cryptically related with each other, we not only compare their call rates, but also compare their age, gender and case control status
    # the logic goes like this :
        # Keep the case and remove the control.
        # If both are cases or controls, then keep the sample with the higher call rate (compare uptill the 4th decimal)
        # If even that ends up being the same, keep the older person
        # If even THAT is same, then do a comparison of the entire call rate value and keep the higher one

CryptRelSamples["ToBeRemoved"] = ""       # the sample IDs that need to be removed will be addded to this column
CryptRelSamples["RemovalReason"] = ""     # the reason for removal will be mentioned
CryptRelSamples["CallRate_FID1"] = 0.0    # this is the placeholder column for mentioning the call rate of the sample..
CryptRelSamples["CallRate_FID2"] = 0.0    # same as above
CryptRelSamples["Age_FID1"] = 0             # age of FID1
CryptRelSamples["Age_FID2"] = 0             # age of FID2
CryptRelSamples["Affected_Status_FID1"] = 0 # the phenotype of FID1
CryptRelSamples["Affected_Status_FID2"] = 0 # the phenotype of FID2

for index, row in CryptRelSamples.iterrows():
    FID1 = row["FID1"]  # getting the IDs of the samples
    FID2 = row["FID2"]
    FID1_CallRate = float(SamplesCallRate[SamplesCallRate["FID"] == FID1]["Call_Rate"]) # gettin the call rates of the samples
    FID2_CallRate = float(SamplesCallRate[SamplesCallRate["FID"] == FID2]["Call_Rate"])
    CryptRelSamples["CallRate_FID1"][index] = FID1_CallRate   # adding the call rates of the samples
    CryptRelSamples["CallRate_FID2"][index] = FID2_CallRate
    # to get the age of the samples
    FID1_Age, FID2_Age = int, int
    # FID1_Age = int(PedigreeFile[PedigreeFile["FID"] == FID1]["Age"])
    # FID2_Age = int(PedigreeFile[PedigreeFile["FID"] == FID2]["Age"])
    FID1_Age = int(PedigreeFile[PedigreeFile["FID"] == FID1]["Age"])
    FID2_Age = int(PedigreeFile[PedigreeFile["FID"] == FID2]["Age"])
    # now to get the case control status
    FID1_Pheno, FID2_Pheno = str, str
    CCTemp_FID1 = int(PedigreeFile[PedigreeFile["FID"] == FID1]["Affected_Status"])
    CCTemp_FID2 = int(PedigreeFile[PedigreeFile["FID"] == FID2]["Affected_Status"])
    if CCTemp_FID1 == 1:
        FID1_Pheno = "Control"
    elif CCTemp_FID1 == 2:
        FID1_Pheno = "Case"
    else:
        FID1_Pheno = "Missing_Pheno"
    # now to do the same for FID2
    if CCTemp_FID2 == 1:
        FID2_Pheno = "Control"
    elif CCTemp_FID2 == 2:
        FID2_Pheno = "Case"
    else:
        FID2_Pheno = "Missing_Pheno"
    # now to set the age and case control status of the samples
    CryptRelSamples["Age_FID1"][index] = FID1_Age
    CryptRelSamples["Age_FID2"][index] = FID2_Age
    CryptRelSamples["Affected_Status_FID1"][index] = FID1_Pheno
    CryptRelSamples["Affected_Status_FID2"][index] = FID2_Pheno
    # now to establish the logic chain as described before the for loop
    if FID1_Pheno == "Case":
        CryptRelSamples["ToBeRemoved"][index] = FID2
        CryptRelSamples["RemovalReason"][index] = f"{FID2} is a control whereas {FID1} is a case"
    elif FID2_Pheno == "Case":
        CryptRelSamples["ToBeRemoved"][index] = FID1
        CryptRelSamples["RemovalReason"][index] = f"{FID1} is a control whereas {FID2} is a case"
    elif (FID1_Pheno == "Case" and FID2_Pheno == "Case") or (FID1_Pheno == "Control" and FID2_Pheno == "Control") or (FID1_Pheno == "Missing_Pheno" and FID2_Pheno == "Missing_Pheno"):
        # both are cases, now we compare the call rates up to the 4th decimal
        FID1_CallRate_Rounded = round(FID1_CallRate,4)
        FID2_CallRate_Rounded = round(FID2_CallRate,4)
        if FID1_CallRate_Rounded > FID2_CallRate_Rounded: # one of the samples has a higher call rate when it's rounded up to the 4th decimal
            CryptRelSamples["ToBeRemoved"][index] = FID2
            CryptRelSamples["RemovalReason"][index] = f"{FID2} has a lower call rate than {FID1}"
        elif FID1_CallRate_Rounded < FID2_CallRate_Rounded:
            # keep FID2 and remove FID1
            CryptRelSamples["ToBeRemoved"][index] = FID1
            CryptRelSamples["RemovalReason"][index] = f"{FID1} has a lower call rate than {FID2}"
        else : # that is to say, that even the call rates are similar, then compare the ages
            if FID1_Age > FID2_Age:
                # keep FID1 and remove FID2
                CryptRelSamples["ToBeRemoved"][index] = FID2
                CryptRelSamples["RemovalReason"][index] = f"Sample {FID1} is older than {FID2}"
            elif FID1_Age < FID2_Age:
                # keep FID2 and remove FID1
                CryptRelSamples["ToBeRemoved"][index] = FID1
                CryptRelSamples["RemovalReason"][index] = f"Sample {FID2} is older than {FID1}"
            else:   # that is to say that if EVEN THE AGE IS SAME, then compare the full call rate number
                if FID1_CallRate > FID2_CallRate:
                    # keep FID1 and remove FID2
                    CryptRelSamples["ToBeRemoved"][index] = FID2
                    CryptRelSamples["RemovalReason"][index] = f"{FID2} has a lower call rate than {FID1}"
                elif FID1_CallRate < FID2_CallRate:
                    # keep FID2 and remove FID1
                    CryptRelSamples["ToBeRemoved"][index] = FID1
                    CryptRelSamples["RemovalReason"][index] = f"{FID1} has a lower call rate than {FID2}"
                else:
                    # just strictly remove the FID1...cuz these samples are pretty damn same anyway
                    CryptRelSamples["ToBeRemoved"][index] = FID1
                    CryptRelSamples["RemovalReason"][index] = f"{FID1} is selected for no particular reason. The call rates are same, age is same, phenotype is same."


##############################                                                                                      ##########################
################################            CASE WHERE 405_R gets contaminated with 577 or 577_R STILL REMAINS      ##########################
################################                                                                                    ##########################

###################################################

# now to create the removal list for cleaning the files of all cryptically related samples

CryptRel_RemovalList = list()
CryptRel_RemovalList.extend(PlannedDuplicates["ToBeRemoved"])
CryptRel_RemovalList.extend(CryptRelSamples["ToBeRemoved"])
CryptRel_RemovalList_Temp = pd.Series(CryptRel_RemovalList).unique()
CryptRel_RemovalList = pd.Series(CryptRel_RemovalList_Temp)

CryptRel_RemovalFile = pd.DataFrame({"FID" : CryptRel_RemovalList, "IID" : CryptRel_RemovalList})
CryptRel_RemovalFile_Len = len(CryptRel_RemovalFile)
CryptRel_RemovalFile.to_csv(f"{FilePath}/{FileName}_CryptRel_RemovalList_{CryptRel_RemovalFile_Len}-Samples.txt", sep="\t", index=False)
PlannedDuplicates_Len = len(PlannedDuplicates)
PlannedDuplicates.to_csv(f"{FilePath}/{FileName}_CryptRel_Planned_{PlannedDuplicates_Len}-PairwiseCombos.csv", index=False)
CryptRelSamples_Len = len(CryptRelSamples)
CryptRelSamples.to_csv(f"{FilePath}/{FileName}_CryptRel_TrueSamples_{CryptRelSamples_Len}-PairwiseCombos.csv", index=False)













cryptRelIID1, cryptRelIID2 = CryptRelDF["IID1"], CryptRelDF["IID2"]
PlannedDuplicates = list()
ContaminatedSamples = list()
ContaminatedDuplicates = list()
Contaminated = list()
FULLREMOVAL_TEMP = list()

# actual contaminations
ContaminationList_1, ContaminationList_2 = list(), list()
for IID1, IID2 in zip(cryptRelIID1, cryptRelIID2):
    # if neither IID1 not IID2 are planned duplis, BUT, they're still related to each other
    if ("_R" not in str(IID1)) and ("_R" not in str(IID2)):
        print(f"Found a contaminated sample.\t\tThis is IID1 : {str(IID1)} and this is IID2 : {str(IID2)}")
        #time.sleep(3)
        ContaminationList_1.append(IID1)
        ContaminationList_2.append(IID2)

if len(ContaminationList_1) <= len(ContaminationList_2):
    print(f"\nThe number of IID1 samples, that are contaminated : {len(ContaminationList_1)}")
    ContaminatedSamples = ContaminationList_1
else :
    print(f"\nThe number of IID2 samples, that are contaminated : {len(ContaminationList_2)}")
    ContaminatedSamples = ContaminationList_2

ContaminatedSamples = pd.Series(pd.Series(ContaminatedSamples).unique())
print("\n\n**********The Following are ContaminatedSamples : \n")
print(ContaminatedSamples)
print(len(ContaminatedSamples))



#
#    Example Sample ID   :   507 and 507_R
#Example Sample ID 2 :   455 and 455_R
#    Example Sample ID 3 :   801
    
#    Scenarios of contamination, ID1, ID2, PI_HAT  :
    
#        1) 507      455     0.99
#       2) 507      507_R   0.99
#      3) 507      455_R   0.70
#        4) 507_R    455_R   0.90
#        5) 507      801     0.50
        
#        Essentially,
#            A regular sample        with        another regular sample
#            A regular sample        with        it's planned duplicate
#            A regular sample        with        another sample's planned duplicate
#        and vice versa,
        
#        Or, to make it more generic :
#            A sample without _R     with        another sample without _R
#            A sample without _R     with        the same sample with _R
#            A sample without _R     with        another sample with _R
    
#    Cases :
#        1) 507      455     :   Compare call rates (Keep higher). If same, keep case. If same, keep balanced sex ratio. If same, keep older sample.
#        2) 507      507_R   :   Compare call rates (Keep higher).
#        3) 507      455_R   :   Do not auto-remove anything, but give the user the chance. Analyse the PI_HAT between the following 4 samples
#                                507     507_R       455     455_R   (Or any of these 4 that are available in the .genome file)
                                
#                                507     with    507_R   :   Use Case 2 above, & analyze the rest
#                                507     with    455     :   Refer Case 1 above
#                                507     with    455_R   :   Mention PI_HAT, & based on the value (Not near 100% or near 50%) suggest contamination ?
#                                507_R   with    455     :   Mention PI_HAT, & based on the value (Not near 100% or near 50%) suggest contamination ?
#                                507_R   with    455_R   :   Mention PI_HAT, & based on the value (Not near 100% or near 50%) suggest contamination ?
#                                455     with    455_R   :   Refer Case 2 above
                                
                                        
#    4) 507_R    455_R   0.90
#    5) 507      801     0.50
