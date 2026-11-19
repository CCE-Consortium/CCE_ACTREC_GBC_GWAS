import os
import pandas as pd
import time
import threading as thr
import random

def split_path_name(FullPath : str) -> list:
    """Given a full path to a file, returns the path to the folder, and the file name

    Args:
        FullPath (str): The full path to a file.

    Returns:
        list: A list containing the full path to the folder containing the file as the 0th element, and the file name as the 1st element.
    """
    PathElements = FullPath.split("/")
    FilePath = ""
    # reconstructing path to the folder
    for p in PathElements[1:len(PathElements)-1]:
        FilePath += f"/{p}"
    # getting the filename
    FileName = PathElements[len(PathElements)-1]
    return [FilePath, FileName]

def threads_list_run(ThreadsList : list) -> None:
    """Runs all threads in a given list, waits till completion, then returns

    Args:
        * ThreadsList (list): A list of all the threads that need to be run till completion
    """
    for t in ThreadsList:
        t.start()   
    for t in ThreadsList:
        t.join()
    return

def execute_command(CommandString : str) -> None:
    os.system(CommandString)

def run_logistic_etc(InputFilePath : str, OutputFilePath : str, Adjustments : list, AdjustmentString : str, PlinkPath : str, ConfidenceInterval = 0.95) -> dict:
    ReturnDict = dict()
    FilePath, FileName = split_path_name(FullPath=InputFilePath)
    # setting the seed for reproducability
    Seed = random.randint(0, 4294967000)
    # setting the permutation batch size & threads number for reproducability
    PermBatchSize = 512
    Threads = 10
    # if there are covariates that need to be adjusted for, then the names of those covariates will be in Adjustments list. Take each term from there, and add it to the output file names
    OutputName_Common = ""
    if len(Adjustments) > 0:
        for term in Adjustments:
            OutputName_Common += f"_{term}"
    else:
        OutputName_Common = "_"
    OutputName_OR = f"{FileName}_Logistic{OutputName_Common}_OR"
    LogisticCommand_OR = f"{PlinkPath} --bfile {InputFilePath} --logistic hide-covar --threads {Threads} --perm-batch-size {PermBatchSize} --seed {Seed} --out {OutputFilePath}/{OutputName_OR} {AdjustmentString}"
    OutputName_SE = f"{FileName}_Logistic{OutputName_Common}_SE"
    LogisticCommand_SE = f"{PlinkPath} --bfile {InputFilePath} --logistic beta hide-covar --ci {ConfidenceInterval} --threads {Threads} --perm-batch-size {PermBatchSize} --seed {Seed} --out {OutputFilePath}/{OutputName_SE} {AdjustmentString}"
    freq_filters_gender, freq_filters_phenos = "", ""
    if "--filter-males" in AdjustmentString:
        freq_filters_gender = "--filter-males"
    elif "--filter-females" in AdjustmentString:
        freq_filters_gender = "--filter-females"
    if "--filter-cases" in AdjustmentString:
        freq_filters_phenos = "--filter-cases"
    elif "--filter-controls" in AdjustmentString:
        freq_filters_phenos = "--filter-controls"
    freq_filter_commands = f"{freq_filters_gender} {freq_filters_phenos}"
    Frequency_Command =                 f"{PlinkPath} --bfile {InputFilePath} {freq_filter_commands} --freq --out {OutputFilePath}/{FileName}_Freq"
    Frequency_Command_CaseControl =     f"{PlinkPath} --bfile {InputFilePath} {freq_filter_commands} --freq case-control --out {OutputFilePath}/{FileName}_Freq_CaseControl"
    Frequency_X_Command =               f"{PlinkPath} --bfile {InputFilePath} {freq_filter_commands} --freqx --out {OutputFilePath}/{FileName}_Freq_X"
    Frequency_Counts_Command =          f"{PlinkPath} --bfile {InputFilePath} {freq_filter_commands} --freq counts --out {OutputFilePath}/{FileName}_Freq_Counts"
    CommandsList = [LogisticCommand_OR, LogisticCommand_SE, Frequency_Command, Frequency_Command_CaseControl, Frequency_X_Command, Frequency_Counts_Command]
    CommandsName = ["LogisticCommand_OR", "LogisticCommand_SE", "Frequency_Command", "Frequency_Command_CaseControl", "Frequency_X_Command", "Frequency_Counts_Command"]
    ThreadsList = list()
    for cmd, cmdName in zip(CommandsList, CommandsName):
        ThreadObject = thr.Thread(target=execute_command, name=cmdName, kwargs={"CommandString" : cmd})
        ThreadsList.append(ThreadObject)
    threads_list_run(ThreadsList=ThreadsList)
    # adding values to the return dictionary
    ReturnDict["Seed"] = Seed
    ReturnDict["PermBatchSize"] = PermBatchSize
    ReturnDict["Threads"] = Threads
    ReturnDict["Confidence_Interval"] = ConfidenceInterval
    ReturnDict["Logistic_OR"] = f"{OutputFilePath}/{OutputName_OR}.assoc.logistic"
    ReturnDict["Logistic_SE"] = f"{OutputFilePath}/{OutputName_SE}.assoc.logistic"
    ReturnDict["Freq"] = f"{OutputFilePath}/{FileName}_Freq.frq"
    ReturnDict["Freq_CaseControl"] = f"{OutputFilePath}/{FileName}_Freq_CaseControl.frq.cc"
    ReturnDict["Freq_X"] = f"{OutputFilePath}/{FileName}_Freq_X.frqx"
    ReturnDict["Freq_Counts"] = f"{OutputFilePath}/{FileName}_Freq_Counts.frq.counts"
    return ReturnDict

def merge_results(OutputFilePaths_Dict : dict) -> str :
    Logistic_OR = OutputFilePaths_Dict["Logistic_OR"]
    Logistic_SE = OutputFilePaths_Dict["Logistic_SE"]
    Freq = OutputFilePaths_Dict["Freq"]
    Freq_CaseControl = OutputFilePaths_Dict["Freq_CaseControl"]
    Freq_X = OutputFilePaths_Dict["Freq_X"]
    Freq_Counts = OutputFilePaths_Dict["Freq_Counts"]
    # fetching the main result file now, the one with the Odds Ratio column and without the additional nonsense
    Logistic_OR_DF = pd.read_csv(Logistic_OR, sep="\s+")
    print("File with Odds ratio has been processed")
    #fetching the frequency report files now
    # and one by one, merging them after renaming the columns
    Freq_DF = pd.read_csv(Freq, sep="\s+",usecols=["SNP","A1","A2","MAF"])
    Freq_DF.columns = ["SNP","A1_MinorAllele","A2_MajorAllele","A1_MinorAlleleFreq"]
    print("File with minor allele frequencies has been processed")
    Logistic_OR_DF = pd.merge(Logistic_OR_DF, Freq_DF, on="SNP")
    del(Freq_DF)
    Freq_CaseControl_DF = pd.read_csv(Freq_CaseControl, sep="\s+", usecols=["SNP", "MAF_A", "MAF_U", "NCHROBS_A", "NCHROBS_U"])
    Freq_CaseControl_DF.columns=["SNP", "MAF_In_Cases", "MAF_In_Controls", "NumOf_Allele_Obs_In_Cases", "NumOf_Allele_Obs_In_Controls"]
    print("File with allele frequencies in cases and controls has been processed")
    Logistic_OR_DF = pd.merge(Logistic_OR_DF, Freq_CaseControl_DF, on="SNP")
    del(Freq_CaseControl_DF)
    Freq_X_DF = pd.read_csv(Freq_X, sep="\t",usecols=["SNP","C(HOM A1)","C(HET)","C(HOM A2)","C(HAP A1)","C(HAP A2)"])
    Freq_X_DF.columns = ["SNP","A1_HomozygCount","HeteroCount", "A2_HomozygCount", "A1_HapCount","A2_HapCount"]
    print("File with X chr allele frequencies has been processed")
    Logistic_OR_DF = pd.merge(Logistic_OR_DF, Freq_X_DF, on="SNP")
    del(Freq_X_DF)
    Freq_Counts_DF = pd.read_csv(Freq_Counts, sep="\s+",usecols=["SNP","C1","C2"])
    Freq_Counts_DF.columns = ["SNP","A1_Count", "A2_Count"]
    print("File with allele counts has been processed")
    Logistic_OR_DF = pd.merge(Logistic_OR_DF, Freq_Counts_DF, on="SNP")
    del(Freq_Counts_DF)
    # fetching the association result files with SE and conf interval etc
    ConfidenceInterval_Header = str(float(OutputFilePaths_Dict["Confidence_Interval"])*100)
    ConfidenceInterval_Header = ConfidenceInterval_Header.replace(".0", "")
    L_Name = f"L{ConfidenceInterval_Header}"
    U_Name = f"U{ConfidenceInterval_Header}"
    Logistic_SE_DF = pd.read_csv(Logistic_SE, sep="\s+", usecols=["SNP", "BETA", "SE", L_Name, U_Name])
    print("File with SE, BETA and confidence intervals has been processed")
    Logistic_OR_DF = pd.merge(Logistic_OR_DF, Logistic_SE_DF, on="SNP")
    print("Final merge completed")
    del(Logistic_SE_DF)
    FilePath, FileName = split_path_name(FullPath=Logistic_OR)
    Logistic_OR_DF.to_csv(f"{FilePath}/{FileName}", sep="\t", index=False)
    print(Logistic_OR_DF.info())
    print(Logistic_OR_DF.head(5))


PlinkPath = str(input("Enter the path to plink executable\t:\t"))
InputBEDFilesPath = str(input("Enter the path to the input BED files for analysis without the extension\t:\t"))
OutputPath = str(input("Enter the full path to the folder where output will go\t:\t"))

Covar_Adjustment_Names = str(input("Enter the covariates - SPACE SEPARATED - (age, gen, PCs etc), gender or case-control filters that are applied to this logistic regression\t:\t"))
Adjustments_List = Covar_Adjustment_Names.split(" ")
CovarCmdLine = str(input("Enter the entire covariate command as a single line\n\t1) --covar <path to covariate file> keep-pheno-on-missing-cov\
    \n\t2) --covar number <the header numbers of the covariate to adjust to>\
        \n\t3) --filter-males OR --filter-females OR --filter-cases OR --filter-controls if any\n\n\t\tEnter the line\t:\t"))

ParamsDict = run_logistic_etc(InputFilePath=InputBEDFilesPath, OutputFilePath=OutputPath, Adjustments=Adjustments_List, AdjustmentString=CovarCmdLine, PlinkPath=PlinkPath)
merge_results(OutputFilePaths_Dict=ParamsDict)
