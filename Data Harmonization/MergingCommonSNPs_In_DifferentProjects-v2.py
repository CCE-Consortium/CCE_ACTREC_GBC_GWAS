import pandas as pd
import os

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

def read_bim(BIMFilePath : str) -> pd.DataFrame :
    df = pd.read_csv(BIMFilePath, sep="\s+", header=None, names = ["CHR", "SNP", "Morgans", "BP", "A1", "A2"])
    return df

# MainBIM_Path = str(input("Enter the path to the main BIM file, with the extension of .bim\t:\t"))
AllBIMs = str(input("Enter the paths to all the BIM files for common SNPs extraction, separated with a COMMA starting with the MAIN bim file\t:\t")
              ).split(",")
OutputPath = str(input("Enter the output path where the EXTRACT file & Common SNPs CSV file will be saved\t:\t"))
PlinkPath = str(input("Enter the full path to the PLINK 1.9 executable\t:\t"))
OutputName = str(input("Enter the name for the output BED fileset containing only the common variants from all the BIM files\t:\t"))

AllBIMLists = list()
AllBIMLists.extend(
    [element.strip() for element in AllBIMs]
    )

print("\n-------------------------------------------------------\n")
Bim_DF_list = list()
for BimPath in AllBIMLists:
    FileSetPath, FileSetName = split_path_name(FullPath=BimPath)
    Bim_DF_list.append(pd.read_csv(BimPath, sep="\s+", header=None, names = ["CHR", "SNP", "BP"], usecols=[0,1,3]))
    # comment the above and uncomment the below if you want to input the entire BIM file while reading
    # Bim_DF_list.append(pd.read_csv(BimPath, sep="\s+", header=None, names = ["CHR", "SNP", "Morgans", "BP", "A1", "A2"]))
    print(f"Processed\t:\t{FileSetName}\tBIM file")
    

SNPsCommonToAll = Bim_DF_list[0]
print("\n-------------------------------------------------------\n")
for index in range(1, len(Bim_DF_list)-1):
    Df_1 = Bim_DF_list[index]
    Df_2 = Bim_DF_list[index+1]
    # the below merge regardless of which file is considered to be the main. Only COMMON SNPs survive
    Merged = pd.merge(Df_1, Df_2, on=["CHR", "BP"], suffixes=("_x", "_y"))
    Merged["SNP"] = Merged["SNP_x"]
    Merged.drop(["SNP_x", "SNP_y"], axis=1, inplace=True)
    SNPsCommonToAll = pd.merge(SNPsCommonToAll, Merged,on=["CHR", "BP"], suffixes=("_x", "_y"))
    SNPsCommonToAll["SNP"] = SNPsCommonToAll["SNP_x"]
    SNPsCommonToAll.drop(["SNP_x", "SNP_y"], axis=1, inplace=True)
    # the below merge considering one file to be the main. The SNPs from those files all survive.
    # Merged = pd.merge(Df_1, Df_2, on=["CHR", "BP"], how="left")
    # SNPsCommonToAll = pd.merge(SNPsCommonToAll, Merged, how="left", on=["CHR", "BP"])
    print(f"Merged df at index 0 with those at index\t{index}\tand\t{index+1}")
    
SNPsCommonToAll.to_csv(f"{OutputPath}/Common-SNPs-BIM-file_{OutputName}.csv", index=False)
SNPsCommonToAll["SNP"].to_csv(f"{OutputPath}/Common-SNPs-EXTRACT-file_{OutputName}.extract", index=False, sep="\t", header=None)

print("\n-------------------------------------------------------\n")
MergeList = list()
for BimPath in AllBIMLists:
    FileSetPath, FileSetName = split_path_name(FullPath=BimPath)
    FileSetName = FileSetName.replace(".bim", "")
    cmd = f"{PlinkPath} --bfile {FileSetPath}/{FileSetName} --extract {OutputPath}/Common-SNPs-EXTRACT-file_{OutputName}.extract --make-bed --out {OutputPath}/CommonSNPsExtracted_{FileSetName}"
    os.system(cmd)
    MergeList.append(f"{OutputPath}/CommonSNPsExtracted_{FileSetName}")

MergeList_DF = pd.DataFrame({"FilesToMerge" : MergeList[1:]})
MergeList_DF.to_csv(f"{OutputPath}/MergeList_ExcludingFirstBIMFile.txt", sep="\t", index=False, header=None)
MergeCMD = f"{PlinkPath} --bfile {MergeList[0]} --allow-extra-chr --merge-list {OutputPath}/MergeList_ExcludingFirstBIMFile.txt --make-bed --out {OutputPath}/All_BED-FileSets_MERGED_Pass1"
os.system(MergeCMD)
