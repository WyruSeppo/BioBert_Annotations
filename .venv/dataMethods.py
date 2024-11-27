from Bio import SeqIO
import requests
import pandas as pd
from classes import *

def read_tsv(file_path):
    try:
        df = pd.read_csv(file_path, sep='\t')
        return df
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return None


def get_uniprot_annotation(protein_id):
    url = f"https://www.uniprot.org/uniprotkb/{protein_id}.txt"
    print(url)
    try:
        response = requests.get(url)
    except Exception as e:
        print(f"An error occurred while calling uniprot API: {protein_id}")
        return None, None, None

    if response.status_code == 200:
        text_data = response.text
        protein_name = None
        function = None
        pfam = None
        endOfFunction = True
        
        for line in text_data.splitlines():
            if line.startswith("ID"):
                protein_name = line.split()[1]
            elif line.startswith("CC   -!- FUNCTION:"):
                endOfFunction = False
                function = line[len("CC   -!- FUNCTION:"):].strip()
            elif line.__contains__("-!-") and not endOfFunction:
                endOfFunction = True
            elif line.startswith("CC") and  not endOfFunction:
                function += line #remove whitespace and stuff at the end
            elif line.startswith("DR   Pfam;"):
                pfam = line.split()[2].removesuffix(';')

        return protein_name, function, pfam
    else:
        return None, None, None

def get_pfam_annotation(protein_id):
    api_url = "https://www.ebi.ac.uk/interpro/api"
    url = f"{api_url}/entry/pfam/{protein_id}"
    url += "?page_size=200"

    try:
        response = requests.get(url)
    except Exception as e:
        print(f"An error occurred while calling pfam API: {protein_id}")
        return None

    if response.status_code == 200:
        description_text = None
        jsonData = response.json()

        if jsonData['metadata']['description'] and isinstance(jsonData['metadata']['description'], list) and jsonData['metadata']['description'] != None and len(jsonData['metadata']['description']) > 0:
            description_text = jsonData['metadata']['description'][0]['text']
            description_text = description_text.removesuffix('</p>')
            description_text = description_text.removeprefix('<p>')
        return description_text
    else:
        return None


def annotate_tsv(tsv_file):
    df = read_tsv(tsv_file)
    annotated_data = []
    counter = 0
    total = 0

    if df is not None:
        total = len(df)
        for protein_id in df.iloc[:, 1]:
            print("\r" + str(counter / total) + "% " + str(counter) + "/" + str(total))
            
            protein_name, function, pfamID = get_uniprot_annotation(protein_id)
            pfamAnnotation = get_pfam_annotation(pfamID)

            annotated_data.append({
                'Protein_ID': protein_id,
                'Pfam_ID':pfamID,
                'Protein_Name': protein_name,
                'Unitprot_Function': function,
                'Pfam_Description': pfamAnnotation,
                'RefSeq':df.iloc[counter,0],
                'Entry':df.iloc[counter,1],
                'Entry_Name':df.iloc[counter,3],
                'Protein_Name':df.iloc[counter,4],
                'Gene_Names':df.iloc[counter,5],
                'Organism':df.iloc[counter,6],
                
            })
            counter += 1
            #print((protein_id or 'NA') + " " + (pfamID or 'NA') + " " + (protein_name or 'NA') + " " + (function or 'NA') + " " + (pfamAnnotation or 'NA'))
    return pd.DataFrame(annotated_data)

def loadAnnotations(filePaths):
    result = []

    # Open the file in read mode
    for filePath in filePaths:
        with open(filePath, 'r') as file:
            # Read each line in the file
            #Protein_ID	Pfam_ID	Protein_Name	Unitprot_Function	Pfam_Description	RefSeq	Entry	Entry_Name	Gene_Names	Organism
            for line in file:
                proteinId = line.split("\t")[0]
                pfamId = line.split("\t")[1]
                proteinName = line.split("\t")[2]
                uniprot_function = line.split("\t")[3]
                pfam_description = line.split("\t")[4]
                refSeqAccession = line.split("\t")[5]
                entry = line.split("\t")[6]
                entryName = line.split("\t")[7]
                geneNames = line.split("\t")[8]
                organism = line.split("\t")[9]
                result.append(AnnotationData(proteinId,pfamId,"","","", refSeqAccession, entry, entryName, proteinName, geneNames, organism, pfam_description,uniprot_function))
    
    return result[1:]#remove the first entry (heading)

def evaluateData(annotationData):
    data = EvaluatedData()
    data.no_sequences = len(annotationData)
    data.pfam_annotation_amount = sum(1 for x in annotationData if x.pfam_description != None)
    
     # Filter descriptions that are not empty
    descriptions = [x.pfam_description for x in annotationData if x.pfam_description != None]

    # Length of non-empty descriptions
    description_lengths = [len(desc) for desc in descriptions]
    
    if description_lengths:  # Check if there are valid descriptions
        data.pfam_annotation_length_avg = sum(description_lengths) / len(description_lengths)
        data.pfam_annotation_length_max = max(description_lengths)
        data.pfam_annotation_length_min = min(description_lengths)
    else:
        data.pfam_annotation_length_avg = 0
        data.pfam_annotation_length_max = 0
        data.pfam_annotation_length_min = 0
    
      # Count missing descriptions
    data.pfam_annotation_missing_amount = sum(1 for x in annotationData if x.pfam_description == "")
    
    # Calculate missing percentage
    if data.no_sequences > 0:
        data.pfam_annotation_missing_percent = (data.pfam_annotation_missing_amount / data.no_sequences) * 100
    else:
        data.pfam_annotation_missing_percent = 0

    # Calculate total number of words across all descriptions
    data.pfam_annotation_no_words = sum(len(desc.split()) for desc in descriptions)
    
    return data

def getAnnotations(inputFilePath, outputFilePath):
    annotated_df = annotate_tsv(inputFilePath)
    annotated_df.to_csv(outputFilePath, sep='\t', index=False)
