from Bio import SeqIO
import time
import requests
import pandas as pd
from classes import *
import logging
from apiMethods import *
import configparser

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("BioBERT")

def read_tsv(file_path):
    try:
        df = pd.read_csv(file_path, sep='\t')
        return df
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return None


def loadAnnotations(filePath):
    result = []

    with open(filePath, 'r') as file:
        
        skipFirstLine = True
            
        for line in file:
            
            if skipFirstLine:
                skipFirstLine = False
                continue
            
            # Split the line by tabs
            columns = line.strip().split('\t')
            
            # Extract each field based on its position
            proteinId = columns[0]
            pfamId = columns[1]
            uniprotId = columns[2]
            pfam_embedding = columns[3]
            uniprot_embedding = columns[4]
            refSeqAccession = columns[5]
            entry = columns[6]
            entryName = columns[7]
            proteinName = columns[8]
            geneNames = columns[9]
            organism = columns[10]
            pfam_description = columns[11]
            uniprot_function = columns[12] if len(columns) > 12 else ""

            # Create an AnnotationData object and append it to the result list
            result.append(
                AnnotationData(
                    id=proteinId,
                    pfam_id=pfamId,
                    uniprot_id=uniprotId,
                    pfam_embedding=pfam_embedding,
                    uniprot_embedding=uniprot_embedding,
                    refSeqAccession=refSeqAccession,
                    entry=entry,
                    entryName=entryName,
                    proteinNames=proteinName,
                    geneNames=geneNames,
                    organism=organism,
                    pfam_description=pfam_description,
                    uniprot_function=uniprot_function
                )
            )

    return result

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

def getAnnotations(annotationData, outputFilePath):
    annotationData = annotate_data(annotationData)
    exportData = pd.DataFrame([obj.to_dict() for obj in annotationData])
    exportData.to_csv(outputFilePath, sep='\t', index=False)
    return annotationData
    
def saveAnnotations(annotationData, outputFilePath):
    exportData = pd.DataFrame([obj.to_dict() for obj in annotationData])
    exportData.to_csv(outputFilePath, sep='\t', index=False)
    

def load_ids_fasta(fasta_file):
    logger.info(f"Loading FASTA file: {fasta_file}")
    ref_seq_ids = [record.id for record in SeqIO.parse(fasta_file, "fasta")]
    return ref_seq_ids

def configIsValid(config):
    #check config for things
    return True

def writeToFile(input, outputfile):
    with open(outputfile, "w") as f:
        f.write(input)
        
def read_config(filePath = 'biobert.ini'):
    #https://www.geeksforgeeks.org/how-to-write-a-configuration-file-in-python/
    
    # Create a ConfigParser object
    config = configparser.ConfigParser()

    # Read the configuration file
    config.read(filePath)

    # Access values from the configuration file
    fasta_file = config.get('General','fasta_file')
    annotation_file_input = config.get('General','annotation_file_input')
    annotation_file_output = config.get('General','annotation_file_output')
    annotation_embedding_file_output = config.get('General','annotation_embedding_file_output')
    data_eval_output = config.get('General','data_eval_output')
    loadAnnotationsFromFile = config.getboolean('General','loadAnnotationsFromFile')
    getPfamEmbeddings = config.getboolean('General','getPfamEmbeddings')
    getUniProtEmbeddings = config.getboolean('General','getUniProtEmbeddings')
    model = config.get('General','model')


    # Return a dictionary with the retrieved values
    config_values = {
        'fasta_file' : fasta_file,
        'annotation_file_input' : annotation_file_input,
        'annotation_file_output' : annotation_file_output,
        'annotation_embedding_file_output' : annotation_embedding_file_output,
        'data_eval_output' : data_eval_output,
        'loadAnnotationsFromFile' : loadAnnotationsFromFile,
        'getPfamEmbeddings' : getPfamEmbeddings,
        'getUniProtEmbeddings' : getUniProtEmbeddings,
        'model' : model
    }

    return config_values
    
    