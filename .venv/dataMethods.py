import os
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
    """
    Reads a TSV file and returns its contents as a Pandas DataFrame.

    Parameters:
        file_path (str): Path to the TSV file.

    Returns:
        pd.DataFrame or None: The DataFrame containing file data, or None if an error occurs.
    """
    try:
        df = pd.read_csv(file_path, sep='\t')
        return df
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return None


def loadAnnotations(filePath):
    """
    Loads annotation data from a TSV file and converts it into a list of AnnotationData objects.

    Parameters:
        filePath (str): Path to the TSV file containing annotation data.

    Returns:
        list[AnnotationData]: A list of AnnotationData objects parsed from the file.
    """
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
            embedding_distance = columns[13] if len(columns) > 13 else ""

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
                    uniprot_function=uniprot_function,
                    embedding_distance=embedding_distance
                )
            )

    return result

def evaluateData(annotationData):
    """
    Evaluates annotation data to compute statistics on Pfam and UniProt annotations.

    Parameters:
        annotationData (list[AnnotationData]): A list of annotation data objects.

    Returns:
        EvaluatedData: An object containing statistical analysis of the annotation data.
    """
    data = EvaluatedData()
    data.no_sequences = len(annotationData)
    data.pfam_annotation_amount = sum(1 for x in annotationData if x.pfam_description != None and x.pfam_description != "None")
    
     # Filter descriptions that are not empty
    descriptions = [x.pfam_description for x in annotationData if x.pfam_description != None and x.pfam_description != "None"]

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
    data.pfam_annotation_missing_amount = sum(1 for x in annotationData if x.pfam_description == "None")
    
    # Calculate missing percentage
    if data.no_sequences > 0:
        data.pfam_annotation_missing_percent = (data.pfam_annotation_missing_amount / data.no_sequences) * 100
    else:
        data.pfam_annotation_missing_percent = 0

    # Calculate total number of words across all descriptions
    data.pfam_annotation_no_words = sum(len(desc.split()) for desc in descriptions)
    
    #again for uniprot
    data.uniprot_annotation_amount = sum(1 for x in annotationData if x.uniprot_function != None and x.uniprot_function != "None")
    
    #Filter descriptions that are not empty
    descriptions = [x.uniprot_function for x in annotationData if x.uniprot_function != None and x.uniprot_function != "None"]

    # Length of descriptions
    description_lengths = [len(desc) for desc in descriptions]
    
    if description_lengths: 
        data.uniprot_annotation_length_avg = sum(description_lengths) / len(description_lengths)
        data.uniprot_annotation_length_max = max(description_lengths)
        data.uniprot_annotation_length_min = min(description_lengths)
    else:
        data.uniprot_annotation_length_avg = 0
        data.uniprot_annotation_length_max = 0
        data.uniprot_annotation_length_min = 0
    
    # Count missing descriptions
    data.uniprot_annotation_missing_amount = sum(1 for x in annotationData if x.uniprot_function == "None")
    
    # Calculate missing percentage
    if data.no_sequences > 0:
        data.uniprot_annotation_missing_percent = (data.uniprot_annotation_missing_amount / data.no_sequences) * 100
    else:
        data.uniprot_annotation_missing_percent = 0

    # Calculate total number of words across all descriptions
    data.uniprot_annotation_no_words = sum(len(desc.split()) for desc in descriptions)
    
    
    return data

def getAnnotations(annotationData, outputFilePath):
    """
    Annotates data using an external method and saves the results to a TSV file.

    Parameters:
        annotationData (list[AnnotationData]): A list of annotation data objects.
        outputFilePath (str): Path to save the annotated data.

    Returns:
        list[AnnotationData]: The updated annotation data list.
    """
    annotationData = annotate_data(annotationData)
    exportData = pd.DataFrame([obj.to_dict() for obj in annotationData])
    exportData.to_csv(outputFilePath, sep='\t', index=False)
    return annotationData
    
def saveAnnotations(annotationData, outputFilePath):
    """
    Saves annotation data to a TSV file.

    Parameters:
        annotationData (list[AnnotationData]): A list of annotation data objects.
        outputFilePath (str): Path to save the data.
    """
    exportData = pd.DataFrame([obj.to_dict() for obj in annotationData])
    exportData.to_csv(outputFilePath, sep='\t', index=False)
    

def load_ids_fasta(fasta_file):
    """
    Loads sequence IDs from a FASTA file.

    Parameters:
        fasta_file (str): Path to the FASTA file.

    Returns:
        list[str]: A list of sequence IDs from the FASTA file.
    """
    logger.info(f"Loading FASTA file: {fasta_file}")
    ref_seq_ids = [record.id for record in SeqIO.parse(fasta_file, "fasta")]
    return ref_seq_ids

def can_save_file(file_path):
    """
    Checks whether a file can be saved in the specified directory.

    Parameters:
        file_path (str): Path to the file.

    Returns:
        bool: True if the file can be saved, False otherwise.
    """
    try:
        # Get the directory part of the path
        directory = os.path.dirname(file_path) or "."
        
        # Check if the directory exists
        if not os.path.exists(directory):
            return False
        
        # Check if the directory is writable
        return os.access(directory, os.W_OK)
    except Exception as e:
        return False

def configIsValid(config):
    """
    Validates the configuration by checking if files can be written to specified paths.

    Parameters:
        config (dict): Dictionary containing file paths from the configuration.

    Returns:
        bool: True if all paths are writable, False otherwise.
    """
    #check if we can write in each of the directories supplied in the config
    valid = True
    valid = valid and can_save_file(config["fasta_file"])
    valid = valid and can_save_file(config["annotation_file_input"])
    valid = valid and can_save_file(config["annotation_file_output"])
    valid = valid and can_save_file(config["annotation_embedding_file_output"])
    valid = valid and can_save_file(config["data_eval_output"])

    return valid

def writeToFile(input, outputfile):
    """
    Writes a string to a file.

    Parameters:
        input (str): The content to be written.
        outputfile (str): Path to the output file.
    """
    with open(outputfile, "w") as f:
        f.write(input)
        
def read_config(filePath='biobert.ini'):
    """
    Reads a configuration file and extracts relevant settings.

    Parameters:
        filePath (str, optional): Path to the configuration file. Defaults to 'biobert.ini'.

    Returns:
        dict: A dictionary containing configuration values.
    """
    #https://www.geeksforgeeks.org/how-to-write-a-configuration-file-in-python/
    
    # Create a ConfigParser object
    config = configparser.ConfigParser()

    # Read the configuration file
    config.read(filePath)

    # Get values from the configuration file
    fasta_file = config.get('General','fasta_file')
    annotation_file_input = config.get('General','annotation_file_input')
    annotation_file_output = config.get('General','annotation_file_output')
    annotation_embedding_file_output = config.get('General','annotation_embedding_file_output')
    data_eval_output = config.get('General','data_eval_output')
    data_eval_output2 = config.get('General','data_eval_output2')
    loadAnnotationsFromFile = config.getboolean('General','loadAnnotationsFromFile')
    getPfamEmbeddings = config.getboolean('General','getPfamEmbeddings')
    getUniProtEmbeddings = config.getboolean('General','getUniProtEmbeddings')
    model = config.get('General','model')
    distances_plot_output = config.get('General','distances_plot_output')
    tsne_plot_output = config.get('General','tsne_plot_output')


    # Return a dictionary with the retrieved values
    config_values = {
        'fasta_file' : fasta_file,
        'annotation_file_input' : annotation_file_input,
        'annotation_file_output' : annotation_file_output,
        'annotation_embedding_file_output' : annotation_embedding_file_output,
        'data_eval_output' : data_eval_output,
        'data_eval_output2' : data_eval_output2,
        'loadAnnotationsFromFile' : loadAnnotationsFromFile,
        'getPfamEmbeddings' : getPfamEmbeddings,
        'getUniProtEmbeddings' : getUniProtEmbeddings,
        'model' : model,
        'distances_plot_output' : distances_plot_output,
        'tsne_plot_output' : tsne_plot_output
    }

    return config_values
    
    