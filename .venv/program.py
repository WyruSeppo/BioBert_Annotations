from dataMethods import *
from bertMethods import *
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("BioBERT")

CONFIG = {
    "fasta_file": "C:\\Users\\SebastianRossboeck\\Desktop\\BioBert4\\spore-formers_test.faa",
    "annotation_file_input": "C:\\Users\\SebastianRossboeck\\Desktop\\BioBert4\\annotationData_output_test.txt",
    "annotation_file_output": "C:\\Users\\SebastianRossboeck\\Desktop\\BioBert4\\annotationData_output_test.txt",
    "annotation_embedding_file_output": "C:\\Users\\SebastianRossboeck\\Desktop\\BioBert4\\annotationData_embedding_output_test.txt",
    "data_eval_output": "C:\\Users\\SebastianRossboeck\\Desktop\\BioBert4\\data_eval_output_test.txt",
    "loadAnnotationsFromFile":True,
    "getPfamEmbeddings":False,
    "getUniProtEmbeddings":False,
    "model":"dmis-lab/biobert-base-cased-v1.1"
}

logger.info("BioBert Annotation Similarity v1 Start") 

#1 evaluate the config to catch invlid inputs early
if configIsValid(CONFIG) == False:
    print("CONFIG is not valid")
    quit()


#2 if there exists a file with annotationData: load it, otherwise fetch annotations from Uniprot/Pfam
if CONFIG["loadAnnotationsFromFile"]:
    #load annotationData from file
    logger.info("loading annotationdata from file")
    try:
        annotationData = loadAnnotations(CONFIG["annotation_file_input"])
    except Exception as e:
        logger.error(f"Error fetching annotations from file: {e}")
        quit()
else :
    logger.info("fetching annotationdata from uniprot/pfam api")
    
    #read refSeqIds from fasta file
    try:
        ref_seq_ids = load_ids_fasta(CONFIG["fasta_file"])
    except Exception as e:
        logger.error(f"Error loading from fasta file: {e}")
    
    try:
        #convert to uniprot and initialize list of AnnotationData objects
        annotationData = getUniProtConversion("RefSeq_Protein","UniProtKB",ref_seq_ids)
    except Exception as e:
        logger.error(f"Error loading conversion from uniprot api: {e}")
        
    #get Annotations
    #this creates a file with the annotationdatawe receive form the apis
    try:
        annotationData = getAnnotations(annotationData, CONFIG["annotation_file_output"])
    except Exception as e:
        logger.error(f"Error getting annotations from api: {e}")


#3 Data Eval on annotations
logger.info("Evaluate Data")
evaluatedData = evaluateData(annotationData)
#evaluatedData.print_data()
writeToFile(evaluatedData.generateString(), CONFIG["data_eval_output"])


#4.1 get Encodings for pfam-description
logger.info("Create Embeddings")
if CONFIG["getPfamEmbeddings"]:
    annotationData = getEmbeddings(annotationData, CONFIG["model"],"pfam")
    

#4.2 get Encodings for uniprot-function
if CONFIG["getUniProtEmbeddings"]:
    annotationData = getEmbeddings(annotationData, CONFIG["model"],"uniprot")
 
#If we created new embeddings: save the data to file
if CONFIG["getUniProtEmbeddings"] or CONFIG["getPfamEmbeddings"]:
    saveAnnotations(annotationData, CONFIG["annotation_embedding_file_output"])
 



#cluster the encodings data

#visualize encoding data
print("BioBert Annotation Similarity v1 End") 
