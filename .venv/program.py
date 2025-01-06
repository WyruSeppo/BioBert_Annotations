from dataMethods import *
from apiMethods import *
from bertMethods import *
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("BioBERT")

CONFIG = read_config('C:\\Users\\SebastianRossboeck\\Desktop\\BioBert4\\.venv\\biobert.ini')

logger.info("BioBert Annotation Similarity v1 Start") 

#1 evaluate the config to catch invalid inputs early
if configIsValid(CONFIG) == False:
    print("CONFIG is not valid")
    quit()

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
    
    #convert to uniprot and initialize list of AnnotationData objects
    try:
        annotationData = getUniProtConversion("RefSeq_Protein","UniProtKB",ref_seq_ids)
    except Exception as e:
        logger.error(f"Error loading conversion from uniprot api: {e}")
        
    #get Annotations
    #this automatically creates a file with the annotationdata we receive from the apis
    try:
        annotationData = getAnnotations(annotationData, CONFIG["annotation_file_output"])
    except Exception as e:
        logger.error(f"Error getting annotations from api: {e}")

#2.5 Clean annotations
annotationData.cleanAnnotations()

#3 Data Eval on annotations
logger.info("Evaluate Data")
writeToFile(evaluateData(annotationData).generateString(), CONFIG["data_eval_output"])


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
