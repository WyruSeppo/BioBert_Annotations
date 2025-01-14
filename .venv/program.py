from dataMethods import *
from apiMethods import *
from bertMethods import *
import logging
import numpy as np

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("BioBERT")

CONFIG = read_config('C:\\Users\\sebastian.rossboeck\\OneDrive - DO & CO Aktiengesellschaft\\Desktop\\biobert\\BioBert_Annotations\\.venv\\biobert.ini')

logger.info("BioBert Annotation Similarity v1 Start") 

#1 evaluate the config to catch invalid inputs early
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
for entry in annotationData:
    entry.cleanAnnotations()

#3 Data Eval on annotations
logger.info("Evaluate Data")
writeToFile(evaluateData(annotationData).generateString(), CONFIG["data_eval_output"])

#3.5 Clean Data: Remove entries with annotations that are longer than 512 words
#remove entries that don't have BOTH a pfam-description AND uniprot-function.
print("before prune1:" + str(len(annotationData)))
annotationData = [
    item for item in annotationData 
    if isinstance(getattr(item, 'pfam_description', None), str) 
    and isinstance(getattr(item, 'uniprot_function', None), str)
    and item.pfam_description.strip() != "None"  # Ensure it's not the string "None"
    and item.uniprot_function.strip() != "None"  # Ensure it's not the string "None"
]

print("before prune1:" + str(len(annotationData)))
annotationData = [item for item in annotationData if len(item.pfam_description.split()) < 512]
print("before prune3:" + str(len(annotationData)))
annotationData = [item for item in annotationData if len(item.uniprot_function.split()) < 512]
print("after prune3:" + str(len(annotationData)))

saveAnnotations(annotationData, CONFIG["annotation_embedding_file_output"])

logger.info("Evaluate Data again")
writeToFile(evaluateData(annotationData).generateString(), CONFIG["data_eval_output2"])

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


#calculate the distance between the embeddings
for entry in annotationData:
    # Check if pfam_embedding and uniprot_embedding are valid (i.e., not 'None' or empty)
    if entry.pfam_embedding != "None" and entry.pfam_embedding:
        pfam_embedding_array = np.array([float(val) for val in entry.pfam_embedding.split()])
    else:
        pfam_embedding_array = np.array([])  # or handle it differently if needed

    if entry.uniprot_embedding != "None" and entry.uniprot_embedding:
        uniprot_embedding_array = np.array([float(val) for val in entry.uniprot_embedding.split()])
    else:
        uniprot_embedding_array = np.array([])  # or handle it differently if needed

    # Now you can check if the embeddings are valid before calculating similarity
    if pfam_embedding_array.size > 0 and uniprot_embedding_array.size > 0:
        cosine_similarity = np.dot(pfam_embedding_array, uniprot_embedding_array) / (
            np.linalg.norm(pfam_embedding_array) * np.linalg.norm(uniprot_embedding_array)
        )
        entry.embedding_distance = 1 - cosine_similarity
        print("----------")
        print(entry.uniprot_function)
        print(entry.pfam_description)
        print("Cosine Similarity:", cosine_similarity)
        print("----------")
    else:
        print("Invalid embeddings detected")
        entry.embedding_distance = 2

saveAnnotations(annotationData, CONFIG["annotation_embedding_file_output"])

#evaluate Embedding-Distances
distances = [item.embedding_distance for item in annotationData]

min_distance = min(distances)
max_distance = max(distances)
avg_distance = sum(distances) / len(distances) if distances else 0 

print("Min Distance:", min_distance)
print("Max Distance:", max_distance)
print("Avg Distance:", avg_distance)




#visualize encoding data
print("BioBert Annotation Similarity v1 End") 
