from dataMethods import *
from apiMethods import *
from bertMethods import *
import logging
import numpy as np
import matplotlib.pyplot as plt

def main():
    """Executes the BioBERT annotation similarity workflow.

    This includes:
        - Reading and validating the configuration.
        - Loading or fetching annotation data.
        - Cleaning and evaluating the data.
        - Generating embeddings.
        - Calculating cosine similarity between embeddings.
        - Saving results to output files.
    """
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("BioBERT")

    CONFIG = read_config('C:\\Users\\SebastianRossboeck\\Desktop\\BioBert4\\.venv\\biobert.ini')

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
    #also, remove duplicates by saving the tuple of annotations in a Set
    #also remove annotations longer than 512 words, those would be truncated by biobert
    seen = set()
    annotationData = [
        item for item in annotationData 
        if isinstance(getattr(item, 'pfam_description', None), str) 
        and isinstance(getattr(item, 'uniprot_function', None), str)
        and item.pfam_description.strip() != "None"
        and item.uniprot_function.strip() != "None"
        and (key := (item.pfam_description.strip(), item.uniprot_function.strip())) not in seen
        and not seen.add(key)
        and len(item.pfam_description.split()) < 512
        and len(item.uniprot_function.split()) < 512
    ]

    saveAnnotations(annotationData, CONFIG["annotation_embedding_file_output"])

    logger.info("Evaluate reduced dataset")
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


    #5 calculate the distance between the embeddings
    for entry in annotationData:
        # Check if pfam_embedding and uniprot_embedding are valid 
        if entry.pfam_embedding != "None" and entry.pfam_embedding:
            pfam_embedding_array = np.array([float(val) for val in entry.pfam_embedding.split()])
        else:
            pfam_embedding_array = np.array([])

        if entry.uniprot_embedding != "None" and entry.uniprot_embedding:
            uniprot_embedding_array = np.array([float(val) for val in entry.uniprot_embedding.split()])
        else:
            uniprot_embedding_array = np.array([])

        if pfam_embedding_array.size > 0 and uniprot_embedding_array.size > 0:
            cosine_similarity = np.dot(pfam_embedding_array, uniprot_embedding_array) / (
                np.linalg.norm(pfam_embedding_array) * np.linalg.norm(uniprot_embedding_array)
            )
            entry.embedding_distance = 1 - cosine_similarity
        else:
            print("Invalid embeddings detected")
            entry.embedding_distance = 2 #the maximum

    saveAnnotations(annotationData, CONFIG["annotation_embedding_file_output"])

    #5.1 evaluate Embedding-Distances
    distances = [item.embedding_distance for item in annotationData]
    min_distance = min(distances)
    max_distance = max(distances)
    avg_distance = sum(distances) / len(distances) if distances else 0 

    #show histogram
    plt.figure(figsize=(8, 5))
    plt.hist(distances, bins=50, edgecolor="black", alpha=0.75)
    plt.xlabel("Cosine Distance")
    plt.ylabel("Frequency")
    plt.title("Distribution of Cosine Distances")
    plt.xlim(0, 2) 
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    
    plt.text(2, plt.ylim()[1] * 0.9, f"Min: {min_distance:.3f}", fontsize=10)
    plt.text(2, plt.ylim()[1] * 0.85, f"Max: {max_distance:.3f}", fontsize=10)
    plt.text(2, plt.ylim()[1] * 0.8, f"Avg: {avg_distance:.3f}", fontsize=10)

    plt.show()
    

    #visualize encoding data
    print("BioBert Annotation Similarity v1 End")
    
if __name__ == "__main__":
    main()