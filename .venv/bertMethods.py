from transformers import AutoTokenizer, AutoModel
import torch
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("BioBERT")

def getEmbeddings(annotationData, modelName = "dmis-lab/biobert-base-cased-v1.1", mode="pfam"):

    # Load BioBERT model and tokenizer
    tokenizer = AutoTokenizer.from_pretrained(modelName)
    model = AutoModel.from_pretrained(modelName)

    count = 0
    max = len(annotationData)
 
    for annotation in annotationData:
        count += 1        
        logger.info(str(count / max) + " %")
        
        # Tokenize sentences. maxlength 512?
        textInput = annotation.pfam_description if mode == "pfam" else annotation.uniprot_function
        inputs = tokenizer(textInput, return_tensors="pt", padding=True, truncation=True, max_length=512)

        # Forward pass through the model
        #the token IDs are embedded into vector space
        #the embeddings are processed by BioBerts transformers
        with torch.no_grad():
            outputs = model(**inputs)

        # Extract embeddings (CLS token representations)
        #the "last_hidden_state" contains the output of the last transformer layer for each token in the sequence
        # Extract embeddings (CLS token representations)
        embeddings = outputs.last_hidden_state[:, 0, :]

        # Convert to numpy
        embeddings_numpy = embeddings.numpy()

        # Flatten to 1D
        embeddings_flat = embeddings_numpy.flatten()

        # Save as space-separated string (no newlines)
        embedding_str = " ".join(map(str, embeddings_flat))

        # Save the embeddings in a text file
        if mode == "pfam":
            annotation.pfam_embedding = embedding_str  # Save as space-separated string
            
        if mode == "uniprot":
            annotation.uniprot_embedding = embedding_str  # Save as space-separated string

        
    return annotationData