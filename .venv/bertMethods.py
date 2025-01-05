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
        with torch.no_grad():
            outputs = model(**inputs)

        # Extract embeddings (CLS token representations)
        embeddings = outputs.last_hidden_state[:, 0, :]

        # Convert to numpy
        embeddings_numpy = embeddings.numpy()
        
        if mode == "pfam":
            annotation.pfam_embedding = embeddings_numpy[0]
            
        if mode == "uniprot":
            annotation.uniprot_embedding = embeddings_numpy[0]
        
    return annotationData