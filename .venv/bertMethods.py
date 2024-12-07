from transformers import AutoTokenizer, AutoModel
import torch

def getEmbeddings(annotationData, modelName = "dmis-lab/biobert-base-cased-v1.1", mode="pfam"):

    # Load BioBERT model and tokenizer
    tokenizer = AutoTokenizer.from_pretrained(modelName)
    model = AutoModel.from_pretrained(modelName)

    #this is for testing! do it in batches
    for annotation in annotationData:
        # Tokenize sentences. maxlength 512?
        inputs = tokenizer(annotation.pfam_description, return_tensors="pt", padding=True, truncation=True, max_length=512)

        # Forward pass through the model
        with torch.no_grad():
            outputs = model(**inputs)

        # Extract embeddings (CLS token representations)
        embeddings = outputs.last_hidden_state[:, 0, :]

        # Convert to numpy
        embeddings_numpy = embeddings.numpy()
        annotation.pfam_embedding = embeddings_numpy[0]
        
    return annotationData