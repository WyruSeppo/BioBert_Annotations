from transformers import AutoTokenizer, AutoModel
import torch

# Load BioBERT model and tokenizer
model_name = "dmis-lab/biobert-base-cased-v1.1"
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModel.from_pretrained(model_name)

# Example sentences
sentences = [
    "Bioinformatics is an interdisciplinary field.",
    "Protein structure prediction is challenging."
]

# Tokenize sentences. maxlength 512?
inputs = tokenizer(sentences, return_tensors="pt", padding=True, truncation=True, max_length=128)

# Forward pass through the model
with torch.no_grad():
    outputs = model(**inputs)

# Extract embeddings (CLS token representations)
embeddings = outputs.last_hidden_state[:, 0, :]

# Convert to numpy
embeddings_numpy = embeddings.numpy()

print(embeddings_numpy.shape) 
print(embeddings_numpy[0])