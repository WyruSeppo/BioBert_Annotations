from transformers import AutoTokenizer, AutoModelForQuestionAnswering
import torch

# Load BioBERT tokenizer and model
tokenizer = AutoTokenizer.from_pretrained("dmis-lab/biobert-base-cased-v1.1-squad")
model = AutoModelForQuestionAnswering.from_pretrained("dmis-lab/biobert-base-cased-v1.1-squad")

# Define a sample question and context
question = "What is the function of hemoglobin?"
context = "Hemoglobin is a protein in red blood cells that carries oxygen throughout the body."


# Tokenize inputs
inputs = tokenizer(question, context, return_tensors="pt")

# Perform inference
with torch.no_grad():
    outputs = model(**inputs)

# Extract answer tokens and convert back to text
answer_start = torch.argmax(outputs.start_logits)
answer_end = torch.argmax(outputs.end_logits) + 1
answer = tokenizer.convert_tokens_to_string(tokenizer.convert_ids_to_tokens(inputs["input_ids"][0][answer_start:answer_end]))

print("Question:", question)
print("Answer:", answer)
