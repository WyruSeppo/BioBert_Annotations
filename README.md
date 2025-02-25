# BioBert_Annotations

## Overview
This program was created during the winter semester 2024/25 at the university of vienna for the course "2024W 053531-1 Softwareentwicklungsprojekt Bioinformatik". by Theodora Beshara and Sebastian Rossböck.

This program reads the Ids of proteins from a file in the fasta format and fetches the annotations for each protein from UniProt and Pfam. With BioBert, it creates the embeddings for each annotation and then compares and visualizes the results.

## Compatibility

This project was set up with python 3.11.9 and Visual Studio. If you chose to run this program in a different environment, please disregard the next section. The libraries from the requirements.txt will be needed to run this program. Please install those in any way you prefer.

## Setting up the environment
When Python and Visual Studio Code are installed, we will need to set up a python environment and install the libraries we will use.
In Visual Studio Code, open the folder that contains repository.
In the terminal (Str+ö or View/Terminal) type ``` python -m venv .venv ```

This sets up the python environment.

After that, install the necesary libraries by typeing ``` pip install -r requirements.txt ``` in the terminal.

## Running the program

In general, steps one through five will create a list of AnnotationData objects that collect the relevant information for each protein. The list of objects will be serialized at various parts of the program. With the config, it is possible to set the program up to skip steps that were already completed. Especially fetching the annotations can be a time-consuming task.

### 1. Config
The biobert.ini file contains multiple parameters to configure the program. At the start, these setting will be loaded.

- fasta_file =  The file location for the fasta file that contains the Ids of proteins we want to process
- annotation_file_input =  The file location for the Annotation Data, when loading the data from file.
- annotation_file_output =  The file location for the Annotation Data, without the embeddings.
- annotation_embedding_file_output = The file location for the Annotation Data, including the embeddings.
- data_eval_output = The file location for the result of the data evaluation
- loadAnnotationsFromFile (boolean): Sets whether or not the annotation data should be loaded from a file or fetched from the APIs 
- getPfamEmbeddings (boolean): Sets whether or not the embeddings for the Pfam annotations should be generated
- getUniProtEmbeddings (boolean): Sets whether or not th embeddings for the UniProt annotations should be generated
- model: the name of the model that will be loaded

### 2.1 Loading annotations
If the corresponding config is set, the annotationdata will be loaded from a file. This file will be created after fetching the data from the APIs. The process of fetching the annotations is quite time consuming, since they are being fetched in serial. It is advised to load the data from file once they have been fetched. "annotationdata_output.txt" contains the annotationData for the spore-former.faa dataset.

### 2.2 Fetching annotations
For each Ref_Seq identifer in the fasta file, we convert those to UniProtKB Ids, then get the annotations from UniProt and Pfam via their respective APIs and create a List of AnnotationData objects. These objects hold the relevant information for the nucleodtid sequences, including: UniprotId, RefSeqId, PfamId, annotations from Uniprot and Pfam, as well as the BioBert embeddings later on.

### 3 Evalutate Data

The annotations are evaluated. Their average length, max, min and missing annotations are noted and the result is saved at the data_eval_output-setting. 

### 3.5 Clean Data

For our purpose, we can only use proteins that have Pfam as well as UniProt annotations. Also, those annotations can not exceed 512 words, since the BioBert model would truncate the annotation after that, and we would be dealing with incomplete data.

For these reasons, we remove all AnnotationData objects that do not fit these criteria in this step.

### 4 Generate Encodings

for each annotation, we generate the embedding by passing the data through the model and then saving the embedding in the respective AnnotationData object.
```
inputs = tokenizer(textInput, return_tensors="pt", padding=True, truncation=True, max_length=512)

# Forward pass through the model
with torch.no_grad():
    outputs = model(**inputs)

# Extract embeddings (CLS token representations)
embeddings = outputs.last_hidden_state[:, 0, :]
```
This data is then saved to the path specified in annotation_embedding_file_output. This file is too big to upload it in github, but If you have created it once, you can load it again and skip the first four steps of the program.

### 5 Calculate the distance between embeddings

For each protein, we compute the distance between the embeddings of the UniProt and Pfam annotations using the cosine similarity. Since the embeddings can be thought of as high-dimensional vectors, we can use the cosine similarity to gauge how similary the "direction" of the two "vectors" is.

The data will be saved to the location designated in the biobert.ini file.

## links
our google doc: https://docs.google.com/document/d/1TD_wkrN5wPjKABs1-eJSVIOPAThNpR0uJz90E76WAeY/edit?tab=t.0
https://huggingface.co/docs/transformers/preprocessing
https://interpro-documentation.readthedocs.io/en/latest/download.html#api

