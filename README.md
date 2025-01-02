# BioBert_Annotations

This project was set up with python 3.11.9 and Visual Studio. These will need to be installed for this guide to work.
When this project was started in November 2024, python version 3.11.9 was chosen due to compatibility issues. Depending on the date, these issues might have been resolved.

## Setting up the environment
When Python and Visual Studio Code are installed, we will need to set up a python environment and install the libraries we will use.
In Visual Studio Code, open the folder that contains repository.
In the terminal (Str+ö or View/Terminal) type ``` python -m venv .venv ```

This sets up the python environment.

After that, install the necesary libraries by typeing ``` pip install -r requirements.txt ``` in the terminal.

# Running the program

## 1. Config
The biobert.ini file contains multiple parameters to configure the program.

- fasta_file =  The file location for the fasta file that contains the Ids of proteins we want to process
- annotation_file_input =  The file location for the Annotation Data, when loading the data from file.
- annotation_file_output =  The file location for the Annotation Data, without the embeddings.
- annotation_embedding_file_output = The file location for the Annotation Data, including the embeddings.
- data_eval_output = The file location for the result of the data evaluation
- loadAnnotationsFromFile (boolean): Sets whether or not the annotation data should be loaded from a file or fetched from the APIs 
- getPfamEmbeddings (boolean): Sets whether or not the embeddings for the Pfam annotations should be generated
- getUniProtEmbeddings (boolean): Sets whether or not th embeddings for the UniProt annotations should be generated
- model: the name of the model that will be loaded

## 2.1 Loading annotations

## 2.2 Fetching annotations

## 3 Evalutate Data

## 4 Generate Encodings


## links
our google doc: https://docs.google.com/document/d/1TD_wkrN5wPjKABs1-eJSVIOPAThNpR0uJz90E76WAeY/edit?tab=t.0
https://huggingface.co/docs/transformers/preprocessing
https://interpro-documentation.readthedocs.io/en/latest/download.html#api

