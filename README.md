# BioBert_Annotations

This project was set up with python 3.11.9 and Visual Studio. If you chose to run this program in a different environment, please disregard the next section. The libraries from the requirements.txt will be needed to run this program. Please install those in any way you prefer.
When this project was started in November 2024, python version 3.11.9 was chosen due to compatibility issues. Depending on the date, these issues might have been resolved.

## Setting up the environment
When Python and Visual Studio Code are installed, we will need to set up a python environment and install the libraries we will use.
In Visual Studio Code, open the folder that contains repository.
In the terminal (Str+ö or View/Terminal) type ``` python -m venv .venv ```

This sets up the python environment.

After that, install the necesary libraries by typeing ``` pip install -r requirements.txt ``` in the terminal.

## Running the program

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


### 4 Generate Encodings


## links
our google doc: https://docs.google.com/document/d/1TD_wkrN5wPjKABs1-eJSVIOPAThNpR0uJz90E76WAeY/edit?tab=t.0
https://huggingface.co/docs/transformers/preprocessing
https://interpro-documentation.readthedocs.io/en/latest/download.html#api

