from dataMethods import *

#init settings
getPfamAnnotations = False
getUniProtAnnotations = False
outputPath = ""
getAnnotationsFromDB = False
annotationFiles = ["C:\\Users\\SebastianRossboeck\\Desktop\\BioBert4\\output_test.txt"]
inputFiles = []
annotationData = []
print("BioBert Annotation Similarity v1 Start") 

#load or get AnnotationData
if annotationFiles != "" and not getAnnotationsFromDB:
    print("loading annotationdata from file")
    annotationData = loadAnnotations(annotationFiles)

if getAnnotationsFromDB and len(inputFiles) > 0:
    print("fetching annotations from api")
    getAnnotations("C:\\Users\\SebastianRossboeck\\Desktop\\BioBert4\\.venv\\spore_formers_proteinIds_test", "C:\\Users\\SebastianRossboeck\\Desktop\\BioBert4\\output_test.txt")
    

for annotation in annotationData:
    print(annotation.pfam_description)

#Data Eval on annotations: put this into a method
evaluatedData = evaluateData(annotationData)
evaluatedData.print_data()
#save to outputPath

#get Encodings for pfam-description
if getPfamAnnotations:
    for protein in annotationData:
        #no, dont do it one by one
        protein.pfam_embedding = ""

#get Encodings for uniprot-function
if getUniProtAnnotations:
    for protein in annotationData:
        #not one by one
        protein.uniprot_embedding = ""

#cluster the encodings data

#visualize encoding data
print("BioBert Annotation Similarity v1 End") 
