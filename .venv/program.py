from dataMethods import *
from bertMethods import *

#init settings
getPfamEmbeddings = True
getUniProtEmbeddings = False
outputPath = ""
getAnnotationsFromDB = False
annotationFile =["C:\\Users\\SebastianRossboeck\\Desktop\\BioBert4\\output_test.txt"]
inputFile = ""
annotationData = []
fastaFile = "C:\\Users\\SebastianRossboeck\\Desktop\\BioBert4\\.venv\\spore-formers.faa"
annotationData_new = []

print("BioBert Annotation Similarity v1 Start") 

#load FASTA file
#for record in SeqIO.parse(fastaFile, "fasta"):
#    annotationData_new.append(AnnotationData("","","","","",record.id,"","","","","","",""))

#get UniProt IDs from RefSeq
#test = getUniProtConversion("RefSeq_Protein","UniProtKB",annotationData_new)


#load or get AnnotationData
if annotationFile != "" and not getAnnotationsFromDB:
    print("loading annotationdata from file")
    annotationData = loadAnnotations(annotationFile)

#no input file was given: get Annotations from API
if getAnnotationsFromDB and inputFile != "":
    print("fetching annotations from api")
    getAnnotations("C:\\Users\\SebastianRossboeck\\Desktop\\BioBert4\\.venv\\spore_formers_proteinIds_test", "C:\\Users\\SebastianRossboeck\\Desktop\\BioBert4\\output_test.txt")
    
for annotation in annotationData:
    print(annotation.pfam_description)

#Data Eval on annotations
evaluatedData = evaluateData(annotationData)
evaluatedData.print_data()
#save to outputPath

#get Encodings for pfam-description
if getPfamEmbeddings:
    annotationData = getEmbeddings(annotationData)
    
for x in annotationData:
    print(x.pfam_embedding)

#get Encodings for uniprot-function
if getUniProtEmbeddings:
    for protein in annotationData:
        #not one by one
        protein.uniprot_embedding = ""
 
#cluster the encodings data

#eine änderung

#visualize encoding data
print("BioBert Annotation Similarity v1 End") 
