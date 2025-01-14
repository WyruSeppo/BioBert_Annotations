#the class AnnotationData bundles all the data we are handling for each protein
class AnnotationData:
    def __init__(self, id: int, pfam_id: str, uniprot_id: str, pfam_embedding, uniprot_embedding, refSeqAccession, entry, entryName, proteinNames, geneNames, organism, pfam_description,uniprot_function, embedding_distance):
        self.id = id
        self.pfam_id = pfam_id
        self.uniprot_id = uniprot_id
        self.pfam_embedding = pfam_embedding
        self.uniprot_embedding = uniprot_embedding
        self.refSeqAccession = refSeqAccession
        self.entry = entry
        self.entry_name = entryName
        self.protein_names = proteinNames
        self.gene_names = geneNames
        self.organism = organism
        self.uniprot_function = uniprot_function
        self.pfam_description = pfam_description
        self.embedding_distance = embedding_distance
        
    def __repr__(self):
        return (f"AnnotationData("
                f"id={self.id}, "
                f"pfam_id='{self.pfam_id}', "
                f"uniprot_id='{self.uniprot_id}', "
                f"pfam_embedding={self.pfam_embedding}, "
                f"uniprot_embedding={self.uniprot_embedding}, "
                f"refSeqAccession='{self.refSeqAccession}', "
                f"entry='{self.entry}', "
                f"entry_name='{self.entry_name}', "
                f"protein_names='{self.protein_names}', "
                f"gene_names='{self.gene_names}', "
                f"organism='{self.organism}', "
                f"pfam_description='{self.pfam_description}', "
                f"uniprot_function='{self.uniprot_function}')")
        
    def to_dict(self):
        return {
            "id": self.id,
            "pfam_id": self.pfam_id,
            "uniprot_id": self.uniprot_id,
            "pfam_embedding": self.pfam_embedding,
            "uniprot_embedding": self.uniprot_embedding,
            "refSeqAccession": self.refSeqAccession,
            "entry": self.entry,
            "entry_name": self.entry_name,
            "protein_names": self.protein_names,
            "gene_names": self.gene_names,
            "organism": self.organism,
            "pfam_description": self.pfam_description,
            "uniprot_function": self.uniprot_function,
            "embedding_distance": self.embedding_distance,
        }
        
    def cleanAnnotations():
        for i in range(len(self.pfam_description)):
            #remove leading/trailing whitespaces
            self.pfam_description[i] = self.pfam_description[i].strip()
            
            # Remove multiple whitespaces
            self.pfam_description[i] = ' '.join(self.pfam_description[i].split())

class EvaluatedData:
    def __init__(self):
        self.no_sequences = 0
        self.uniprot_annotation_amount = 0
        self.uniprot_annotation_length_min = 0
        self.uniprot_annotation_length_max = 0
        self.uniprot_annotation_length_avg = 0
        self.uniprot_annotation_no_words = 0
        self.uniprot_annotation_missing_amount = 0
        self.uniprot_annotation_missing_percent = 0
        self.pfam_annotation_amount = 0
        self.pfam_annotation_length_min = 0
        self.pfam_annotation_length_max = 0
        self.pfam_annotation_length_avg = 0
        self.pfam_annotation_no_words = 0
        self.pfam_annotation_missing_amount = 0
        self.pfam_annotation_missing_percent = 0

    def print_data(self):
        print("Evaluated Data:")
        print(f"  Number of Sequences: {self.no_sequences}")
        print(f"  UniProt Annotations:")
        print(f"    Amount: {self.uniprot_annotation_amount}")
        print(f"    Annotation-Length (Min/Max/Avg): {self.uniprot_annotation_length_min} / {self.uniprot_annotation_length_max} / {self.uniprot_annotation_length_avg}")
        print(f"    Number of Words: {self.uniprot_annotation_no_words}")
        print(f"    Missing Amount: {self.uniprot_annotation_missing_amount}")
        print(f"    Missing Percent: {self.uniprot_annotation_missing_percent}%")
        print(f"  Pfam Annotations:")
        print(f"    Amount: {self.pfam_annotation_amount}")
        print(f"    Annotation-Length (Min/Max/Avg): {self.pfam_annotation_length_min} / {self.pfam_annotation_length_max} / {self.pfam_annotation_length_avg}")
        print(f"    Number of Words: {self.pfam_annotation_no_words}")
        print(f"    Missing Amount: {self.pfam_annotation_missing_amount}")
        print(f"    Missing Percent: {self.pfam_annotation_missing_percent}%")

    def __repr__(self):
        return (
            f"EvaluatedData("
            f"no_sequences={self.no_sequences}, "
            f"uniprot_annotation_amount={self.uniprot_annotation_amount}, "
            f"uniprot_annotation_length_min={self.uniprot_annotation_length_min}, "
            f"uniprot_annotation_length_max={self.uniprot_annotation_length_max}, "
            f"uniprot_annotation_length_avg={self.uniprot_annotation_length_avg}, "
            f"uniprot_annotation_no_words={self.uniprot_annotation_no_words}, "
            f"uniprot_annotation_missing_amount={self.uniprot_annotation_missing_amount}, "
            f"uniprot_annotation_missing_percent={self.uniprot_annotation_missing_percent}, "
            f"pfam_annotation_amount={self.pfam_annotation_amount}, "
            f"pfam_annotation_length_min={self.pfam_annotation_length_min}, "
            f"pfam_annotation_length_max={self.pfam_annotation_length_max}, "
            f"pfam_annotation_length_avg={self.pfam_annotation_length_avg}, "
            f"pfam_annotation_no_words={self.pfam_annotation_no_words}, "
            f"pfam_annotation_missing_amount={self.pfam_annotation_missing_amount}, "
            f"pfam_annotation_missing_percent={self.pfam_annotation_missing_percent}"
            f")"
        )
        
    def generateString(self):
        return (
            f"Evaluated Data:\n"
            f"  Number of Sequences: {self.no_sequences}\n"
            f"  UniProt Annotations:\n"
            f"    Amount: {self.uniprot_annotation_amount}\n"
            f"    Annotation-Length (Min/Max/Avg): "
            f"{self.uniprot_annotation_length_min} / "
            f"{self.uniprot_annotation_length_max} / "
            f"{self.uniprot_annotation_length_avg}\n"
            f"    Number of Words: {self.uniprot_annotation_no_words}\n"
            f"    Missing Amount: {self.uniprot_annotation_missing_amount}\n"
            f"    Missing Percent: {self.uniprot_annotation_missing_percent}%\n"
            f"  Pfam Annotations:\n"
            f"    Amount: {self.pfam_annotation_amount}\n"
            f"    Annotation-Length (Min/Max/Avg): "
            f"{self.pfam_annotation_length_min} / "
            f"{self.pfam_annotation_length_max} / "
            f"{self.pfam_annotation_length_avg}\n"
            f"    Number of Words: {self.pfam_annotation_no_words}\n"
            f"    Missing Amount: {self.pfam_annotation_missing_amount}\n"
            f"    Missing Percent: {self.pfam_annotation_missing_percent}%"
        )