#the class AnnotationData bundles all the data we are handling for each protein
class AnnotationData:
    def __init__(self, id: int, pfam_id: str, uniprot_id: str, pfam_embedding, uniprot_embedding, refSeqAccession, entry, entryName, proteinNames, geneNames, organism, pfam_description,uniprot_function):
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
        
    def __repr__(self):
        return f"AnnotationData(id={self.id}, pfam_id='{self.pfam_id}', uniprot_id='{self.uniprot_id}')"

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
