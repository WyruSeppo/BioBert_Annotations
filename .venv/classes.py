#the class AnnotationData bundles all the data we are handling for each protein
class AnnotationData:
    """
    Represents annotation data for a protein, including Pfam and UniProt annotations.
    
    Attributes:
        id (int): Unique identifier for the annotation.
        pfam_id (str): Pfam identifier for the protein.
        uniprot_id (str): UniProt identifier for the protein.
        pfam_embedding: Embedding representation of the Pfam annotation.
        uniprot_embedding: Embedding representation of the UniProt annotation.
        refSeqAccession (str): RefSeq accession number.
        entry (str): Entry name from UniProt.
        entry_name (str): Alternative entry name.
        protein_names (str): Names of the protein.
        gene_names (str): Names of the associated genes.
        organism (str): Organism from which the protein originates.
        pfam_description (str): Description of the Pfam annotation.
        uniprot_function (str): Functional description from UniProt.
        embedding_distance: Distance metric between Pfam and UniProt embeddings.
    """
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
        """Returns a string representation of the AnnotationData object."""
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
        """Converts the AnnotationData object into a dictionary."""
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
        
    def cleanAnnotations(self):
        """Cleans annotation descriptions by removing excess whitespace."""
        #remove leading/trailing whitespaces
        self.pfam_description = self.pfam_description.strip()
        
        # Remove multiple whitespaces
        self.pfam_description = ' '.join(self.pfam_description.split())

class EvaluatedData:
    """
    Stores statistical evaluations of annotation data, such as counts and lengths of annotations.
    
    Attributes:
        no_sequences (int): Number of sequences analyzed.
        uniprot_annotation_amount (int): Count of UniProt annotations.
        uniprot_annotation_length_min (int): Minimum length of UniProt annotations.
        uniprot_annotation_length_max (int): Maximum length of UniProt annotations.
        uniprot_annotation_length_avg (int): Average length of UniProt annotations.
        uniprot_annotation_no_words (int): Number of words in UniProt annotations.
        uniprot_annotation_missing_amount (int): Number of missing UniProt annotations.
        uniprot_annotation_missing_percent (float): Percentage of missing UniProt annotations.
        pfam_annotation_amount (int): Count of Pfam annotations.
        pfam_annotation_length_min (int): Minimum length of Pfam annotations.
        pfam_annotation_length_max (int): Maximum length of Pfam annotations.
        pfam_annotation_length_avg (int): Average length of Pfam annotations.
        pfam_annotation_no_words (int): Number of words in Pfam annotations.
        pfam_annotation_missing_amount (int): Number of missing Pfam annotations.
        pfam_annotation_missing_percent (float): Percentage of missing Pfam annotations.
    """
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
        """Prints a formatted summary of the evaluation data."""
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
        """Returns a formatted string summary of the evaluation data."""
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