from dataMethods import *
import os

tsv_file_path = "spore_formers_proteinIds_1"
tsv_file_path2 = "spore_formers_proteinIds_2"

base_dir = os.path.dirname(__file__)
file_path = os.path.join(base_dir, "spore_formers_proteinIds_1")
file_path2 = os.path.join(base_dir, "spore_formers_proteinIds_2")

annotated_df = annotate_tsv(file_path)
annotated_df.to_csv('output1.txt', sep='\t', index=False)

annotated_df = annotate_tsv(file_path2)
annotated_df.to_csv('output2.txt', sep='\t', index=False)