from dataMethods import *
from Bio import SeqIO
import time
import requests
import pandas as pd
from classes import *
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("BioBERT")

def get_uniprot_annotation(protein_id):
    url = f"https://www.uniprot.org/uniprotkb/{protein_id}.txt"
    try:
        response = requests.get(url)
    except Exception as e:
        print(f"An error occurred while calling uniprot API: {protein_id}")
        return None, None, None

    if response.status_code == 200:
        text_data = response.text
        protein_name = None
        function = None
        pfam = None
        endOfFunction = True
        
        for line in text_data.splitlines():
            if line.startswith("ID"):
                protein_name = line.split()[1]
            elif line.startswith("CC   -!- FUNCTION:"):
                endOfFunction = False
                function = line[len("CC   -!- FUNCTION:"):].strip()
            elif line.__contains__("-!-") and not endOfFunction:
                endOfFunction = True
            elif line.startswith("CC") and  not endOfFunction:
                function += line[2:].lstrip() #remove whitespace and stuff at the end
            elif line.startswith("DR   Pfam;"):
                pfam = line.split()[2].removesuffix(';')

        return protein_name, function, pfam
    else:
        return None, None, None

def get_pfam_annotation(protein_id):
    api_url = "https://www.ebi.ac.uk/interpro/api"
    url = f"{api_url}/entry/pfam/{protein_id}"
    url += "?page_size=200"

    try:
        response = requests.get(url)
    except Exception as e:
        print(f"An error occurred while calling pfam API: {protein_id}")
        return None

    if response.status_code == 200:
        description_text = None
        jsonData = response.json()

        if jsonData['metadata']['description'] and isinstance(jsonData['metadata']['description'], list) and jsonData['metadata']['description'] != None and len(jsonData['metadata']['description']) > 0:
            description_text = jsonData['metadata']['description'][0]['text']
            description_text = description_text.removesuffix('</p>')
            description_text = description_text.removeprefix('<p>')
        return description_text
    else:
        return None

#no longer in use, but we might need it
'''
def annotate_tsv(tsv_file):
    df = read_tsv(tsv_file)
    annotated_data = []
    counter = 0
    total = 0

    if df is not None:
        total = len(df)
        for protein_id in df.iloc[:, 1]:
            print("\r" + str(counter / total) + "% " + str(counter) + "/" + str(total))
            
            protein_name, function, pfamID = get_uniprot_annotation(protein_id)
            pfamAnnotation = get_pfam_annotation(pfamID)

            annotated_data.append({
                'Protein_ID': protein_id,
                'Pfam_ID':pfamID,
                'Protein_Name': protein_name,
                'Unitprot_Function': function,
                'Pfam_Description': pfamAnnotation,
                'RefSeq':df.iloc[counter,0],
                'Entry':df.iloc[counter,1],
                'Entry_Name':df.iloc[counter,3],
                'Protein_Name':df.iloc[counter,4],
                'Gene_Names':df.iloc[counter,5],
                'Organism':df.iloc[counter,6],
                
            })
            counter += 1
            #print((protein_id or 'NA') + " " + (pfamID or 'NA') + " " + (protein_name or 'NA') + " " + (function or 'NA') + " " + (pfamAnnotation or 'NA'))
    return pd.DataFrame(annotated_data)
'''

def annotate_data(annotationData, showProgress = False):
    counter = 0
    total = len(annotationData)

    for entry in annotationData:
        print("\r" + str(counter / total) + "% " + str(counter) + "/" + str(total))
        
        protein_name, function, pfamID = get_uniprot_annotation(entry.uniprot_id)
        pfamAnnotation = get_pfam_annotation(pfamID)

        entry.pfam_id = pfamID or "None"
        entry.protein_names = protein_name or "None"
        entry.uniprot_function = function or "None"
        entry.pfam_description = pfamAnnotation or "None"

        # Set other attributes to "None" if needed
        for attr in ['id', 'pfam_embedding', 'uniprot_embedding', 'entry', 'entry_name', 'gene_names', 'organism', 'embedding_distance']:
            setattr(entry, attr, "None")

        counter += 1
        
    return annotationData
    
def getUniProtConversion(ffrom,to,refseqIds):

    ids = ",".join(refseqIds)
    url = "https://rest.uniprot.org/idmapping/run"
    data = {
        "ids": ids,
        "from": ffrom,
        "to": to
    }

    response = requests.post(url, data=data)

    #error handling
    print(response.status_code)
    jobId = response.json()['jobId']
    jobDone = False
    
    while(not jobDone):
        url = "https://rest.uniprot.org/idmapping/status/" + str(jobId)
        response = requests.get(url)
        status = response.json().get('jobStatus')
        
        if(status):
            print(status)
            jobDone = status == "FINISHED"
        else:
            jobDone = True
        # Wait two minutes
        time.sleep(60)


    #get all the results form idmapping job

    print("job done")
    #put UniprotIds into annotationData 
    result = []
    #for entry in response.json()["results"]:
    #   result.append(AnnotationData("","",entry["to"],"","",entry["from"],"","","","","","",""))
    
    # Base URL for the UniProt API
    base_url = "https://rest.uniprot.org/idmapping/results/"

    # Your job ID
    job_id = jobId

    # Pagination settings
    page_size = 500  # Number of results per page

    # Collect all results
    all_results = []
    url = f"{base_url}{job_id}"

    while url:
        # Make the request
        response = requests.get(url)
        response.raise_for_status()
        
        # Parse the JSON response
        data = response.json()
        all_results.extend(data.get("results", []))
        
        # Get the next page URL from the headers
        url = response.headers.get("Link")
        if url:
            # The "Link" header may contain multiple URLs, ensure you parse the `rel="next"` URL
            links = [link.strip() for link in url.split(",")]
            next_link = [link for link in links if 'rel="next"' in link]
            if next_link:
                # Extract the actual URL (it will be enclosed in < >)
                url = next_link[0].split(";")[0].strip("<>")
            else:
                url = None  # No "next" link available
        else:
            url = None  # No "Link" header

        # Use the collected results
        print(f"Fetched {len(all_results)} results")
        
    for entry in all_results:
       result.append(AnnotationData("","",entry["to"],"","",entry["from"],"","","","","","","",""))
    
    
    return result    