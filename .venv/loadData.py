#load all annotatios from output file
#source: pfam = 3, uniprot = 2
def loadAnnotations(fileName="", source=3):
    annotations = []
    file = open(fileName, "r")

    for line in file:
        annotations.append(line.split("\t")[source])

    return annotations