import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ZS_GFF3IO import GFF3 # expose this to any callers

def _build_data_to_yield(thisFeature):
    # Parse out variables from mRNA line
    contig, _, _, start, end, \
        _, strand, _, attributes \
        = thisFeature[1]
    
    # Start building a data structure to yield
    dataDict = {"contig": contig, "start": int(start), "end": int(end), "strand": strand}
    
    # Get the GMAP alignment quality attributes
    for attribute in attributes.split(";"):
        key, value = attribute.split("=")
        if key in ["identity", "coverage", "indels"]:
            dataDict[key] = value
        elif key == "Name":
            dataDict["Name"] = value
    
    # Get the exons
    dataDict["exons"] = []
    for _sl in thisFeature[2:]:
        contig, _, _, start, end, \
            _, strand, _, attributes \
            = _sl
        dataDict["exons"].append([int(start), int(end)])
    dataDict["exons"] = sorted(dataDict["exons"])
    
    return dataDict

def iterate_gmap_gff3(gff3File):
    '''
    Provides a simple iterator for a GMAP GFF3 file which yields only the
    relevant details for BINge. Assumes the file is sorted as is the case
    with GMAP output.
    
    Parameters:
        gff3File -- a string indicating the location of a GMAP GFF3 formatted
                    file that is sorted.
    '''
    thisAlignment = []
    with open(gff3File, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\t\r\n ").split("\t")
            
            # Skip comment / irrelevant lines
            if line.startswith("#") or len(sl) != 9: # gff3 lines are always 9 long
                continue
            
            # Handle content lines
            else:
                featureType = sl[2]
                
                # Yield a completed gene data dictionary
                if featureType == "gene" and thisAlignment != []:
                    dataDict = _build_data_to_yield(thisAlignment)
                    
                    # Reset our feature storage, and yield result now
                    thisAlignment = []
                    yield dataDict
                
                # Build an ongoing feature
                thisAlignment.append(sl)
    
    # Yield the last feature in the GFF3
    if thisAlignment != []:
        dataDict = _build_data_to_yield(thisAlignment)
    else:
        raise Exception(f"'{gff3File}' does not appear to be a valid GFF3 file!")
    yield dataDict
