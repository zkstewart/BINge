import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Various_scripts import ZS_GFF3IO

from ZS_GFF3IO import GFF3 # expose this to any callers

def _create_feature_from_sl(sl):
    # Extract details from sl
    contig, source, featureType, start, end, \
        score, strand, frame, attributes \
        = sl
    
    # Build feature with main details
    feature = ZS_GFF3IO.Feature()
    feature.add_attributes({
        "contig": contig, "source": source, "type": featureType,
        "start": int(start), "end": int(end), "coords": [int(start), int(end)],
        "score": score, "strand": strand, "frame": frame
    })
    
    # Add attributes
    splitAttributes = []
    for a in attributes.split("="):
        if ";" in a:
            splitAttributes += a.rsplit(";", maxsplit=1)
        else:
            splitAttributes.append(a)
    attributesDict = {splitAttributes[i]: splitAttributes[i+1] for i in range(0, len(splitAttributes)-(len(splitAttributes)%2), 2)}
    feature.add_attributes(attributesDict)
    
    return feature

def _build_feature_to_yield(thisFeature):
    feature = _create_feature_from_sl(thisFeature[0])
    mrnaFeature = _create_feature_from_sl(thisFeature[1])
    feature.add_child(mrnaFeature)
    
    # Add subfeatures to base
    for _sl in thisFeature[2:]:
        _feature = _create_feature_from_sl(_sl)
        mrnaFeature.add_child(_feature)
    
    return feature

def iterate_through_gff3(gff3File):
    '''
    Provides a simple iterator for a GFF3 file which yields GFF3
    features. Only works if the GFF3 file is sorted.
    
    Parameters:
        gff3File -- a string indicating the location of a GFF3 formatted
                    file that is sorted.
    '''
    
    
    thisFeature = []
    with open(gff3File, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\t\r\n ").split("\t")
            
            # Skip comment / irrelevant lines
            if line.startswith("#") or len(sl) != 9: # gff3 lines are always 9 long
                continue
            
            # Handle content lines
            else:
                featureType = sl[2]
                
                # Yield a completed gene feature
                if featureType == "gene" and thisFeature != []:
                    feature = _build_feature_to_yield(thisFeature)
                    
                    # Reset our feature storage, and yield result now
                    thisFeature = []
                    yield feature
                
                # Build an ongoing feature
                thisFeature.append(sl)
    
    # Yield the last feature in the GFF3
    if thisFeature != []:
        feature = _build_feature_to_yield(thisFeature)
    else:
        raise Exception(f"'{gff3File}' does not appear to be a valid GFF3 file!")
    yield feature
