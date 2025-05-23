import os, requests, re
from itertools import groupby
from goatools import obo_parser

from .go import fix_obsoletions

def init_table(repToClust, blastFile, evalueCutoff, numHits, outputFileName,
               databaseTag=".", largeTable=False):
    '''
    Parses an outfmt6 file, initialising our output annotation table with the BLAST
    results for the best numHits values. The input BLAST file can be sorted or
    unsorted (in terms of E-values) as long as query sequences and their results
    appear in blocks. Randomly assorted files will fail horribly.
    
    It will also hold onto the target IDs for these best hits, generating a dictionary
    with keys representing values we want to get ID mappings for, and empty (None)
    values to be set when parsing the idmapping_selected.tab file.
    
    Parameters:
        repToClust -- a dictionary with representative sequence IDs as keys and cluster IDs
                      as values.
        blastFile -- a string indicating the location of a BLAST result file in
                     outfmt6 format.
        evalueCutoff -- a float indicating the maximum allowed E-value for a hit to be
                        parsed and stored.
        numHits -- an integer indicating the maximum number of hits to be stored
                   (if more are available).
        outputFileName -- a string indicating the location to write the annotation
                          TSV to.
        databaseTag -- a string indicating what database the hits come from (if any).
        largeTable -- a boolean indicating whether all outfmt6 fields should be
                      written to file or just the "important" ones i.e., identity,
                      E-value, and bitscore.
    Returns:
        hitMapDict -- a dictionary with structure like:
                      {
                          'tid1': None,
                          'tid2': None,
                          ...
                      }
    '''
    hitMapDict = {}
    
    # Iterate through the blastTab file looking at each group of hits to a single query sequence
    grouper = lambda x: x.split("\t")[0]
    with open(blastFile, 'r') as fileIn, open(outputFileName, "w") as fileOut:
        # Write header line for output
        if largeTable == True:
            fileOut.write("#{0}\n".format("\t".join([
                "Query", "Source", "Target_accession", "Percentage_identity", "Alignment_length",
                "Mismatches", "Gap_opens", "Query_start", "Query_end", "Target_start",
                "Target_end", "Expect_value", "Bit_score"
            ])))
        else:
            fileOut.write("#{0}\n".format("\t".join([
                "Query", "Source", "Target_accession", "Percentage_identity", 
                "Expect_value", "Bit_score"
            ])))
        
        # Continue iterating through the blast file
        for queryID, value in groupby(fileIn, grouper):
            # Skip if this sequence isn't a representative sequence
            if queryID not in repToClust:
                continue
            
            # Get the values in sorted order
            value = sorted(
                [ v.rstrip("\r\n ").split("\t") for v in value ],
                key = lambda x: (float(x[10]), -float(x[11]))
            )
            
            # Get the numHits best hits for this query (that also pass E-value cutoff)
            bestHits = [ v for v in value if float(v[10]) <= evalueCutoff ][0: numHits]
            
            # Fix up UniRef##_ prefixes if necessary
            for i in range(len(bestHits)):
                if bestHits[i][1].startswith("UniRef"): # [1] column in outfmt6 == target
                    bestHits[i][1] = bestHits[i][1].split("_", maxsplit=1)[1]
            
            # Format hits for output
            formattedList = []
            
            ## > Large table formatting
            if largeTable == True:
                # Encapsulate things in square brackets if they're not the first hit
                for i in range(len(bestHits)):
                    if i == 0:
                        for x in range(1, 12): # there's 11 columns in the initial output
                            formattedList.append([bestHits[i][x]])
                    else:
                        for x in range(1, 12):
                            formattedList[x-1].append('[' + bestHits[i][x] + ']')
                
                # Format values into a string
                spacer = " " if len(bestHits) > 1 else "" # add a space inbetween first and second values if relevant
                for i in range(len(formattedList)):
                    formattedList[i] = ''.join([formattedList[i][0], spacer, *formattedList[i][1:]])
                
                # Handle database tagging and fix blank formattedList if no hits were found
                if formattedList == []:
                    hitTag = "."
                    formattedList = ["."]*11
                else:
                    hitTag = databaseTag
            
            ## > Slim table formatting
            else:
                # Encapsulate things in square brackets if they're not the first hit
                slimTargetIndices = [1, 2, 10, 11]
                for i in range(len(bestHits)):
                    if i == 0:
                        for x in slimTargetIndices:
                            formattedList.append([bestHits[i][x]])
                    else:
                        for x in slimTargetIndices:
                            formattedList[slimTargetIndices.index(x)].append('[' + bestHits[i][x] + ']')
                
                # Format values into a string
                spacer = " " if len(bestHits) > 1 else "" # add a space inbetween first and second values if relevant
                for i in range(len(formattedList)):
                    formattedList[i] = ''.join([formattedList[i][0], spacer, *formattedList[i][1:]])
                
                # Handle database tagging and fix blank formattedList if no hits were found
                if formattedList == []:
                    hitTag = "."
                    formattedList = ["."]*4
                else:
                    hitTag = databaseTag
            
            # Write to file
            clusterID = repToClust[queryID] # substitute queryID (representative sequence) with cluster ID
            fileOut.write("{0}\t{1}\t{2}\n".format(clusterID, hitTag,"\t".join(formattedList)))
            
            # Store values in our hitMapDict
            for hit in bestHits:
                hitMapDict.setdefault(hit[1], None) # hit[1] corresponds to target ID
    return hitMapDict

def parse_idmap(idmappingFile, hitMapDict):
    '''
    Parameters:
        idmappingFile -- a string indicating the location of the idmapping_selected.tab
                         file for parsing.
        hitMapDict -- a dictionary with structure like:
                      {
                          'dbSeqID1': None,
                          'dbseqID2': None,
                          ...
                      }
    Modifies:
        hitMapDict -- the dictionary will now have values instead of None, with structure like:
                      {
                          'dbSeqID1': 'GO:0046782'
                          'dbseqID2': 'GO:0033644; GO:0016020',
                          'dbseqID3': '.',
                          'dbseqID4': None, # if this was not in the .tab file it will be None
                          ...
                      }
    '''
    with open(idmappingFile, "r") as fileIn:
        for line in fileIn:
            sl = line.split("\t")
            if len(sl) < 10:
                raise ValueError("The idmapping file is not formatted as expected; it should have " + 
                                 "at least 10 tab-delimited columns.")
            
            # Extract relevant data from table
            uniprotKbAccession = sl[0]
            refseqAccession = sl[3].split("; ") # can have multiple
            uniref100Accession = sl[7].split('_', maxsplit=1)[1] if sl[7] != "" else None
            uniref90Accession = sl[8].split('_', maxsplit=1)[1] if sl[8] != "" else None
            uniref50Accession = sl[9].split('_', maxsplit=1)[1] if sl[9] != "" else None
            
            geneOntologies = sl[6] if sl[6] != "" else "."
            
            # Associate gene ontologies to any IDs in the idmapDict
            for id in set([uniprotKbAccession, *refseqAccession, uniref100Accession, \
            uniref90Accession, uniref50Accession]):
                if id in hitMapDict:
                    hitMapDict[id] = geneOntologies

def update_table_with_gos(originalTable, newTable, hitMapDict, goObo):
    '''
    Updates the originalTable file to include three additional columns appended
    onto the file. The first is a column indicating the best hit with an ID
    mapping, the second are the GOs obtained from that ID mapping, and the third
    are those GOs expanded to include ancestor terms.
    
    This function will attempt to fix GO term obsoletions by use of the EMBL-EBI
    QuickGO API.
    
    Parameters:
        originalTable -- a string indicating the location of the file table which
                         will be updated herein.
        newTable -- a string indicating the location of the new updated table
                    which we will write to.
        hitMapDict -- a dictionary resulting from init_table() and parse_idmap()
                      with structure like:
                      {
                          'dbSeqID1': 'GO:0046782'
                          'dbseqID2': 'GO:0033644; GO:0016020',
                          'dbseqID3': '.',
                          'dbseqID4': None,
                          ...
                      }
        goObo -- a obo_parser.GODag object which has parsed a go.obo file.
    '''
    queriedGOs = {}
    
    with open(originalTable, "r") as fileIn, open(newTable, "w") as fileOut:
        for line in fileIn:
            # Handle header line
            if line.startswith("#"):
                fileOut.write(line.rstrip("\r\n ") + "\t{0}\n".format("\t".join([
                    "Best_hit_with_idmapping", "Best_mapped_GOs",
                    "Best_mapped_GOs_+_parents"
                ])))
            
            # Handle content lines
            else:
                # Extract relevant information from this line
                sl = line.rstrip("\r\n ").split("\t")
                hitAccessions = sl[2]
                
                # Skip any handling if there's no hit
                if hitAccessions == ".":
                    outputValues = ["."]*3
                # Find the best ID with an ID mapping (if any)
                else:
                    hitList = hitAccessions.replace(" ", "").replace("]", "").split("[")
                    mappedHits = [ hit for hit in hitList if hitMapDict[hit] != None ]
                    
                    # If no hits were found, handle that
                    if mappedHits == []:
                        outputValues = ["."]*3
                    # If there were hits, get the best one and process it further
                    else:
                        bestHitID = mappedHits[0]
                        bestHitGOs = [ bhg.strip(" ") for bhg in hitMapDict[bestHitID].split("; ") ]
                        bestHitGOs = [ bhg for bhg in bestHitGOs if bhg != "" and bhg != "." ]
                        
                        # Try to fix any obsoletions (if any)
                        fixedGOs = fix_obsoletions(bestHitGOs, goObo, queriedGOs)
                        
                        # If we ended up with no hits, handle that
                        if fixedGOs == []:
                            outputValues = ["."]*3 # technically we have a bestHitID, but it led nowhere
                        
                        # Otherwise, try to get ancestor terms
                        else:
                            ancestorGOs = set(fixedGOs)
                            for term in fixedGOs:
                                try:
                                    "This may fail if we encounter a 'not really obsoleted' condition during API query"
                                    ancestorGOs = ancestorGOs.union(goObo[term].get_all_parents())
                                except:
                                    pass
                            
                            # Format our output values
                            outputValues = [bestHitID, "; ".join(fixedGOs), "; ".join(ancestorGOs)]
                
                # Write to file
                fileOut.write("{0}\t{1}\n".format("\t".join(sl), "\t".join(outputValues)))

def update_table_with_seq_details(originalTable, newTable, hitMapDict, 
                                  representativesFasta, blastFasta):
    '''
    Receives a table file which is to be modified with sequence details as obtained from
    the FASTA file that was BLASTed against. That FASTA file is expected to have RefSeq
    sequence formatting like:
    
        >NP_001158986.1 uncharacterized protein LOC100303959 [Zea mays]
    
    Or UniRef sequence formatting like:

        >UniRef90_A0A5A9P0L4 peptidylprolyl isomerase n=1 Tax=Triplophysa tibetana TaxID=1572043 RepID=A0A5A9P0L4_9TELE
    
    Any other sequence formats will be assumed to be the reference genome sequence and
    they won't have any details like gene names obtained, because it's impossible to
    anticipate every possible FASTA IDs format.
    
    Parameters:
        originalTable -- a string indicating the location of the file table which
                         will be updated herein.
        newTable -- a string indicating the location of the new updated table
                    which we will write to.
        hitMapDict -- a dictionary with string keys for the sequences that received
                      BLAST hits and hence we want to obtain sequence details from
                      the blastFasta file. We don't care what its values are, since
                      they will be overwritten.
        representativesFasta -- a string indicating the location of the BINge
                                representatives FASTA file.
        blastFasta -- a string indicating the location of the FASTA file that was used
                      for BLAST; it's expected to RefSeq, UniRef, and/or reference
                      genome sequences (which won't be handled in any special way)
    '''
    nuisanceStrings = ["PREDICTED: ", "LOW QUALITY PROTEIN: "]
    refseqRegex = re.compile(r">.+?\s(.+?)\sn=\d+?\sTax=(.+?)\sTaxID=")
    unirefRegex = re.compile(r">.+?\s(.+?)\s\[(.+?)\]")
    
    # Parse the blastFasta file to get details for relevant sequences
    with open(blastFasta, "r") as fileIn:
        for line in fileIn:
            if line.startswith(">"):
                # Extract sequence ID from FASTA line
                seqID = line[1:].rstrip("\r\n ").split(" ")[0]
                seqID = seqID.split("_", maxsplit=1)[1] if seqID.startswith("UniRef") else seqID
                
                # Store value in dict
                if seqID in hitMapDict:
                    # See if this sequence has information we can obtain
                    refseqMatch = refseqRegex.match(line)
                    unirefMatch = unirefRegex.match(line)
                    
                    # Extract any possible information
                    if refseqMatch != None:
                        geneName, taxaName = refseqMatch.group(1), refseqMatch.group(2)
                    elif unirefMatch != None:
                        geneName, taxaName = unirefMatch.group(1), unirefMatch.group(2)
                    else:
                        geneName, taxaName = ".", "."
                    
                    # Remove any nuisance strings
                    for string in nuisanceStrings:
                        geneName = geneName.replace(string, "")
                    if geneName.strip(" ") == "":
                        "I don't think this is necessary; it's just a precaution"
                        geneName = "."
                    
                    # Store data
                    hitMapDict[seqID] = [geneName, taxaName]
    
    # Parse the original table file into a dictionary of lines
    "This lets us write an ordered output later"
    originalLines = {}
    with open(originalTable, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            
            # Handle header lines
            if line.startswith("#"):
                header = sl
            # Handle content lines
            else:
                originalLines[sl[0]] = sl
    
    # Piece things together into a final annotation table
    with open(representativesFasta, "r") as fileIn, open(newTable, "w") as fileOut:
        # Write updated header
        originalLength = len(header)
        
        accessionIndex = header.index("Target_accession")
        for text in ["Taxa_names", "Gene_names"]:
            header.insert(accessionIndex + 1, text)
        fileOut.write("{0}\n".format("\t".join(header)))
        
        # Write content lines
        for line in fileIn:
            if line.startswith(">"):
                clusterID = line[1:].rstrip("\r\n").split(" ")[0]
                
                # Get an annotation line for this sequence
                if clusterID in originalLines:
                    sl = originalLines[clusterID]
                else:
                    sl = [clusterID] + (["."]*(originalLength-1))
                
                # If we have BLAST hits, get their information
                if sl[accessionIndex] != ".":
                    hitList = sl[accessionIndex].replace(" ", "").replace("]", "").split("[")
                    seqInfo = [ hitMapDict[hit] for hit in hitList ]
                    
                    # Format it for insertion into the split line
                    geneNames = "".join([
                        seqInfo[i][0] if i == 0 and len(seqInfo) == 1
                        else f"{seqInfo[i][0]} " if i == 0
                        else f"[{seqInfo[i][0]}]"
                        for i in range(len(seqInfo))
                    ])
                    taxaNames = "".join([
                        seqInfo[i][1] if i == 0 and len(seqInfo) == 1
                        else f"{seqInfo[i][1]} " if i == 0
                        else f"[{seqInfo[i][1]}]"
                        for i in range(len(seqInfo))
                    ])
                # Otherwise, just make some dummy information
                else:
                    geneNames = "."
                    taxaNames = "."
                
                # Insert data into split line
                sl.insert(accessionIndex + 1, taxaNames)
                sl.insert(accessionIndex + 1, geneNames)
                
                # Write to file
                fileOut.write("{0}\n".format("\t".join(sl)))
