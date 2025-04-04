import re, requests

def fix_obsoletions(goList, goObo, queriedGOs):
    '''
    Helper function to receive a list of GO terms, alongside the parsed go.obo
    and a dictionary potentially containing previous API queries. Using this, it
    will check if the GOs in the GO list are obsolete and, if so, attempt to fix
    them via API query.
    
    Parameters:
        goList -- a list of string GO terms e.g., ['GO:0033644', 'GO:0016020'].
        goObo -- a obo_parser.GODag object which has parsed a go.obo file.
        queriedGOs -- a dictionary with previous API query hits (if any) that
                      will be modified by this function; has structure like:
                      {
                          'goID1': ['replacementID1'],
                          'goID2': ['replacementID1', 'replacementID2'],
                          'goID3': None, # if no replacement is available
                          ...
                      }
    '''
    fixedGOs = []
    for goTerm in goList:
        
        # If we have a hit, everything's good!
        if goTerm in goObo:
            fixedGOs.append(goTerm)
        
        # Otherwise, if we found an obsoletion, handle it
        else:
            # Perform an API query if we haven't seen this GO yet
            if not goTerm in queriedGOs:
                replacements = query_go_api(goTerm)
                queriedGOs[goTerm] = replacements
            
            # Get the API query result
            apiGOs = queriedGOs[goTerm]
            
            # Store result if we were successful
            if apiGOs != None:
                fixedGOs += apiGOs
    return fixedGOs

def query_go_api(goTerm):
    '''
    Parameters:
        goTerm -- a GO term to query for replacement to its obsoletion
    Returns:
        replacementList -- a list containing GO terms which replace this one, OR
                           None if no replacements were found
    '''
    # Query the API
    apiUrl = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO%3A{0}/history".format(
        goTerm.split(":")[1]
    )
    response = requests.get(apiUrl)
    json = response.json()
    
    # Handle errors
    if "messages" in json:
        message = json["messages"][0]
        print(f"# WARNING: Querying the GO term '{goTerm}' resulted in an API error")
        print(f"# > The return message was: '{message}'")
        print("# > No result will be available for this GO term")
        return None
    
    # Handle null hits
    if json["results"] == []:
        print(f"# WARNING: Querying the GO term '{goTerm}' returned no results")
        print("# > No result will be available for this GO term")
        return None
    
    # Parse results out of the JSON response
    resultsDict = json["results"][0]
    
    # Handle situations where the GO isn't actually obsolete
    if resultsDict["isObsolete"] != True:
        print(f"# WARNING: The GO term '{goTerm}' does not appear to be obsolete")
        print("# > However, it is not found in the go.obo file?")
        print("# > Ancestor terms will be unavailable for this GO term")
        return [goTerm]
    
    # Process results and try to find replacements or considers
    else:
        replacedBys = set()
        considers = set()
        for historyDict in resultsDict["history"]:
            historyText = historyDict["text"]
            if "consider GO:" in historyText:
                considerMatch = re.search(r"consider (GO:\d{7})", historyText)
                if considerMatch != None:
                    considers.add(considerMatch.group(1))
            if "replaced_by GO:" in historyText:
                replacedByMatch = re.search(r"replaced_by (GO:\d{7})", historyText)
                if replacedByMatch != None:
                    replacedBys.add(replacedByMatch.group(1))
    
    # Return the best value, or None if no replacements were possible
    if len(replacedBys) > 0:
        return list(replacedBys)
    elif len(considers) > 0:
        return list(considers)
    else:
        return None
