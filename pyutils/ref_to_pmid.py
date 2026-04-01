"""
Script: ref_to_pmid.py
Description: Automates the process of finding PMIDs for a list of bibliographic references 
             and generates formatted records for literatureAnnotation.org.

Usage:
    python ref_to_pmid.py <references_file>

Arguments:
    <references_file> : A text file containing one reference per line.

Outputs:
    pmid_mapping.json     : A JSON file mapping original reference strings to found PMIDs.
    formatted_records.txt : A text file containing literatureAnnotation-compatible records 
                            (generated using entrez_esearch.py).

Example:
    1. Create a file 'refs.txt' with:
       Belsky, D. W. et al. 2022. "DunedinPACE..."
    2. Run:
       python ref_to_pmid.py refs.txt
"""

import urllib.request
import urllib.parse
import json
import time
import sys
import os

def search_pubmed(query):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
    params = {
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": 1
    }
    url = base_url + urllib.parse.urlencode(params)
    try:
        with urllib.request.urlopen(url) as response:
            data = json.loads(response.read().decode())
            idlist = data.get("esearchresult", {}).get("idlist", [])
            return idlist[0] if idlist else None
    except Exception:
        return None

def get_record(pmid):
    # This calls the existing lab tool to get the formatted record
    cmd = f"python -Wignore ~/repo/labtools/pyutils/entrez_esearch.py {pmid}"
    try:
        output = os.popen(cmd).read().strip()
        return output
    except Exception:
        return f"{{{pmid}}} (Error fetching record)"

def main():
    if len(sys.argv) < 2:
        print("Usage: python ref_to_pmid.py <references_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    if not os.path.exists(input_file):
        print(f"File not found: {input_file}")
        sys.exit(1)

    with open(input_file, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    results = []
    mapping = {}

    for line in lines:
        # Remove leading numbers/dots if present (e.g., "1. Author...")
        clean_ref = re.sub(r"^\d+[\.\)]\s*", "", line)
        
        # Try searching by a few words from the middle/title to avoid issues with formatting
        # Often the title starts after the year. 
        # For simplicity, we search the whole cleaned string first.
        pmid = search_pubmed(clean_ref)
        
        if not pmid:
            # Fallback: try first 100 characters
            pmid = search_pubmed(clean_ref[:100])

        if pmid:
            mapping[line] = pmid
            record = get_record(pmid)
            print(f"Found: {pmid} for reference starting with: {line[:50]}...")
            results.append(record)
        else:
            print(f"NOT FOUND: {line[:50]}...")
            results.append(f"{{UNFOUND}} {line}")
        
        time.sleep(0.4) # Respect NCBI rate limits

    # Output mapping to a JSON for programmatic use
    with open("pmid_mapping.json", "w") as f:
        json.dump(mapping, f, indent=2)

    # Output formatted records for literatureAnnotation.org
    with open("formatted_records.txt", "w") as f:
        for res in results:
            f.write(res + "\n")

    print("\nProcessing complete.")
    print("PMID mapping saved to pmid_mapping.json")
    print("Formatted records saved to formatted_records.txt")

if __name__ == "__main__":
    import re
    main()
