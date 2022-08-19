from bs4 import BeautifulSoup
from Bio.Entrez import efetch
from Bio.Entrez import elink
from Bio.Entrez import esearch
from Bio import Entrez
Entrez.email = 'zhouwanding@gmail.com'

## this searches papers that matches the term, Note that the term can also be a PMID

## https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&linkname=pubmed_pubmed_refs&id=24752654
## see documentation of BeautifulSoup4: https://www.crummy.com/software/BeautifulSoup/bs4/doc/
## Bio.Entrez: https://biopython.org/docs/1.76/api/Bio.Entrez.html
## eFetch example: https://www.biostars.org/p/493725/
## https://stackoverflow.com/questions/24146466/pubmed-id-to-author-list-citation-python
## https://www.ncbi.nlm.nih.gov/pmc/tools/cites-citedby/

def main(search_term):

    xml_data = esearch(db="pubmed", retmax=10, term=search_term, idtype="acc", sort="relevance").read()
    soup = BeautifulSoup(xml_data, features="lxml")
    
    pmids = ",".join([tag.string for tag in soup.idlist.findAll("id")])
    xml_data1 = efetch(db='pubmed', id=pmids, retmode='xml').read()
    soup1 = BeautifulSoup(xml_data1, features="lxml")
    for tag in soup1.findAll("pubmedarticle"):
        authors = tag.findAll("author")
        try:
            if len(authors) == 0:
                print("{%s}" % tag.pmid.string)
                continue

            first_au = authors[0]
            last_au = authors[-1]

            if tag.journal is None or tag.journal.year is None:
                year = "NA"
            else:
                year = tag.journal.year.string

            if tag.journal is None or tag.journal.isoabbreviation is None:
                journal = "NA"
            else:
                journal = tag.journal.isoabbreviation.string

            print("{%s} (%s %s %s) (%s) %s" % (
                tag.pmid.string,
                first_au.collectivename.text if first_au.lastname is None else first_au.lastname.text,
                year,
                journal,
                last_au.collectivename.text if last_au.lastname is None else last_au.lastname.text,
                tag.articletitle.string))
        except:
            print("{%s}" % tag.pmid.string)
    return

import sys
main(sys.argv[1])


