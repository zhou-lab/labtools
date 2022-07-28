#This script can take in list of PMIDs and find (most of the time) impact factor of journal and # of times study was cited

import re
import json
import bs4
import requests


#input tsv that has a column of PMIDs
file = "/Users/davidgoldberg/Desktop/20220519_EWASmaster.tsv"


#loop through lines of file and store PMID #'s in list. For NA's, keep track of titles
pmids = []
missing_pmids = []
#import file containing PMIDs
with open(file, 'r') as tsv:
	next(tsv)
	for line in tsv:
		pmid = line.split("\t")[4].strip()
		if pmid in pmids:
			next
		elif pmid == "NA":
			missing_pmids.append(line.split("\t")[7].strip())
		else:
			pmids.append(pmid)


#Loop through PMIDs, request citation info from ncbi API, make dictionary
# with PMID IDs as keys, study title and journal as values. Also track studies where no data pulled
missing_data = []
pmid_dict = {}
for pmid in pmids[901:1064]:   #note this line run repeatedly in blocks to not make too many requests at once. I missed two entries by not getting the overlapping index.. some issues at end of 701:901
	try:
		print(pmid)
		url = "https://api.ncbi.nlm.nih.gov/lit/ctxp/v1/pubmed/?format=citation&contenttype=json&id=" + pmid
		response = requests.get(url)
		ama = response.json()["ama"]
		info = ama["orig"]
		pmid_dict[pmid] = info.split(".")[1:3]
	except:
		print("error")
		missing_data.append(pmid)


#Now get impact factor journals and # of citations for study
for pmid in pmid_dict.keys():
	try:
		print(pmid_dict[pmid][1])
		impact_f_url = "https://www.google.com/search?q=" + pmid_dict[pmid][1].strip() + "+impact+factor"
		response_if = requests.get(impact_f_url)
		soup = bs4.BeautifulSoup(response_if.text, "html.parser")
		imp_fac = soup.find_all("div", class_="BNeawe iBp4i AP7Wnd")[0].get_text()
		pmid_dict[pmid].append(imp_fac)
	except:
		print("error (IF)")
		pmid_dict[pmid].append("IF not found")
	try:
		cits_url = "https://www.google.com/search?q=" + pmid_dict[pmid][0].strip().replace(" ", "+")
		cits_response = requests.get(cits_url)
		cits_soup = bs4.BeautifulSoup(cits_response.text, "html.parser")
		cit_num = cits_soup.find_all(text=re.compile("Cited by"))[0]
		print(cit_num)
		pmid_dict[pmid].append(cit_num)
	except:
		print("error (citations)")
		pmid_dict[pmid].append("citations not found")

    
#Write results to file 
with open('/Users/davidgoldberg/Desktop/pmid_info.txt', 'a') as f:
    for key, value in pmid_dict.items():
        f.write('%s:%s\n' % (key, value))
f.close()


#file moved to /Users/davidgoldberg/Dropbox/Family Room/Lab Gallery/goldbergd/temp/pmid_info.txt 


