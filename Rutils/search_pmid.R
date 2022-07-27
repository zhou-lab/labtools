## source: https://stackoverflow.com/questions/69156614/fetching-multiple-pubmed-abstract-using-r-httr

library(XML)
library(httr)
library(glue)
library(dplyr)
####
####



query = 'asthma[mesh]+AND+leukotrienes[mesh]+AND+2009[pdat]'


reqq = glue ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={query}')


reqq = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=science[journal]+AND+breast+cancer+AND+2008[pdat]&usehistory=y"

op = GET(reqq)

content(op)


df_op <- op %>% xml2::read_xml() %>% xml2::as_list()

pmids <- df_op$eSearchResult$IdList %>% unlist(use.names = FALSE)

reqq1 = glue("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=pmids&rettype=abstract&retmode=xml")

op1 = GET(reqq1)

content(op1)
