#!/usr/bin/env python

import argparse
import sys
import scipy.stats as stats
import numpy as np
import re

def opengz(fn):
    
    if fn.endswith('.gz'):
        import gzip
        fh = gzip.open(fn)
    else:
        fh = open(fn)

    return fh

class GOTerm():

    def __init__(self, termid):

        self.id = termid
        self.cnt = 0
        self.refcnt = 0
        self.children = []

    def __repr__(self):

        return '<%s: %s>' % (self.id, self.name)

def parse_go_terms(obofn):

    fh = open(obofn)
    for line in fh:
        if line.startswith('[Term]'):
            break

    keys = []
    terms = {}
    for line in fh:
        if line.strip() == '':
            continue
        if line.strip() == '[Term]':
            continue
        if line.strip() == '[Typedef]':
            break
        if line.startswith('id:'):
            termid = line.strip()[4:]
            assert termid.startswith('GO:')
            if termid in terms:
                ob = terms[termid]
            else:
                ob = GOTerm(termid)
                terms[termid] = ob
        else:
            k, v = line.strip().split(': ',1)

            if k == 'name':
                ob.name = v
                continue
            if k == 'def':
                ob.tdef = v
                continue
            if k == 'alt_id':
                terms[v] = ob
                continue

            if v.startswith('GO:'):
                _termid = v.split()[0]
                if _termid in terms:
                    _ob = terms[_termid]
                else:
                    _ob = GOTerm(_termid)
                    terms[_termid] = _ob
                v = _ob
                    
            if hasattr(ob, k):
                if isinstance(getattr(ob, k), list):
                    getattr(ob, k).append(v)
                else:
                    setattr(ob, k, [getattr(ob, k), v])
            else:
                setattr(ob, k, v)

    # for term in terms.itervalues():
    #     print
    #     print '============'
    #     print term.id
    #     print term.name
    #     print term.tdef
    #     if hasattr(term, 'is_a'):
    #         print 'is a', term.is_a

    return terms

def parse_association(assoc, terms):

    gene2go = {}
    for line in opengz(assoc):

        if line.startswith('!'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 3:
            continue
        assert(fields[4].startswith('GO:'))
        g = fields[2].upper()
        t = terms[fields[4]]
        if g not in gene2go:
            gene2go[g] = []
        gene2go[g].append(t)
        
    return gene2go

def inc_refcnt(g, gene2go):

    left = gene2go[g][:]
    taboo = set()
    while left:
        t = left.pop()
        if t in taboo:
            continue
        t.refcnt += 1
        taboo.add(t)
        if hasattr(t, 'is_a'):
            if isinstance(t.is_a, list):
                left.extend(t.is_a)
            else:
                left.append(t.is_a)

def inc_cnt(g, gene2go):

    left = gene2go[g][:]
    taboo = set()
    while left:
        t = left.pop()
        if t in taboo:
            continue
        t.cnt += 1
        taboo.add(t)
        if hasattr(t, 'is_a'):
            if isinstance(t.is_a, list):
                left.extend(t.is_a)
            else:
                left.append(t.is_a)

def mark_counts(genes, gene2go):

    i = 0
    for i, g in enumerate(genes):
        if g not in gene2go:
            sys.stderr.write('gene: %s not recognized, skip\n' % g)
            continue
        inc_cnt(g, gene2go)

def main_split_class(args):

    terms = parse_go_terms(args.obo)
    gene2go = parse_association(args.a, terms)

    roots = []
    for t in terms.itervalues():
        if hasattr(t, 'is_a'):
            if isinstance(t.is_a, list):
                for _t in t.is_a:
                    if t not in _t.children:
                        _t.children.append(t)
            else:
                _t = t.is_a
                if t not in _t.children:
                    _t.children.append(t)
        else:
            roots.append(t)

    for g in gene2go:
        inc_refcnt(g, gene2go)

    for t in terms.itervalues():
        if t.refcnt > args.t:
            print '%s\t%s\t%s\t%d' % (t.id, t.namespace, re.sub(' ', "_", t.name), t.refcnt)
            
    # tovisit = roots
    # visited = set()
    # cats = []
    # while tovisit:
    #     t = tovisit.pop()
    #     visited.add(t)
    #     if t.refcnt < args.t and t.refcnt > 0:
    #         if t not in cats:
    #             cats.append(t)
    #     else:
    #         # tovisit.extend(t.children)
    #         for _ in t.children:
    #             if _ not in visited:
    #                 tovisit.append(_)

    # for t in cats:
    #     print t, t.refcnt
    return

def main_enrich(args):

    terms = parse_go_terms(args.obo)
    gene2go = parse_association(args.a, terms)

    # background reference genes
    if args.reftable is None:
        for g in gene2go:
            inc_refcnt(g, gene2go)
        M = len(gene2go)
    else:
        M = 0
        for line in open(args.reftable):
            fields = line.strip().split('\t')
            g = fields[args.rc-1].upper()
            inc_refcnt(g, gene2go)
            M += 1
    sys.stderr.write('load %d reference genes\n' % M)

    # target genes
    genes = []
    for line in args.table:
        fields = line.strip().split('\t')
        gene = fields[args.c-1].upper()
        genes.append(gene)
    N = len(genes)
    mark_counts(genes, gene2go)

    res = []
    for t in terms.itervalues():
        if t.cnt > 0:
            x = t.cnt
            n = t.refcnt
            a = (-stats.hypergeom.logsf(x, M, n, N), float(x)/N/(float(n)/M), t, x, n)
            res.append(a)

    print 'pval\tenrichment_ratio\t#hit\t#target\t#category\t#total\tid\tname'
    for mlogp, ratio, t, x, n in sorted(res, reverse=True):
        if mlogp > args.t:
            print '%1.2f\t%1.2f\t%d\t%d\t%d\t%d\t%s\t%s' % (mlogp, ratio, x, N, n, M, t.id, t.name)

def add_go_args(parser):
    parser.add_argument('-reftable', help="reference gene table [if none, use all genes]", default=None)
    parser.add_argument('-obo', help="GO obo file (term definition)")
    parser.add_argument('-a', help="annotation association file for the organism")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='GO analysis')
    subparsers = parser.add_subparsers()

    s = subparsers.add_parser("split", help="split gene set based on go category")
    s.add_argument('-t', type=int, default=1000, help='threshold')
    add_go_args(s)
    s.set_defaults(func=main_split_class)

    s = subparsers.add_parser('enrich', help='find enriched GO terms, ex. wzgo enrich -obo ~/projects/pj-mm/2015-04-23-alu/go/go.obo -a gene_association.goa_human.gz pluripotential_genes | les')
    add_go_args(s)
    s.add_argument('table', help="target gene table", type = argparse.FileType('r'), default='-')
    s.add_argument('-t', type=int, default=10, help='loglikelihood threshold (default, 10)')
    s.add_argument('-c', type=int, default=1, help="columns of gene name, 1-based. (default 1)")
    s.set_defaults(func=main_enrich)

    args = parser.parse_args()
    args.func(args)
