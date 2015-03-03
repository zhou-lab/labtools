#include <ctype.h>
#include "wqueue.h"
#include "encode.h"
#include "sam.h"
#include "refseq.h"
#include "kstring.h"
#include "wvec.h"

typedef struct {
	int32_t tid;
  uint32_t beg, end;
} window_t;

DEFINE_WQUEUE(window, window_t)

typedef struct {
	int step;
	int n_threads;
	uint32_t min_base_qual;
	uint32_t max_retention;
	uint32_t min_read_len;
	uint8_t min_dist_end;
	uint8_t max_nm;
	uint8_t filter_secondary:1;
	uint8_t print_retention_only:1;
} conf_t;

typedef enum {BSS_RETENTION, BSS_CONVERSION, BSS_OTHER} bsstate_t;
typedef struct {
	uint8_t base:7;
	uint8_t bsstrand:1;
	uint8_t qual:7;
	uint8_t strand:1;
	uint16_t qpos;
	bsstate_t bsstate;
	uint8_t cnt_ret;
	uint16_t rlen;								/* read length */
} __attribute__((__packed__)) pileup_data_t;

DEFINE_VECTOR(pileup_data_v, pileup_data_t)

typedef struct {
	int n;												/* number of sites */
	pileup_data_v **data;
} pileup_t;

pileup_t *init_pileup(int n) {
	pileup_t *p = malloc(sizeof(pileup_t));
	p->n = n;
	p->data = calloc(n, sizeof(pileup_data_v*));
	return p;
}

void destroy_pileup(pileup_t *p) {
	int i;
	for (i=0; i<p->n; ++i) {
		if (p->data[i]) free_pileup_data_v(p->data[i]);
	}
	free(p->data);
	free(p);
}

DEFINE_WQUEUE(record, char*)

typedef struct {
	char *bam_fn;									/* on stack */
	char *ref_fn;									/* on stack */
  wqueue_t(window) *q;
	wqueue_t(record) *rq;
	conf_t *conf;
} result_t;

typedef struct {
	wqueue_t(record) *q;
	char *outfn;
} writer_conf_t;

void *write_func(void *data) {
	writer_conf_t *c = (writer_conf_t*) data;
	FILE *out = fopen(c->outfn, "w");
	while (1) {
		char *rec;
		wqueue_get(record, c->q, &rec);
		if(!rec) break;
		fputs(rec, out);
		free(rec);									/* rec is alloc-ed */
	}
	fclose(out);
	return 0;
}

uint32_t plp_cnt_retention(pileup_data_v *plp_data) {
	uint32_t i, cnt = 0;
	for (i=0; i<plp_data->size; ++i) {
		pileup_data_t *d = ref_pileup_data_v(plp_data, i);
		if (d->bsstate == BSS_RETENTION) cnt++;
	}
	return cnt;
}

char *plp_format(refseq_t *rs, char *chrm, uint32_t rpos, pileup_data_v *dv, conf_t *conf) {
	uint32_t i;
	kstring_t s;
	s.l = s.m = 0; s.s = 0;
	char rb = toupper(getbase_refseq(rs, rpos));

	ksprintf(&s, "%s\t%u\t%u\t%c", chrm, rpos, rpos, rb);

	/* context */
	char trinuc[3];
	if (rpos == 1) {
		subseq_refseq2(rs, 1, trinuc, 2);
		trinuc[0] = 'N';
	} else if (rpos == (unsigned) rs->seqlen) {
		subseq_refseq2(rs, rpos-1, trinuc, 2);
		trinuc[2] = 'N';
	} else {
		subseq_refseq2(rs, rpos-1, trinuc, 3);
	}
	if (rb == 'C') {
		ksprintf(&s, "\t%.3s", trinuc);
	} else if (rb == 'G') {
		char trinuc_r[3];
		_nt256char_rev(trinuc_r, trinuc, 3);
		ksprintf(&s, "\t%.3s", trinuc_r);
	} else {
		fprintf(stderr, "[%s:%d] This is a bug.\n", __func__, __LINE__);
		exit(1);
	}

	/* coverage */
	ksprintf(&s, "\t%u", dv->size);

	/* count retention and conversion */
	int cnt_retention=0, cnt_conversion=0;
	for (i=0; i<dv->size; ++i) {
		pileup_data_t *d = ref_pileup_data_v(dv,i);
		if (d->qual < conf->min_base_qual) continue;
		if (d->qpos < conf->min_dist_end ||
				d->rlen < d->qpos + conf->min_dist_end) continue;
		if (d->bsstate == BSS_RETENTION) cnt_retention++;
		if (d->bsstate == BSS_CONVERSION) cnt_conversion++;
	}
	ksprintf(&s, "\t%d\t%d", cnt_retention, cnt_conversion);

	kputc('\t', &s);
	uint8_t nf = 0;
	for (i=0; i<dv->size; ++i) {
		pileup_data_t *d = ref_pileup_data_v(dv,i);
		if (conf->print_retention_only && d->bsstate != BSS_RETENTION) continue;
		if (nf) kputc(',', &s);
		else nf = 1;
		kputuw(d->cnt_ret, &s);
	}
	
	kputc('\t', &s);
	for (i=0; i<dv->size; ++i) {
		pileup_data_t *d = ref_pileup_data_v(dv,i);
		if (conf->print_retention_only && d->bsstate != BSS_RETENTION) continue;
		kputc(d->strand?'-':'+', &s);
	}
	
	kputc('\t', &s);
	nf = 0;
	for (i=0; i<dv->size; ++i) {
		pileup_data_t *d = ref_pileup_data_v(dv,i);
		if (conf->print_retention_only && d->bsstate != BSS_RETENTION) continue;
		if (nf) kputc(',', &s);
		else nf = 1;
		kputuw(d->qpos, &s);
	}

	kputc('\t', &s);
	for (i=0; i<dv->size; ++i) {
		pileup_data_t *d = ref_pileup_data_v(dv,i);
		if (conf->print_retention_only && d->bsstate != BSS_RETENTION) continue;
		kputc(d->qual+33, &s);
	}

	kputc('\t', &s);
	for (i=0; i<dv->size; ++i) {
		pileup_data_t *d = ref_pileup_data_v(dv,i);
		if (conf->print_retention_only && d->bsstate != BSS_RETENTION) continue;
		kputc(d->bsstrand?'-':'+', &s);
	}

	kputc('\n', &s);
	
	return s.s;
}

#define bscall(b, pos) bam_nt16_rev_table[bam1_seqi(bam1_seq(b), pos)]

/* return -1 if abnormal (missing bsstrand) */
int cnt_retention(refseq_t *rs, bam1_t *b) {
	int cnt = 0;
	uint8_t *bsstrand = bam_aux_get(b, "ZS");
	if (!bsstrand) return -1;
	bsstrand++;

	bam1_core_t *c = &b->core;
	uint32_t rpos = c->pos+1, qpos = 0;
	uint32_t op, oplen;
	char rb, qb;
	int i; unsigned j;
	for (i=0; i<c->n_cigar; ++i) {
		op = bam_cigar_op(bam1_cigar(b)[i]);
		oplen = bam_cigar_oplen(bam1_cigar(b)[i]);
		switch(op) {
		case BAM_CMATCH:
			for (j=0; j<oplen; ++j) {
				rb = toupper(getbase_refseq(rs, rpos+j));
				qb = bscall(b, qpos+j);
				if (bsstrand[0] == '+' && rb == 'C' && qb == 'C') cnt++;
				else if (bsstrand[0] == '-' && rb == 'G' && qb == 'G') cnt++;
				else continue;
			}
			rpos += oplen;
			qpos += oplen;
			break;
		case BAM_CINS:
			qpos += oplen;
			break;
		case BAM_CDEL:
			rpos += oplen;
			break;
		case BAM_CSOFT_CLIP:
			qpos += oplen;
			break;
		default:
			fprintf(stderr, "Unknown cigar, %u\n", op);
			abort();
		}
	}

	return cnt;
}

void *process_func(void *data) {

  result_t *res = (result_t*) data;
	conf_t *conf = (conf_t*) res->conf;
	samfile_t *in = samopen(res->bam_fn, "rb", 0);
	bam_index_t *idx = bam_index_load(res->bam_fn);
	refseq_t *rs = init_refseq(res->ref_fn, 1000, 1000);

	window_t w;
	int i; uint32_t j;
  while (1) {

    wqueue_get(window, res->q, &w);
    if (w.tid == -1) break;

		pileup_t *plp = init_pileup(w.end - w.beg);
		
		char *chrm = in->header->target_name[w.tid];
		fetch_refseq(rs, chrm, w.beg>100?w.beg-100:1, w.end+100);
		bam_iter_t iter = bam_iter_query(idx, w.tid, w.beg, w.end);
		bam1_t *b = bam_init1();
		int ret;
		char qb, rb;
		while ((ret = bam_iter_read(in->x.bam, iter, b))>0) {

			uint8_t *bsstrand = bam_aux_get(b, "ZS");
			if (!bsstrand) continue;
			bsstrand++;

			bam1_core_t *c = &b->core;

			if (c->l_qseq < 0 || (unsigned) c->l_qseq < conf->min_read_len) continue;
			
			uint32_t rpos = c->pos+1, qpos = 0;
			if (conf->filter_secondary && c->flag & BAM_FSECONDARY) continue;
			uint8_t *nm = bam_aux_get(b, "NM");
			if (nm && bam_aux2i(nm)>conf->max_nm) continue;

			int cnt_ret = cnt_retention(rs, b);
			if (cnt_ret < 0 || (unsigned) cnt_ret > conf->max_retention) continue;

			rpos = c->pos+1; qpos = 0;
			for (i=0; i<c->n_cigar; ++i) {
				uint32_t op = bam_cigar_op(bam1_cigar(b)[i]);
				uint32_t oplen = bam_cigar_oplen(bam1_cigar(b)[i]);
				switch(op) {
				case BAM_CMATCH:
					for (j=0; j<oplen; ++j) {
						if (rpos+j<w.beg || rpos+j>=w.end) continue; /* include begin but not end */
						rb = toupper(getbase_refseq(rs, rpos+j));
						if (rb != 'C' && rb != 'G') continue;
						qb = bscall(b, qpos+j);
						pileup_data_v **plp_data_vec = plp->data+rpos+j-w.beg;
						if (!*plp_data_vec) *plp_data_vec = init_pileup_data_v(2);
						pileup_data_t *plp_data = next_ref_pileup_data_v(*plp_data_vec);
						plp_data->base = qb;
						plp_data->qual = bam1_qual(b)[qpos+j];
						if (bsstrand[0] == '-') plp_data->bsstrand = 1;
						else if (bsstrand[0] == '+') plp_data->bsstrand = 0;
						else continue;
						if (rb == 'C' && !plp_data->bsstrand) {
							if (qb == 'C') plp_data->bsstate = BSS_RETENTION;
							else if (qb == 'T') plp_data->bsstate = BSS_CONVERSION;
							else plp_data->bsstate = BSS_OTHER;
						} else if (rb == 'G' && plp_data->bsstrand) {
							if (qb == 'G') plp_data->bsstate = BSS_RETENTION;
							else if (qb == 'A') plp_data->bsstate = BSS_CONVERSION;
							else plp_data->bsstate = BSS_OTHER;
						} else plp_data->bsstate = BSS_OTHER;
						plp_data->cnt_ret = (unsigned) cnt_ret;
						plp_data->strand = (c->flag&BAM_FREVERSE)?1:0;
						plp_data->qpos = qpos+j;
						plp_data->rlen = c->l_qseq;
					}
					rpos += oplen;
					qpos += oplen;
					break;
				case BAM_CINS:
					qpos += oplen;
					break;
				case BAM_CDEL:
					rpos += oplen;
					break;
				case BAM_CSOFT_CLIP:
					qpos += oplen;
					break;
				default:
					fprintf(stderr, "Unknown cigar, %u\n", op);
					abort();
				}
			}
		}

		for (j=w.beg; j<w.end; ++j) {
			rb = getbase_refseq(rs, j);
			pileup_data_v *plp_data = plp->data[j-w.beg];
			if (plp_data && plp_cnt_retention(plp_data)>0) {
				wqueue_put2(record, res->rq, plp_format(rs, chrm, j, plp_data, res->conf));
			}
		}

		destroy_pileup(plp);
		bam_destroy1(b);
		bam_iter_destroy(iter);
  }
	free_refseq(rs);
	samclose(in);
	bam_index_destroy(idx);
	return 0;
}


int pileup_cytosine_main(char *reffn, char *infn, char *outfn, char *reg, conf_t *conf) {

  wqueue_t(window) *wq = wqueue_init(window, 100000);
  pthread_t *processors = calloc(conf->n_threads, sizeof(pthread_t));
  result_t *results = calloc(conf->n_threads, sizeof(result_t));
	int i;
	samfile_t *in = samopen(infn, "rb", 0);

	pthread_t writer;
	writer_conf_t writer_conf = {
		.q = wqueue_init(record, 100000),
		.outfn = outfn,
	};
	pthread_create(&writer, NULL, write_func, &writer_conf);
  for (i=0; i<conf->n_threads; ++i) {
		results[i].q = wq;
		results[i].rq = writer_conf.q;
		results[i].ref_fn = reffn;
		results[i].bam_fn = infn;
		results[i].conf = conf;
    pthread_create(&processors[i], NULL, process_func, &results[i]);
  }

	window_t w; memset(&w, 0, sizeof(window_t));
	uint32_t wbeg;
	if (reg) {
		int tid;
		uint32_t beg, end;
		bam_parse_region(in->header, reg, &tid, (int*) &beg, (int*) &end);
		/* chromosome are assumed to be less than 2**29 */
		beg++; end++;
		if (beg<=0) beg = 1;
		if (end>in->header->target_len[tid]) end = in->header->target_len[tid];
		for (wbeg = beg; wbeg < end; wbeg += conf->step) {
			w.tid = tid;
			w.beg = wbeg;
			w.end = wbeg + conf->step;
			if (w.end > end) w.end = end;
			wqueue_put(window, wq, &w);
		}
	} else {
		for (i=0; i<in->header->n_targets; ++i) {
			uint32_t target_len = in->header->target_len[i];
			for (wbeg = 1; wbeg < target_len; wbeg += conf->step) {
				w.tid = i;
				w.beg = wbeg;
				w.end = wbeg+conf->step;
				if (w.end > target_len) w.end = target_len;
				wqueue_put(window, wq, &w);
			}
		}
	}
  for (i=0; i<conf->n_threads; ++i) {
    w.tid = -1;
    wqueue_put(window, wq, &w);
  }

  for (i=0; i<conf->n_threads; ++i) {
    pthread_join(processors[i], NULL);
	}

	wqueue_put2(record, writer_conf.q, NULL);
	pthread_join(writer, NULL);
	wqueue_destroy(record, writer_conf.q);

	free(results);
	free(processors);
  wqueue_destroy(window, wq);
	samclose(in);

	return 0;
}

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: pileup_cytosine [options] -r [ref.fa] -i [in.bam] -o [out.pileup] -g [chr1:123-234]\n");
	fprintf(stderr, "output format: chrm, pos, pos, refbase, context, coverage, filtered_retention, filtered_conversion, num_retention_in_reads, strands, position_on_reads, quality, bsstrand\n");
  fprintf(stderr, "Input options:\n");
	fprintf(stderr, "     -i        input bam.\n");
  fprintf(stderr, "     -r        reference in fasta.\n");
	fprintf(stderr, "     -g        region (optional, if not specified the whole bam will be processed).\n");
	fprintf(stderr, "     -o        pileup output file\n");
	fprintf(stderr, "     -s        step of window dispatching [100000].\n");
	fprintf(stderr, "     -q        number of threads [3] recommend 20.\n");
	fprintf(stderr, "     -b        min base quality [10].\n");
	fprintf(stderr, "     -t        max retention in a read [5].\n");
	fprintf(stderr, "     -l        minimum read length [60].\n");
	fprintf(stderr, "     -e        minimum distance to end of a read [3].\n");
	fprintf(stderr, "     -c        turn off filtering secondary mapping.\n");
	fprintf(stderr, "     -n        maximum nm tag [2].\n");
	fprintf(stderr, "     -a        print all, not just retention reads.\n");
  fprintf(stderr, "     -h        this help.\n");
  fprintf(stderr, "\n");
  return 1;
}

int main(int argc, char *argv[]) {

	int c;
	char *reffn = 0;
	char *reg = 0;
	char *infn = 0;
	char *outfn = 0;
	conf_t conf = {
		.step = 100000,
		.n_threads = 3,
		.min_base_qual = 10,
		.max_retention = 5,
		.min_read_len = 60,
		.filter_secondary = 1,
		.print_retention_only = 1,
		.min_dist_end = 3,
		.max_nm = 2,
	};


	if (argc<2) return usage();
	while ((c=getopt(argc, argv, "i:o:r:g:q:e:b:t:l:cah"))>=0) {
		switch (c) {
		case 'i': infn = optarg; break;
		case 'r': reffn = optarg; break;
		case 'g': reg = optarg; break;
		case 'o': outfn = optarg; break;
		case 's': conf.step = atoi(optarg); break;
		case 'q': conf.n_threads = atoi(optarg); break;
		case 'b': conf.min_base_qual = atoi(optarg); break;
		case 't': conf.max_retention = atoi(optarg); break;
		case 'l': conf.min_read_len = atoi(optarg); break;
		case 'e': conf.min_dist_end = atoi(optarg); break;
		case 'c': conf.filter_secondary = 0; break;
		case 'a': conf.print_retention_only = 0; break;
		case 'n': conf.max_nm = atoi(optarg); break;
		case 'h': return usage();
    default:
      fprintf(stderr, "[%s:%d] Unrecognized command: %c.\n", __func__, __LINE__, c);
      exit(1);
      break;
    }
	}

	if (!infn || !reffn || !outfn) usage();

	pileup_cytosine_main(reffn, infn, outfn, reg, &conf);
}
