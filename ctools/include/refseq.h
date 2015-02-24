#ifndef _WZ_REFSEQ_H_
#define _WZ_REFSEQ_H_

#include "faidx.h"

typedef struct {
	faidx_t *fai;
	char *chrm;
	uint32_t beg;
	uint32_t end;
	char *seq;
	uint32_t flank1;
	uint32_t flank2;
} refseq_t;

static inline refseq_t* init_refseq(char *ref_fn, uint32_t flank1, uint32_t flank2) {
	
	refseq_t *rs = calloc(1, sizeof(refseq_t));
  rs->fai = fai_load(ref_fn);
	rs->flank1 = flank1;
	rs->flank2 = flank2;
	return rs;
}


/* beg and end are 1-based */
static inline void fetch_refseq(refseq_t *rs, char *chrm, uint32_t beg, uint32_t end) {

	if (rs->chrm != 0
			&& strcmp(chrm, rs->chrm) == 0
			&& rs->beg <= beg
			&& rs->end >= end) return;
	else {

		/* get sequence length */
		int seqlen = faidx_seqlen(rs->fai, chrm);
		if (seqlen < 0) {
			fprintf(stderr, "[%s:%d] Error, cannot retrieve reference %s:%u-%u.\n",
							__func__, __LINE__, chrm, rs->beg, rs->end);
			exit(1);
		}

		/* beg and end */
		if (rs->flank1 > beg + 1) rs->beg = 1;
		else rs->beg = beg - rs->flank1;
		if (end + rs->flank2 > (unsigned) seqlen) rs->end = seqlen;
		else rs->end = end + rs->flank2;

		rs->chrm = realloc(rs->chrm, strlen(chrm)+1);
		strcpy(rs->chrm, chrm);
		if (rs->seq) free(rs->seq);
		int l;
		rs->seq = faidx_fetch_seq(rs->fai, rs->chrm, rs->beg-1, rs->end-1, &l);
		if ((unsigned) l != rs->end-rs->beg+1){
			fprintf(stderr, "[%s:%d] Error, cannot retrieve reference %s:%u-%u.\n",
							__func__, __LINE__, chrm, rs->beg, rs->end);
			exit(1);
		}
	}
}

static inline void free_refseq(refseq_t *rs) {

	if (rs->seq) free(rs->seq);
	if (rs->chrm) free(rs->chrm);
	fai_destroy(rs->fai);
  free(rs);

}

/* rpos is 1-based */
static inline char getbase_refseq(refseq_t *rs, uint32_t rpos) {
	return rs->seq[rpos-rs->beg];
}

/* rpos is 1-based */
static inline char *subseq_refseq(refseq_t *rs, uint32_t rpos) {
	return rs->seq+rpos-rs->beg;
}
#endif /* _WZ_REFSEQ_H_ */
