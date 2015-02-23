#include "wqueue.h"

DEFINE_WQUEUE(window, window_t);
char trinuc[16][3] = {"ACA", "ACC", ""};

typedef struct {
  uint16_t tid;
  uint32_t beg, end;
} window_t;

typedef struct {
  int n_targets;
  long **total;
  long **retained;
  wqueue_t(window) *q;
} result_t;

result_t *init_result(int n_targets, wqueue_t(window) *q) {
  result_t *result = (result_t*) malloc(sizeof(result_t));
  result->n_targets = n_targets;
  result->total = (long**) malloc(n_targets*sizeof(long*));
  result->retained = (long**) malloc(n_targets*sizeof(long*));
  for (i=0; i<n_targets; ++i) {
    result->total[i] = (long*) calloc(16, sizeof(long));
    result->retained[i] = (long*) calloc(16, sizeof(long));
  }
  result->q = q;
  return result;
}

void destroy_result(result_t *result) {
  for (i=0; i<result->n_targets; ++i) {
    free(result->total[i]);
    free(result->retained[i]);
  }
  free(result->total);
  free(result->retained);
  free(result);
}

/* process windows */
void *process_func(void *data) {

  result_t *res = (result_t*) data;
  while (1) {
    wqueue_get(window, res->q, &w);
    if (w.chrm == -1) break;
    bam_iter_query(window_t, w.chrm, w.beg, w.end);
    res->total[trinuc2inds()]++;
    res->retained[trinuc2inds()]++
  }
}

void merge_results() {

}

int main() {

  wqueue_t(window) *wq = wqueue_init(window, 100000);
  pthread_t *processor = (pthread_t*) calloc(n_threads, sizeof(pthread_t));
  result_t **results = (result_t**) malloc(n_threads*sizeof(result_t*));
  for (i=0; i<n_threads; ++i) {
    results[i] = init_result(n_targets, wq);
    pthread_create(&processor[i], NULL, process_func, results[i]);
  }

  for (i=0; i<in->n_targets; ++i) {
    for (wbeg = 0; wbeg < end; wbeg += step) {
      w.chrm = targets[i];
      w.beg = wbeg;
      w.end = wbeg+step;
      if (w.end > end) w.end = end;
    }
    wqueue_put(window, wq, &w);
  }
  for (i=0; i<n_threads; ++i) {
    w.chrm = -1;
    wqueue.put(window, wq, &w);
  }

  for (i=0; i<n_threads; ++i)
    pthread_join(&processor[i], NULL);

  result_t final_result = init_result(n_targets, NULL);
  for (i=0; i<n_threads; ++i) {
    merge_results(&final_result, &results[i]);
    destroy_result(results[i]);
  }
  free(results);

  kstring_t s;
  s.l = s.m = 0; s.s = 0;
  kputs("trinuc", &s);
  for (j=0; j<16; ++j) {
    kputc('\t', &s);
    kputsn(trinucs[j], 3, &s);
  }
  puts(s.s);
  for (i=0; i<n_targets; ++i) {
    s.l = s.m = 0; free(s.s);
    kputs(targets[i], &s);
    for (j=0; j<16; ++j) ksprintf(&s, "\t%d", final_result->total[i][j]);
    for (j=0; j<16; ++j) ksprintf(&s, "\t%d", final_result->retained[i][j]);
    puts(s.s);
  }

  wqueue_destroy(window, wq);
}
