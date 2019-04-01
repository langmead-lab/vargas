#ifndef KTHREAD_H
#define KTHREAD_H

#ifndef KTHREAD_STORAGE
#define KTHREAD_STORAGE
#endif

#ifdef __cplusplus
extern "C" {
#endif

KTHREAD_STORAGE void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
KTHREAD_STORAGE void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

KTHREAD_STORAGE void *kt_forpool_init(int n_threads);
KTHREAD_STORAGE void kt_forpool_destroy(void *_fp);
KTHREAD_STORAGE void kt_forpool(void *_fp, void (*func)(void*,long,int), void *data, long n);

#ifdef __cplusplus
}
#endif

#endif
