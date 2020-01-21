#pragma once
#define KTHREAD_STORAGE static inline
#include "./kthread.h"
#include "./kthread.c"

namespace rg {
struct ForPool {
    void *fp_;
    ForPool(int nthreads): fp_(kt_forpool_init(nthreads)) {}
    ForPool(const ForPool &) = delete;
    ForPool(ForPool &&) = delete;
    void forpool(void (*func)(void*,long,int), void *data, long n) {
        kt_forpool(fp_, func, data, n);
    }
    ~ForPool() {
        kt_forpool_destroy(fp_);
    }
};
}
