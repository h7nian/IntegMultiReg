#ifndef IMR_CV_VALUE_CACHE_H
#define IMR_CV_VALUE_CACHE_H

#include <stddef.h>

typedef struct {
    void *data;
    size_t bytes;
    int previous;
    int next;
} postfit_value_cache_entry;

typedef struct {
    postfit_value_cache_entry *entries;
    size_t n_entries;
    size_t budget_bytes;
    size_t external_bytes;
    size_t entry_bytes;
    size_t used_bytes;
    int first;
    int last;
    void *(*allocate)(size_t);
    void (*release)(void *);
} postfit_value_cache;

/* Budget counts externally owned key/index bytes plus this store's requested
 * entry and payload bytes. Allocation failure disables caching or yields a miss.
 * Returned payload pointers remain valid only until the next put or destroy.
 * Keys are immutable model identities within one fold/subgroup. */
int postfit_value_cache_init(postfit_value_cache *cache, size_t n_entries,
                            size_t budget_bytes, size_t external_bytes,
                            void *(*allocate)(size_t), void (*release)(void *));
void *postfit_value_cache_get(const postfit_value_cache *cache, size_t key);
void *postfit_value_cache_put(postfit_value_cache *cache, size_t key, size_t bytes);
void postfit_value_cache_destroy(postfit_value_cache *cache);

#endif
