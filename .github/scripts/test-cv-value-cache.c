#include "cv_value_cache.h"
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static void *live[32];
static size_t n_live, calls, fail_at;
static void *test_allocate(size_t bytes)
{
    if (++calls == fail_at) return NULL;
    void *pointer = malloc(bytes);
    assert(pointer && n_live < 32);
    live[n_live++] = pointer;
    return pointer;
}
static void test_release(void *pointer)
{
    size_t index = 0;
    while (index < n_live && live[index] != pointer) ++index;
    assert(index < n_live);
    live[index] = live[--n_live];
    free(pointer);
}
static void check_accounting(const postfit_value_cache *cache)
{
    size_t bytes = cache->external_bytes + cache->entry_bytes;
    if (cache->entries)
        for (size_t key = 0; key < cache->n_entries; ++key)
            bytes += cache->entries[key].bytes;
    assert(bytes == cache->used_bytes && bytes <= cache->budget_bytes);
}
int main(void)
{
    const size_t budget = 32 + 8 * sizeof(postfit_value_cache_entry) + 32;
    postfit_value_cache cache;
    assert(postfit_value_cache_init(&cache, 8, budget, 32, test_allocate, test_release));
    for (size_t key = 0; key < 8; ++key) {
        assert(postfit_value_cache_put(&cache, key, 16));
        check_accounting(&cache);
    }
    assert(!postfit_value_cache_get(&cache, 0));
    assert(postfit_value_cache_get(&cache, 6));
    assert(postfit_value_cache_put(&cache, 0, 16));
    assert(!postfit_value_cache_get(&cache, 6)); /* Hits do not refresh FIFO. */
    assert(postfit_value_cache_put(&cache, 7, 32));
    assert(!postfit_value_cache_get(&cache, 0));
    assert(!postfit_value_cache_put(&cache, 7, 33));
    assert(postfit_value_cache_get(&cache, 7)); /* Oversize misses leave old data. */
    fail_at = calls + 1;
    assert(!postfit_value_cache_put(&cache, 2, 16));
    assert(!postfit_value_cache_get(&cache, 7));
    check_accounting(&cache);
    fail_at = 0;
    assert(postfit_value_cache_put(&cache, 3, 16));
    assert(!postfit_value_cache_put(&cache, 8, 16));
    assert(!postfit_value_cache_put(&cache, 3, 0));
    postfit_value_cache_destroy(&cache);
    postfit_value_cache_destroy(&cache);
    assert(n_live == 0);
    for (size_t failure = 0; failure < 30; ++failure) {
        calls = 0;
        fail_at = failure;
        assert(postfit_value_cache_init(&cache, 8, budget, 32, test_allocate, test_release));
        for (size_t iteration = 0; iteration < 30; ++iteration) {
            (void)postfit_value_cache_put(&cache, iteration % 8, 1 + iteration % 32);
            check_accounting(&cache);
        }
        postfit_value_cache_destroy(&cache);
        assert(n_live == 0);
    }
    assert(postfit_value_cache_init(&cache, SIZE_MAX, 64, 0, malloc, free));
    assert(!postfit_value_cache_put(&cache, 0, 1));
    postfit_value_cache_destroy(&cache);
    assert(postfit_value_cache_init(&cache, 8, 0, 0, malloc, free));
    assert(!postfit_value_cache_put(&cache, 0, 1));
    postfit_value_cache_destroy(&cache);
    assert(!postfit_value_cache_init(&cache, 8, 1, 2, malloc, free));
    postfit_value_cache_destroy(&cache);
    puts("FIFO eviction, accounting, allocation-failure recovery and overflow guards pass.");
    return 0;
}
