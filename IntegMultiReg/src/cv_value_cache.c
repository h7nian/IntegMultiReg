#include "cv_value_cache.h"
#include <limits.h>
#include <stdint.h>
#include <string.h>

static void postfit_value_cache_remove(postfit_value_cache *cache, int key)
{
    postfit_value_cache_entry *entry = &cache->entries[key];
    if (entry->previous >= 0) cache->entries[entry->previous].next = entry->next;
    else cache->first = entry->next;
    if (entry->next >= 0) cache->entries[entry->next].previous = entry->previous;
    else cache->last = entry->previous;
    cache->used_bytes -= entry->bytes;
    cache->release(entry->data);
    entry->data = NULL;
    entry->bytes = 0;
    entry->previous = entry->next = -1;
}

int postfit_value_cache_init(postfit_value_cache *cache, size_t n_entries,
                            size_t budget_bytes, size_t external_bytes,
                            void *(*allocate)(size_t), void (*release)(void *))
{
    memset(cache, 0, sizeof(*cache));
    cache->first = cache->last = -1;
    if (!allocate || !release || external_bytes > budget_bytes) return 0;
    cache->n_entries = n_entries;
    cache->budget_bytes = budget_bytes;
    cache->external_bytes = cache->used_bytes = external_bytes;
    cache->allocate = allocate;
    cache->release = release;
    if (!n_entries || n_entries > INT_MAX ||
        n_entries > SIZE_MAX / sizeof(postfit_value_cache_entry)) return 1;
    size_t entry_bytes = n_entries * sizeof(postfit_value_cache_entry);
    if (entry_bytes > budget_bytes - external_bytes) return 1;
    cache->entries = allocate(entry_bytes);
    if (!cache->entries) return 1;
    cache->entry_bytes = entry_bytes;
    cache->used_bytes += entry_bytes;
    for (size_t key = 0; key < n_entries; ++key) {
        cache->entries[key].data = NULL;
        cache->entries[key].bytes = 0;
        cache->entries[key].previous = cache->entries[key].next = -1;
    }
    return 1;
}

void *postfit_value_cache_get(const postfit_value_cache *cache, size_t key)
{
    if (!cache->entries || key >= cache->n_entries) return NULL;
    return cache->entries[key].data;
}

void *postfit_value_cache_put(postfit_value_cache *cache, size_t key, size_t bytes)
{
    if (!cache->entries || key >= cache->n_entries || !bytes ||
        bytes > cache->budget_bytes - cache->external_bytes - cache->entry_bytes)
        return NULL;
    if (cache->entries[key].data) postfit_value_cache_remove(cache, (int)key);
    /* Evict before allocating so live requested bytes never exceed the budget.
     * A hit does not change insertion order: eviction is deterministic FIFO. */
    while (bytes > cache->budget_bytes - cache->used_bytes)
        postfit_value_cache_remove(cache, cache->first);
    void *data = cache->allocate(bytes);
    if (!data) return NULL;
    postfit_value_cache_entry *entry = &cache->entries[key];
    entry->data = data;
    entry->bytes = bytes;
    entry->previous = cache->last;
    entry->next = -1;
    if (cache->last >= 0) cache->entries[cache->last].next = (int)key;
    else cache->first = (int)key;
    cache->last = (int)key;
    cache->used_bytes += bytes;
    return data;
}

void postfit_value_cache_destroy(postfit_value_cache *cache)
{
    if (!cache->entries) return;
    while (cache->first >= 0) postfit_value_cache_remove(cache, cache->first);
    cache->release(cache->entries);
    cache->entries = NULL;
    cache->entry_bytes = 0;
    cache->used_bytes = cache->external_bytes;
}
