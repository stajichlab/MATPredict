/* glib-lite implementation: see glib.h. Semantics follow GLib 2.x for the functions
 * exonerate 2.4.0 calls. */
#include "glib.h"
#include <ctype.h>
#include <stdarg.h>
#include <strings.h>

/* ---- allocation ---- */

void *g_malloc(gsize n) {
    void *p = n ? malloc(n) : NULL;
    if (n && !p) { fprintf(stderr, "glib-lite: out of memory (%zu bytes)\n", n); abort(); }
    return p;
}

void *g_malloc0(gsize n) {
    void *p = n ? calloc(1, n) : NULL;
    if (n && !p) { fprintf(stderr, "glib-lite: out of memory (%zu bytes)\n", n); abort(); }
    return p;
}

void *g_realloc(void *p, gsize n) {
    if (!n) { free(p); return NULL; }
    void *q = realloc(p, n);
    if (!q) { fprintf(stderr, "glib-lite: out of memory (%zu bytes)\n", n); abort(); }
    return q;
}

void g_free(void *p) { free(p); }

/* ---- diagnostics ---- */

void g_print(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vprintf(fmt, ap);
    va_end(ap);
}

static void vlog(const char *prefix, const char *fmt, va_list ap) {
    fprintf(stderr, "%s", prefix);
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
}

void g_message(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vlog("** Message: ", fmt, ap);
    va_end(ap);
}

void g_warning(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vlog("** WARNING **: ", fmt, ap);
    va_end(ap);
}

void g_critical(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vlog("** CRITICAL **: ", fmt, ap);
    va_end(ap);
}

void g_error(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vlog("** ERROR **: ", fmt, ap);
    va_end(ap);
    abort(); /* GLib aborts on g_error */
}

void g_assertion_fail(const char *expr, const char *file, int line) {
    fprintf(stderr, "%s:%d: assertion failed: (%s)\n", file, line, expr);
    abort();
}

unsigned int g_log_set_handler(const char *domain, int levels, void *handler, void *data) {
    (void)domain; (void)levels; (void)handler; (void)data;
    return 0;
}

void g_on_error_stack_trace(const char *prg_name) { (void)prg_name; }

/* exonerate is single-threaded here; GLib's thread setup is a no-op. */
int g_thread_supported(void) { return 1; }
void g_thread_init(void *vtable) { (void)vtable; }

/* ---- strings ---- */

gchar *g_strdup(const gchar *s) {
    if (!s) return NULL;
    gsize n = strlen(s) + 1;
    gchar *d = g_malloc(n);
    memcpy(d, s, n);
    return d;
}

gchar *g_strndup(const gchar *s, gsize n) {
    if (!s) return NULL;
    gsize len = 0;
    while (len < n && s[len]) ++len;
    gchar *d = g_malloc(len + 1);
    memcpy(d, s, len);
    d[len] = 0;
    return d;
}

gchar *g_strdup_vprintf(const gchar *fmt, va_list ap) {
    va_list cp;
    va_copy(cp, ap);
    int n = vsnprintf(NULL, 0, fmt, cp);
    va_end(cp);
    gchar *d = g_malloc((gsize)n + 1);
    vsnprintf(d, (gsize)n + 1, fmt, ap);
    return d;
}

gchar *g_strdup_printf(const gchar *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    gchar *d = g_strdup_vprintf(fmt, ap);
    va_end(ap);
    return d;
}

gchar *g_strconcat(const gchar *first, ...) {
    gsize total = 0;
    va_list ap;
    va_start(ap, first);
    for (const gchar *s = first; s; s = va_arg(ap, const gchar *)) total += strlen(s);
    va_end(ap);
    gchar *d = g_malloc(total + 1);
    gchar *p = d;
    va_start(ap, first);
    for (const gchar *s = first; s; s = va_arg(ap, const gchar *)) {
        size_t n = strlen(s);
        memcpy(p, s, n);
        p += n;
    }
    va_end(ap);
    *p = 0;
    return d;
}

gchar *g_strnfill(gsize len, gchar fill) {
    gchar *d = g_malloc(len + 1);
    memset(d, fill, len);
    d[len] = 0;
    return d;
}

gchar *g_strup(gchar *s) {
    for (gchar *p = s; *p; ++p) *p = (gchar)toupper((unsigned char)*p);
    return s;
}

gchar *g_strdown(gchar *s) {
    for (gchar *p = s; *p; ++p) *p = (gchar)tolower((unsigned char)*p);
    return s;
}

gchar *g_strstrip(gchar *s) {
    gchar *start = s;
    while (*start && isspace((unsigned char)*start)) ++start;
    gsize n = strlen(start);
    while (n > 0 && isspace((unsigned char)start[n - 1])) --n;
    memmove(s, start, n);
    s[n] = 0;
    return s;
}

gint g_strcasecmp(const gchar *a, const gchar *b) { return strcasecmp(a, b); }

gchar *g_strjoinv(const gchar *sep, gchar **strv) {
    if (!sep) sep = "";
    gsize total = 0, seplen = strlen(sep), count = 0;
    for (gchar **v = strv; v && *v; ++v, ++count) total += strlen(*v);
    if (count > 1) total += seplen * (count - 1);
    gchar *d = g_malloc(total + 1);
    gchar *p = d;
    for (gchar **v = strv; v && *v; ++v) {
        if (v != strv) { memcpy(p, sep, seplen); p += seplen; }
        size_t n = strlen(*v);
        memcpy(p, *v, n);
        p += n;
    }
    *p = 0;
    return d;
}

/* g_strsplit: delim is a whole string, not a set of characters. With max_tokens > 0,
 * at most max_tokens pieces are returned and the last one holds the rest. Empty
 * pieces are kept. The result is NULL-terminated. */
gchar **g_strsplit(const gchar *s, const gchar *delim, gint max_tokens) {
    gsize dlen = strlen(delim);
    gint cap = 4, n = 0;
    gchar **out = g_new(gchar *, cap + 1);
    const gchar *p = s;
    for (;;) {
        const gchar *hit = NULL;
        if (dlen && !(max_tokens > 0 && n == max_tokens - 1)) hit = strstr(p, delim);
        gsize len = hit ? (gsize)(hit - p) : strlen(p);
        if (n == cap) { cap *= 2; out = g_renew(gchar *, out, cap + 1); }
        out[n++] = g_strndup(p, len);
        if (!hit) break;
        p = hit + dlen;
    }
    out[n] = NULL;
    return out;
}

void g_strfreev(gchar **strv) {
    if (!strv) return;
    for (gchar **v = strv; *v; ++v) g_free(*v);
    g_free(strv);
}

const gchar *g_getenv(const gchar *name) { return getenv(name); }

/* ---- GString ---- */

GString *g_string_sized_new(gsize dfl_size) {
    GString *s = g_new(GString, 1);
    s->allocated_len = dfl_size < 2 ? 2 : dfl_size;
    s->str = g_malloc(s->allocated_len);
    s->len = 0;
    s->str[0] = 0;
    return s;
}

GString *g_string_new(const gchar *init) {
    GString *s = g_string_sized_new(init ? strlen(init) + 1 : 2);
    if (init) g_string_append(s, init);
    return s;
}

static void gs_reserve(GString *s, gsize extra) {
    gsize need = s->len + extra + 1;
    if (need > s->allocated_len) {
        gsize cap = s->allocated_len;
        while (cap < need) cap *= 2;
        s->str = g_realloc(s->str, cap);
        s->allocated_len = cap;
    }
}

GString *g_string_append(GString *s, const gchar *val) {
    gsize n = strlen(val);
    gs_reserve(s, n);
    memcpy(s->str + s->len, val, n + 1);
    s->len += n;
    return s;
}

GString *g_string_append_c(GString *s, gchar c) {
    gs_reserve(s, 1);
    s->str[s->len++] = c;
    s->str[s->len] = 0;
    return s;
}

void g_string_sprintfa(GString *s, const gchar *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    gchar *d = g_strdup_vprintf(fmt, ap);
    va_end(ap);
    g_string_append(s, d);
    g_free(d);
}

GString *g_string_truncate(GString *s, gsize len) {
    if (len < s->len) {
        s->len = len;
        s->str[len] = 0;
    }
    return s;
}

gchar *g_string_free(GString *s, gboolean free_segment) {
    gchar *str = s->str;
    g_free(s);
    if (free_segment) { g_free(str); return NULL; }
    return str;
}

/* ---- GString chunks ---- */

struct _GStringChunk {
    gchar **strings;
    gsize n, cap;
};

GStringChunk *g_string_chunk_new(gsize size) {
    (void)size;
    GStringChunk *c = g_new0(GStringChunk, 1);
    return c;
}

gchar *g_string_chunk_insert(GStringChunk *chunk, const gchar *s) {
    if (chunk->n == chunk->cap) {
        chunk->cap = chunk->cap ? chunk->cap * 2 : 64;
        chunk->strings = g_renew(gchar *, chunk->strings, chunk->cap);
    }
    chunk->strings[chunk->n] = g_strdup(s);
    return chunk->strings[chunk->n++];
}

void g_string_chunk_free(GStringChunk *chunk) {
    if (!chunk) return;
    for (gsize i = 0; i < chunk->n; ++i) g_free(chunk->strings[i]);
    g_free(chunk->strings);
    g_free(chunk);
}

/* ---- GArray ---- */

GArray *g_array_new(gboolean zero_terminated, gboolean clear, guint element_size) {
    GArray *a = g_new0(GArray, 1);
    a->zero_terminated = zero_terminated;
    a->clear = clear;
    a->elt_size = element_size;
    a->alloc = 0;
    a->data = NULL;
    return a;
}

static void ga_reserve(GArray *a, guint len) {
    if (len + 1 > a->alloc) {
        guint cap = a->alloc ? a->alloc : 4;
        while (cap < len + 1) cap *= 2;
        a->data = g_realloc(a->data, (gsize)cap * a->elt_size);
        memset(a->data + (gsize)a->alloc * a->elt_size, 0, (gsize)(cap - a->alloc) * a->elt_size);
        a->alloc = cap;
    }
}

GArray *g_array_append_vals(GArray *a, gconstpointer data, guint len) {
    ga_reserve(a, a->len + len);
    memcpy(a->data + (gsize)a->len * a->elt_size, data, (gsize)len * a->elt_size);
    a->len += len;
    return a;
}

GArray *g_array_set_size(GArray *a, guint len) {
    ga_reserve(a, len);
    if (len > a->len)
        memset(a->data + (gsize)a->len * a->elt_size, 0, (gsize)(len - a->len) * a->elt_size);
    a->len = len;
    return a;
}

gchar *g_array_free(GArray *a, gboolean free_segment) {
    gchar *data = a->data;
    g_free(a);
    if (free_segment) { g_free(data); return NULL; }
    return data;
}

/* ---- GPtrArray ---- */

GPtrArray *g_ptr_array_new(void) {
    GPtrArray *a = g_new0(GPtrArray, 1);
    return a;
}

static void pa_reserve(GPtrArray *a, guint len) {
    if (len > a->alloc) {
        guint cap = a->alloc ? a->alloc : 8;
        while (cap < len) cap *= 2;
        a->pdata = g_renew(gpointer, a->pdata, cap);
        for (guint i = a->alloc; i < cap; ++i) a->pdata[i] = NULL;
        a->alloc = cap;
    }
}

void g_ptr_array_add(GPtrArray *a, gpointer p) {
    pa_reserve(a, a->len + 1);
    a->pdata[a->len++] = p;
}

void g_ptr_array_set_size(GPtrArray *a, gint len) {
    pa_reserve(a, (guint)len);
    for (guint i = a->len; i < (guint)len; ++i) a->pdata[i] = NULL;
    a->len = (guint)len;
}

/* Removes the first match by moving the last element into its place, as GLib does. */
gboolean g_ptr_array_remove_fast(GPtrArray *a, gpointer p) {
    for (guint i = 0; i < a->len; ++i) {
        if (a->pdata[i] == p) {
            a->pdata[i] = a->pdata[a->len - 1];
            a->len--;
            return TRUE;
        }
    }
    return FALSE;
}

gpointer *g_ptr_array_free(GPtrArray *a, gboolean free_segment) {
    gpointer *data = a->pdata;
    g_free(a);
    if (free_segment) { g_free(data); return NULL; }
    return data;
}

/* ---- GMemChunk ---- */

struct _GMemChunk {
    gsize atom_size;
    gint type;
    gpointer *atoms;
    gsize n, cap;
};

GMemChunk *g_mem_chunk_new(const gchar *name, gint atom_size, gsize area_size, gint type) {
    (void)name; (void)area_size;
    GMemChunk *m = g_new0(GMemChunk, 1);
    m->atom_size = (gsize)atom_size;
    m->type = type;
    return m;
}

gpointer g_mem_chunk_alloc(GMemChunk *m) {
    if (m->n == m->cap) {
        m->cap = m->cap ? m->cap * 2 : 256;
        m->atoms = g_renew(gpointer, m->atoms, m->cap);
    }
    gpointer p = g_malloc0(m->atom_size);
    m->atoms[m->n++] = p;
    return p;
}

/* G_ALLOC_ONLY atoms stay valid until the chunk is destroyed; freeing them is a no-op,
 * as in GLib. */
void g_mem_chunk_free(GMemChunk *m, gpointer mem) {
    if (m->type == G_ALLOC_AND_FREE) {
        for (gsize i = 0; i < m->n; ++i)
            if (m->atoms[i] == mem) {
                g_free(mem);
                m->atoms[i] = m->atoms[--m->n];
                return;
            }
    }
}

void g_mem_chunk_destroy(GMemChunk *m) {
    if (!m) return;
    for (gsize i = 0; i < m->n; ++i) g_free(m->atoms[i]);
    g_free(m->atoms);
    g_free(m);
}

/* ---- paths and environment ---- */

const gchar *g_basename(const gchar *path) {
    const gchar *slash = strrchr(path, '/');
    return slash ? slash + 1 : path;
}

const gchar *g_get_tmp_dir(void) {
    const gchar *dir = getenv("TMPDIR");
    return (dir && *dir) ? dir : "/tmp";
}

/* ---- GTree: a sorted array of key/value pairs, binary search ---- */

typedef struct {
    gpointer key;
    gpointer value;
} gt_pair;

struct _GTree {
    GCompareFunc compare;
    gt_pair *pairs;
    gsize n, cap;
};

GTree *g_tree_new(GCompareFunc key_compare_func) {
    GTree *t = g_new0(GTree, 1);
    t->compare = key_compare_func;
    return t;
}

/* Index of the first pair whose key is >= key. */
static gsize gt_lower(GTree *t, gconstpointer key, gboolean *found) {
    gsize lo = 0, hi = t->n;
    while (lo < hi) {
        gsize mid = lo + (hi - lo) / 2;
        int c = t->compare(t->pairs[mid].key, key);
        if (c < 0) lo = mid + 1;
        else hi = mid;
    }
    *found = lo < t->n && t->compare(t->pairs[lo].key, key) == 0;
    return lo;
}

void g_tree_insert(GTree *t, gpointer key, gpointer value) {
    gboolean found;
    gsize at = gt_lower(t, key, &found);
    if (found) {
        t->pairs[at].value = value; /* GLib keeps one entry per key and replaces the value */
        return;
    }
    if (t->n == t->cap) {
        t->cap = t->cap ? t->cap * 2 : 32;
        t->pairs = g_renew(gt_pair, t->pairs, t->cap);
    }
    memmove(&t->pairs[at + 1], &t->pairs[at], (t->n - at) * sizeof(gt_pair));
    t->pairs[at].key = key;
    t->pairs[at].value = value;
    t->n++;
}

gpointer g_tree_lookup(GTree *t, gconstpointer key) {
    gboolean found;
    gsize at = gt_lower(t, key, &found);
    return found ? t->pairs[at].value : NULL;
}

gboolean g_tree_remove(GTree *t, gconstpointer key) {
    gboolean found;
    gsize at = gt_lower(t, key, &found);
    if (!found) return FALSE;
    memmove(&t->pairs[at], &t->pairs[at + 1], (t->n - at - 1) * sizeof(gt_pair));
    t->n--;
    return TRUE;
}

gint g_tree_nnodes(GTree *t) { return (gint)t->n; }

/* Only in-order traversal is implemented: it is the only order exonerate uses. A callback
 * that returns TRUE stops the walk, as in GLib. */
void g_tree_traverse(GTree *t, GTraverseFunc func, GTraverseType order, gpointer data) {
    if (order != G_IN_ORDER) { fprintf(stderr, "glib-lite: only G_IN_ORDER traversal\n"); abort(); }
    for (gsize i = 0; i < t->n; ++i)
        if (func(t->pairs[i].key, t->pairs[i].value, data)) break;
}

void g_tree_destroy(GTree *t) {
    if (!t) return;
    g_free(t->pairs);
    g_free(t);
}
