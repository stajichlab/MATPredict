/* glib-lite: the subset of GLib that exonerate 2.4.0 uses, for builds where the real
 * GLib cannot be cross-compiled (it always builds GIO and GObject, which need libffi).
 * Behaviour follows GLib 2.x for every function listed in glib-lite.c. Not a general
 * replacement for GLib. Checked by building exonerate against both libraries and
 * comparing the output. */
#ifndef GLIB_LITE_H
#define GLIB_LITE_H

#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

/* ---- basic types ---- */
typedef int gboolean;
typedef char gchar;
typedef unsigned char guchar;
typedef short gshort;
typedef unsigned short gushort;
typedef int gint;
typedef unsigned int guint;
typedef long glong;
typedef unsigned long gulong;
typedef float gfloat;
typedef double gdouble;
typedef int8_t gint8;
typedef uint8_t guint8;
typedef int16_t gint16;
typedef uint16_t guint16;
typedef int32_t gint32;
typedef uint32_t guint32;
typedef int64_t gint64;
typedef uint64_t guint64;
typedef size_t gsize;
typedef ssize_t gssize;
typedef void *gpointer;
typedef const void *gconstpointer;
typedef int (*GCompareFunc)(gconstpointer a, gconstpointer b);
typedef void (*GDestroyNotify)(gpointer data);

#define TRUE 1
#define FALSE 0
#define G_GNUC_EXTENSION
#define G_GINT64_CONSTANT(val) ((gint64)(val##LL))
#define G_HAVE_GINT64 1
#define G_DIR_SEPARATOR '/'
#define G_DIR_SEPARATOR_S "/"

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define ABS(a) (((a) < 0) ? -(a) : (a))
/* Byte-order conversion. Every target here is little-endian (x86-64, arm64, wasm32), so
 * the big-endian forms are byte swaps. */
#define GUINT16_TO_BE(x) ((guint16)__builtin_bswap16((guint16)(x)))
#define GUINT32_TO_BE(x) ((guint32)__builtin_bswap32((guint32)(x)))
#define GUINT64_TO_BE(x) ((guint64)__builtin_bswap64((guint64)(x)))
#define GUINT16_FROM_BE(x) GUINT16_TO_BE(x)
#define GUINT32_FROM_BE(x) GUINT32_TO_BE(x)
#define GUINT64_FROM_BE(x) GUINT64_TO_BE(x)
#define GINT_TO_POINTER(i) ((gpointer)(glong)(i))
#define GPOINTER_TO_INT(p) ((gint)(glong)(p))
#define GUINT_TO_POINTER(u) ((gpointer)(gulong)(u))
#define GPOINTER_TO_UINT(p) ((guint)(gulong)(p))

#define GLIB_MAJOR_VERSION 2
#define GLIB_MINOR_VERSION 82
#define GLIB_MICRO_VERSION 0
typedef int GLogLevelFlags;
#define G_LOG_LEVEL_ERROR (1 << 2)
#define G_LOG_FLAG_FATAL (1 << 4)

/* ---- allocation ---- */
void *g_malloc(gsize n);
void *g_malloc0(gsize n);
void *g_realloc(void *p, gsize n);
void g_free(void *p);

#define g_new(type, n) ((type *)g_malloc(sizeof(type) * (gsize)(n)))
#define g_new0(type, n) ((type *)g_malloc0(sizeof(type) * (gsize)(n)))
#define g_renew(type, p, n) ((type *)g_realloc((p), sizeof(type) * (gsize)(n)))

/* ---- diagnostics ---- */
void g_print(const char *fmt, ...);
void g_message(const char *fmt, ...);
void g_warning(const char *fmt, ...);
void g_critical(const char *fmt, ...);
void g_error(const char *fmt, ...);
void g_assertion_fail(const char *expr, const char *file, int line);
unsigned int g_log_set_handler(const char *domain, int levels, void *handler, void *data);
void g_on_error_stack_trace(const char *prg_name);
int g_thread_supported(void);
void g_thread_init(void *vtable);

#ifndef G_DISABLE_ASSERT
#define g_assert(expr) do { if (!(expr)) g_assertion_fail(#expr, __FILE__, __LINE__); } while (0)
#define g_assert_not_reached() g_assertion_fail("code should not be reached", __FILE__, __LINE__)
#else
#define g_assert(expr) ((void)0)
#define g_assert_not_reached() ((void)0)
#endif

/* ---- strings ---- */
gchar *g_strdup(const gchar *s);
gchar *g_strndup(const gchar *s, gsize n);
gchar *g_strdup_printf(const gchar *fmt, ...);
gchar *g_strdup_vprintf(const gchar *fmt, va_list ap);
gchar *g_strconcat(const gchar *first, ...);
gchar *g_strnfill(gsize len, gchar fill);
gchar *g_strup(gchar *s);
gchar *g_strdown(gchar *s);
gchar *g_strstrip(gchar *s);
gint g_strcasecmp(const gchar *a, const gchar *b);
gchar *g_strjoinv(const gchar *sep, gchar **strv);
gchar **g_strsplit(const gchar *s, const gchar *delim, gint max_tokens);
void g_strfreev(gchar **strv);
const gchar *g_getenv(const gchar *name);

/* ---- GString ---- */
typedef struct {
    gchar *str;
    gsize len;
    gsize allocated_len;
} GString;

GString *g_string_new(const gchar *init);
GString *g_string_sized_new(gsize dfl_size);
GString *g_string_append(GString *s, const gchar *val);
GString *g_string_append_c(GString *s, gchar c);
void g_string_sprintfa(GString *s, const gchar *fmt, ...);
GString *g_string_truncate(GString *s, gsize len);
gchar *g_string_free(GString *s, gboolean free_segment);

/* ---- GString chunks: strings owned by a chunk, all freed together ---- */
typedef struct _GStringChunk GStringChunk;
GStringChunk *g_string_chunk_new(gsize size);
gchar *g_string_chunk_insert(GStringChunk *chunk, const gchar *s);
void g_string_chunk_free(GStringChunk *chunk);

/* ---- GArray: element_size bytes per element, data in ->data ---- */
typedef struct {
    gchar *data;
    guint len;
    gboolean zero_terminated;
    gboolean clear;
    guint elt_size;
    guint alloc;
} GArray;

GArray *g_array_new(gboolean zero_terminated, gboolean clear, guint element_size);
GArray *g_array_append_vals(GArray *a, gconstpointer data, guint len);
GArray *g_array_set_size(GArray *a, guint len);
gchar *g_array_free(GArray *a, gboolean free_segment);

#define g_array_append_val(a, v) g_array_append_vals((a), &(v), 1)
#define g_array_index(a, t, i) (((t *)(void *)(a)->data)[(i)])

/* ---- GPtrArray: pointers in ->pdata ---- */
typedef struct {
    gpointer *pdata;
    guint len;
    guint alloc;
} GPtrArray;

GPtrArray *g_ptr_array_new(void);
void g_ptr_array_add(GPtrArray *a, gpointer p);
void g_ptr_array_set_size(GPtrArray *a, gint len);
gboolean g_ptr_array_remove_fast(GPtrArray *a, gpointer p);
gpointer *g_ptr_array_free(GPtrArray *a, gboolean free_segment);

/* ---- GMemChunk: fixed-size atom allocation. Atoms are allocated with malloc and tracked,
 * so g_mem_chunk_destroy frees them all. ---- */
typedef enum {
    G_ALLOC_ONLY,
    G_ALLOC_AND_FREE
} GAllocatorType;
typedef struct _GMemChunk GMemChunk;
GMemChunk *g_mem_chunk_new(const gchar *name, gint atom_size, gsize area_size, gint type);
gpointer g_mem_chunk_alloc(GMemChunk *mem_chunk);
void g_mem_chunk_free(GMemChunk *mem_chunk, gpointer mem);
void g_mem_chunk_destroy(GMemChunk *mem_chunk);
#define g_chunk_new(type, chunk) ((type *)g_mem_chunk_alloc(chunk))

/* ---- paths and environment ---- */
const gchar *g_basename(const gchar *path);
const gchar *g_get_tmp_dir(void);

/* ---- GTree: sorted by a compare function; kept as a sorted array ---- */
typedef struct _GTree GTree;
GTree *g_tree_new(GCompareFunc key_compare_func);
void g_tree_insert(GTree *tree, gpointer key, gpointer value);
gpointer g_tree_lookup(GTree *tree, gconstpointer key);
void g_tree_destroy(GTree *tree);
gboolean g_tree_remove(GTree *tree, gconstpointer key);
gint g_tree_nnodes(GTree *tree);

typedef enum {
    G_IN_ORDER,
    G_PRE_ORDER,
    G_POST_ORDER,
    G_LEVEL_ORDER
} GTraverseType;
typedef gboolean (*GTraverseFunc)(gpointer key, gpointer value, gpointer data);
void g_tree_traverse(GTree *tree, GTraverseFunc func, GTraverseType order, gpointer data);

#endif
