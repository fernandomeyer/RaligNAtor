#ifndef GT_INTERFACE_H
#define GT_INTERFACE_H
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#define GT_DBL_MAX_ABS_ERROR 1.0E-100
#define GT_DBL_MAX_REL_ERROR 1.0E-8

#ifndef MAX
#define MAX(X,Y) (((X) > (Y)) ? (X) : (Y))
#endif

#define GT_EXIT_PROGRAMMING_ERROR  2

#define gt_assert(expression) assert(expression)
#define gt_free(ptr) free(ptr)
#define gt_malloc(size) malloc(size)
#define gt_realloc(ptr,size) realloc(ptr,size)
#define gt_realloc_mem(ptr,size,file,line) realloc(ptr,size)
#define gt_free_func free

typedef uint8_t       GtUchar;
typedef unsigned long GtUlong;
typedef FILE GtError;
typedef FILE GtLogger;

/* Unused function arguments should be annotated with this macro to get rid of
   compiler warnings. */
#define GT_UNUSED \
        __attribute__ ((unused)) /*@unused@*/

static inline void gt_error_set(GtError *err, const char *format, ...)
{
  va_list ap;
  if (!err) return;
  va_start(ap, format);
  (void) vfprintf(err, format, ap);
  va_end(ap);
}

static inline void gt_logger_log(GtLogger *err, const char *format, ...)
{
  va_list ap;
  if (!err) return;
  va_start(ap, format);
  (void) vfprintf(err, format, ap);
  va_end(ap);
}

static inline bool gt_double_relative_equal(double d1, double d2)
{
  double relerr;
  if (fabs(d1 - d2) < GT_DBL_MAX_ABS_ERROR)
    return true;
  if (fabs(d2) > fabs(d1))
    relerr = fabs((d1 - d2) / d2);
  else
    relerr = fabs((d1 - d2) / d1);
  if (relerr <= GT_DBL_MAX_REL_ERROR)
    return true;
  return false;
}

static inline bool gt_double_equals_double(double d1, double d2)
{
  return gt_double_relative_equal(d1, d2);
}

#endif
