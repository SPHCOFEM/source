#ifndef _MEMUTILS_H_
#define _MEMUTILS_H_

/*!
  @file
  @short General usability things dealing with the memory - header file.
  Inspired by PETSC memory routines.
*/

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>

#include "types.h"

#define AL_CookieValue   0xf0e0d0c9
#define AL_AlreadyFreed  0x0f0e0d9c

/*!
  @par Revision history:
  - 25.05.2001, c
*/
typedef struct _AllocSpace {
    unsigned long   size;
    int             id;
    int             lineNo;
    char            *fileName;
    char            *funName;
    char            *dirName;
    unsigned long   cookie;
    struct _AllocSpace *next,*prev;
} AllocSpace;

#define AL_HeaderDoubles      sizeof(AllocSpace)/sizeof(double)+1

/*!
  This union is used to insure that the block passed to the user is
  aligned on a double boundary

  @par Revision history:
  - 25.05.2001, c
*/
typedef union {
    AllocSpace sp;
    double  v[AL_HeaderDoubles];
} AllocSpaceAlign;

/*! @enum AllocMode
  @par Revision history:
  - 15.04.2001, c
*/
typedef enum {
  AL_Alloc, AL_Free, AL_Realloc
} AllocMode;

void *mem_allocMem( size_t size, int lineNo, char *funName,
		    char *fileName, char *dirName );
void *mem_reallocMem( void *pp, size_t size, int lineNo, char *funName,
		      char *fileName, char *dirName );
void mem_freeMem( void *pp, int lineNo, char *funName,
		  char *fileName, char *dirName );
void mem_checkIntegrity();
void mem_statistics( int lineNo, char *funName,
		     char *fileName, char *dirName );
int32 mem_print( FILE *file, int32 mode );
int32 mem_printSome( FILE *file, int32 mode, int32 num );
int32 mem_freeGarbage();

void swapPointers( void **p1, void **p2 );
void dummy_mem_stats();

/*!
  @name Allocator wrapping macros
  @par Revision history:
  - 25.09.2000, 0.26.0, c
  - 11.10.2000, 0.27.0
  - 05.01.2001
  - 22.02.2001
  - 25.05.2001
  - 18.09.2001
*/
/*@{*/
#define freeMem( p ) do {\
  mem_freeMem( p, __LINE__, __FUNC__, __FILE__, __SDIR__ ); p = 0; } while (0)
#define allocMem( size )\
  (mem_allocMem( size, __LINE__, __FUNC__, __FILE__, __SDIR__ ))
#define reallocMem( p, size ) (\
  mem_reallocMem( p, size, __LINE__, __FUNC__, __FILE__, __SDIR__ ))
#define printMemStats()\
  (mem_statistics( __LINE__, __FUNC__, __FILE__, __SDIR__ ))
#define checkMemoryIntegrity()\
  (mem_checkIntegrity( __LINE__, __FUNC__, __FILE__, __SDIR__ ))
#define createMem( Type ) (Type *) allocMem( sizeof( Type ) )
#define createMemMore( Type, num ) (Type *) allocMem( (num) * sizeof( Type ) )
#define destroyMem( p ) freeMem( p )
/*@}*/

#ifdef __cplusplus
}
#endif

/*!
  @name Delete wrapping macro
  @par Revision history:
  - 08.09.2000, 0.25.1, c
  - 12.09.2000, 0.25.1, c
*/
/*@{*/
#define freez(x) do { if (x) { delete (x); x = 0; } } while (0)
/*@}*/

#endif /* Header */
