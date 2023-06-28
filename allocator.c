/*!
  @file
  @short General usability things dealing with the memory.
*/

#include <stdio.h>
#include <stdlib.h>

#include "allocator.h"
#include "errutils.h"
#include "ioutils.h"

/*!
  @par Revision history:
  - 25.05.2001, c
*/
//static int32 al_expUsage;
static int32 al_curUsage;
static int32 al_maxUsage;
static int32 al_frags;
static AllocSpace *al_head = 0;


#undef __FUNC__
#define __FUNC__ "mem_allocMem"
/*!
  Allocates the memory and updates memory usage statistics.

  @par Revision history:
  - 05.01.2001, c
  - 25.05.2001
*/
void *mem_allocMem( size_t size, int lineNo, char *funName,
		    char *fileName, char *dirName )
{
  char *p;
  int32 hsize = sizeof( AllocSpaceAlign );
  int32 tsize;
  float64 *endptr;
  AllocSpace *head;

  if (size == 0) {
    errput( "%s, %s, %s, %d: zero allocation!\n",
	    dirName, fileName, funName, lineNo );
    ERR_SetAndDie( ERR_Alloc );
  }

  size = (size / 8) * 8 + 8;
  
  tsize = size + hsize + sizeof( float64 );
  if ((p = (char *) malloc( tsize )) == 0) {
    errput( "%s, %s, %s, %d: error allocating %d bytes (current: %d).\n",
	    dirName, fileName, funName, lineNo, size, al_curUsage );
    ERR_SetAndDie( ERR_Alloc );
  }
  head = (AllocSpace *) p;
  p += hsize;

  if (al_head) al_head->prev = head;
  head->next     = al_head;
  al_head         = head;
  head->prev     = 0;
  head->size     = size;
  head->id       = 1234567;
  head->lineNo   = lineNo;
  head->fileName = fileName;
  head->funName  = funName;
  head->dirName  = dirName;
  head->cookie   = AL_CookieValue;
  endptr         = (float64 *) ((((char *) p) + size));
  endptr[0]      = (float64) AL_CookieValue;

/*    output( "%d _> ptr: %p, head: %p, end %p\n", hsize, p, head, endptr ); */

  al_curUsage += size;
  if (al_curUsage > al_maxUsage) {
    al_maxUsage = al_curUsage;
  }
  al_frags++;

  return( (void *) p );
}

#undef __FUNC__
#define __FUNC__ "mem_reallocMem"
/*!
  Reallocates the memory and updates memory usage statistics.

  @par Revision history:
  - 05.01.2001, c
  - 25.05.2001
  - 07.06.2001
  - 05.11.2001
*/
void *mem_reallocMem( void *pp, size_t size, int lineNo, char *funName,
		      char *fileName, char *dirName )
{
  char *p = (char *) pp;
  char *p1;
  int32 hsize = sizeof( AllocSpaceAlign );
  int32 tsize;
  float64 *endptr;
  AllocSpace *head;
  char *phead;

  if (pp) {
    phead = p - hsize;
    head = (AllocSpace *) phead;
    if (head->cookie != AL_CookieValue) {
      errput( "%s, %s, %s, %d: ptr: %p, cookie: %d\n",
	      dirName, fileName, funName, lineNo,
	      p, head->cookie );
      if (head->cookie == AL_AlreadyFreed) {
	errput( "memory was already freed!\n" );
      }
      ERR_SetAndDie( ERR_CorruptedMemory );
    }

    endptr = (float64 *) (p + head->size);
/*      output( "%d _> ptr: %p, head: %p, end %p\n", hsize, p, head, endptr ); */
    if (endptr[0] != AL_CookieValue) {
      errput( "%s, %s, %s, %d:\n",
	      dirName, fileName, funName, lineNo );
      if (endptr[0] == AL_AlreadyFreed) {
	errput( "already freed!\n" );
      } else {
	errput( "damaged tail!\n" );
      }
      ERR_SetAndDie( ERR_CorruptedMemory );
    }
    endptr[0] = AL_AlreadyFreed;
    al_curUsage -= head->size;

    if (head->prev) head->prev->next = head->next;
    else al_head = head->next;

    if (head->next) head->next->prev = head->prev;

    if (size == 0) {
      al_frags--;
      
      free( phead );
      p1 = 0;
    } else {
      tsize = size + hsize + sizeof( float64 );
      if ((p1 = (char *) realloc( phead, tsize )) == 0) {
	errput( "%s, %s, %s, %d: error reallocating %d bytes (current: %d).\n",
		dirName, fileName, funName, lineNo, size, al_curUsage );
	ERR_SetAndDie( ERR_Alloc );
      }
      head = (AllocSpace *) p1;
      p1 += hsize;

      if (al_head) al_head->prev = head;
      head->next     = al_head;
      al_head         = head;
      head->prev     = 0;
      head->size     = size;
      head->id       = 1234567;
      head->lineNo   = lineNo;
      head->fileName = fileName;
      head->funName  = funName;
      head->dirName  = dirName;
      head->cookie   = AL_CookieValue;
      endptr         = (float64 *) (p1 + size);
      endptr[0]      = (float64) AL_CookieValue;

/*        output( "%d _> ptr: %p, head: %p, end %p\n", hsize, p1, head, endptr ); */

      al_curUsage += size;
      if (al_curUsage > al_maxUsage) {
	al_maxUsage = al_curUsage;
      }
    }
  } else {
    p1 = (char *) mem_allocMem( size, lineNo, funName,
				fileName, dirName );
  }

  return( (void *) p1 );
}

#undef __FUNC__
#define __FUNC__ "mem_freeMem"
/*!
  Frees the memory and updates memory usage statistics.

  @par Revision history:
  - 05.01.2001, c
  - 25.05.2001
  - 05.11.2001
*/
void mem_freeMem( void *pp, int lineNo, char *funName,
		  char *fileName, char *dirName )
{
  char *p = (char *) pp;
  int32 hsize = sizeof( AllocSpaceAlign );
  float64 *endptr;
  AllocSpace *head;
  char *phead;

  if (p == 0) return;

  phead = p - hsize;
  head = (AllocSpace *) phead;
  if (head->cookie != AL_CookieValue) {
    errput( "%s, %s, %s, %d: ptr: %p, cookie: %d\n",
	    dirName, fileName, funName, lineNo,
	    p, head->cookie );
    if (head->cookie == AL_AlreadyFreed) {
      errput( "memory was already freed!\n" );
    }
    ERR_SetAndDie( ERR_CorruptedMemory );
  }
  head->cookie = AL_AlreadyFreed;

  endptr = (float64 *) (p + head->size);
  if (endptr[0] != AL_CookieValue) {
    errput( "%s %s %s %d:\n",
	    dirName, fileName, funName, lineNo );
    if (endptr[0] == AL_AlreadyFreed) {
      errput( "already freed!\n" );
    } else {
      errput( "damaged tail!\n" );
    }
    ERR_SetAndDie( ERR_CorruptedMemory );
  }

  endptr[0] = (float64) AL_AlreadyFreed;

  al_curUsage -= head->size;
  al_frags--;

  if (head->prev) head->prev->next = head->next;
  else al_head = head->next;

  if (head->next) head->next->prev = head->prev;

  free( phead );
}

#undef __FUNC__
#define __FUNC__ "mem_checkIntegrity"
/*!
  Frees the memory and updates memory usage statistics.

  @par Revision history:
  - 07.06.2001, c
  - 18.09.2001
  - 05.11.2001
*/
void mem_checkIntegrity( int lineNo, char *funName,
			 char *fileName, char *dirName )
{
  char *p, *pp;
  int32 cnt, allocated;
  int32 hsize = sizeof( AllocSpaceAlign );
  float64 *endptr;
  AllocSpace *head = al_head;

  output( "checking memory integrity in\n" );
  output( "%s, %s, %s(), %d:\n",
	  dirName, fileName, funName, lineNo, al_maxUsage, al_curUsage );
  output( "allocated memory: %d records, usage: %d, max: %d\n",
	  al_frags, al_curUsage, al_maxUsage );
  if (head == 0) {
    goto end_label;
  }

  cnt = 0;
  allocated = 0;
  while (head) {
    p = (char *) head;

    pp = p + hsize;
    if (head->cookie != AL_CookieValue) {
      errput( "ptr: %p, ptrhead: %p, cookie: %d\n",
	      pp, p , head->cookie );
      if (head->cookie == AL_AlreadyFreed) {
	errput( "memory was already freed!\n" );
      }
      ERR_SetAndDie( ERR_CorruptedMemory );
    }

    endptr = (float64 *) (pp + head->size);
    if (endptr[0] != AL_CookieValue) {
      output( "  %s, %s, %s, %d: size: %d, ptr: %p\n",
	      head->dirName, head->fileName, head->funName, head->lineNo,
	      head->size, pp );
      if (endptr[0] == AL_AlreadyFreed) {
	errput( "already freed!\n" );
      } else {
	errput( "damaged tail!\n" );
      }
      ERR_SetAndDie( ERR_CorruptedMemory );
    }
    cnt++;
    allocated += head->size;
    if (cnt > al_frags) {
      errput( "damaged allocation record (overrun)!\n" );
      ERR_SetAndDie( ERR_CorruptedMemory );
    }
    head = head->next;
  }
  if (cnt < al_frags) {
    errput( "damaged allocation record (underrun)!\n" );
    ERR_SetAndDie( ERR_CorruptedMemory );
  }
  if (allocated != al_curUsage) {
    errput( "memory leak!? (%d == %d)\n", allocated, al_curUsage );
    ERR_SetAndDie( ERR_CorruptedMemory );
  }

 end_label:
  output( "memory OK.\n" );
}

#undef __FUNC__
#define __FUNC__ "mem_statistics"
/*!
  Prints memory usage statistics.
  @par Revision history:
  - 25.05.2001
  - 18.09.2001
*/
void mem_statistics( int lineNo, char *funName,
		     char *fileName, char *dirName )
{
  output( "%s, %s, %s(), %d: memory max: %d, current: %d\n",
	  dirName, fileName, funName, lineNo, al_maxUsage, al_curUsage );
}

#undef __FUNC__
#define __FUNC__ "mem_print"
/*!
  Prints memory usage statistics.
  @par Revision history:
  - 25.05.2001, c
*/
int32 mem_print( FILE *file, int32 mode )
{
  int32 cnt = 0;
  int32 hsize = sizeof( AllocSpaceAlign );
  AllocSpace *head = al_head;
  char *p;

  mode = 0;
  output( "allocated memory: %d records, usage: %d, max: %d\n",
	  al_frags, al_curUsage, al_maxUsage );
  if (head == 0) {
    goto end_label;
  }

  while (head) {
    p = (char *) head;
/*      output( "%d _> head: %p, end %p\n", hsize, p, */
/*  	    p + hsize + head->size ); */
    output( "  %s, %s, %s, %d: size: %d, ptr: %p\n",
	    head->dirName, head->fileName, head->funName, head->lineNo,
	    head->size, p + hsize );
    cnt++;
    if (cnt > al_frags) {
      errput( "damaged allocation record (overrun)!\n" );
      ERR_SetAndDie( ERR_CorruptedMemory );
    }
    head = head->next;
  }
  if (cnt < al_frags) {
    errput( "damaged allocation record (underrun)!\n" );
    ERR_SetAndDie( ERR_CorruptedMemory );
  }

 end_label:
  output( "done.\n" );

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "mem_printSome"
/*!
  Prints memory usage statistics.
  @par Revision history:
  - 25.05.2001, c
*/
int32 mem_printSome( FILE *file, int32 mode, int32 num )
{
  int32 cnt = 0;
  int32 hsize = sizeof( AllocSpaceAlign );
  AllocSpace *head = al_head;
  char *p;

  mode = 0;
  output( "allocated memory: %d records, usage: %d, max: %d\n",
	  al_frags, al_curUsage, al_maxUsage );
  output( "printing max: %d\n", num );
  if (head == 0) {
    goto end_label;
  }

  while (head) {
    p = (char *) head;
/*      output( "%d _> head: %p, end %p\n", hsize, p, */
/*  	    p + hsize + head->size ); */
    output( "  %s, %s, %s, %d: size: %d, ptr: %p\n",
	    head->dirName, head->fileName, head->funName, head->lineNo,
	    head->size, p + hsize );
    cnt++;
    if (cnt > al_frags) {
      errput( "damaged allocation record (overrun)!\n" );
      ERR_SetAndDie( ERR_CorruptedMemory );
    }
    if (cnt == num) break;
    head = head->next;
  }

 end_label:
  output( "done.\n" );

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "mem_freeGarbage"
/*!
  Free all memory records.
  @par Revision history:
  - 27.05.2001, c
*/
int32 mem_freeGarbage()
{
  int32 cnt = 0, frags = al_frags;
  int32 hsize = sizeof( AllocSpaceAlign );
  char *p;

  output( "freeing garbage.\n" );
  while (al_head) {
    p = (char *) al_head + hsize;
    freeMem( p );
/*      output( "  %s, %s, %s, %d: size: %d, ptr: %p\n", */
/*  	    head->dirName, head->fileName, head->funName, head->lineNo, */
/*  	    head->size, p + hsize ); */
    cnt++;
/*      printf( "%d %d\n", cnt, frags ); */
    if (cnt > frags) {
      errput( "damaged allocation record (overrun)!\n" );
      ERR_SetAndDie( ERR_CorruptedMemory );
    }
  }
  if (cnt < frags) {
    errput( "damaged allocation record (underrun)!\n" );
    ERR_SetAndDie( ERR_CorruptedMemory );
  }

  return( RET_OK );
}

void swapPointers( void **p1, void **p2 )
{
  void *p;

  p = *p1;
  *p1 = *p2;
  *p2 = p;
}

#ifdef WITH_PYTHON
#include "pym.h"
#endif

#undef __FUNC__
#define __FUNC__ "exitfun"
/*!
  @par Revision history:
  - 10.01.2001, c
  - 18.01.2001
  - 15.05.2001
  - 27.05.2001
  - 07.06.2001
*/
void dummy_mem_stats()
{
  if (g_error == ERR_CorruptedMemory)
    output( "memory will NOT be printed and freed!\n"  );
  else {
    mem_printSome( stdout, 0, 10 );
    mem_freeGarbage();
  }

  printMemStats();
  io_free( &g_ioc );
#ifdef WITH_PYTHON
  Py_Finalize();
#endif
}
