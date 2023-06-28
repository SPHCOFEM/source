#ifndef _ERRUTILS_H_
#define _ERRUTILS_H_

#ifdef __cplusplus
extern "C" {
#endif

//#include "types.h"

#include <stdio.h>

#ifndef __SDIR__
#define __SDIR__ "unknowndirectory/"
#endif

#ifndef __FUNC__
#define __FUNC__ "unknownfunction"
#endif

/*!

  @par Revision history:
  - 03.01.2001, c
  - 04.01.2001
  - 05.01.2001
  - 16.01.2001
  - 05.02.2001
  - 06.02.2001
  - 07.02.2001
  - 09.02.2001
  - 15.05.2001
  - 25.05.2001
  - 05.06.2001
  - 18.12.2001
  - 13.01.2002
*/
typedef enum {
  ERR_None,
  ERR_CorruptedFile,
  ERR_EndOfFile,
  ERR_UnknownKeyword,
  ERR_Alloc,
  ERR_NotAllocated,
  ERR_BadMatch,
  ERR_BadLimit,
  ERR_FileOpen,
  ERR_DataRead,
  ERR_BadObj,
  ERR_NullPointer,
  ERR_VerificationFail,
  ERR_Compare,
  ERR_Switch,
  ERR_Search,
  ERR_Numeric,
  ERR_Python,
  ERR_Expat,
  ERR_External,
  ERR_CorruptedMemory,
  ERR_Parse,
  ERR_OperationFail,
} ErrorName;

typedef enum {
  RET_OK,
  RET_Fail
} ReturnStatus;

extern int g_error;
extern int g_nerr;
extern char *g_errList[];

void printError( const char *s, const char *sourcedir, const char *fname,
		 const char *funname, int line );
#define PrintError( s ) (\
  printError( (s), __SDIR__, __FILE__, __FUNC__, __LINE__ ))
#define PrintErrorAndDie( s ) do {\
  printError( (s), __SDIR__, __FILE__, __FUNC__, __LINE__ );\
  printf( "Exiting...\n" );\
  exit( g_error );\
} while (0);

#define ERR_GotoEnd( i ) do { g_error = (i); goto end_label; } while (0)
#define ERR_Check( ret ) do {\
  (ret) = (g_error == ERR_None) ? RET_OK : RET_Fail;\
} while (0)
#define ERR_SetAndDie( i ) do { g_error = (i); PrintErrorAndDie( 0 );\
} while (0)
#define ERR_CheckDie do { if (g_error != ERR_None) PrintErrorAndDie( 0 );\
} while (0)
#define ERR_Set( i ) (g_error = (i))
#define ERR_SetPrintAndReturn( i ) do { g_error = (i); PrintError( 0 );\
  return( RET_Fail );\
} while (0)


#ifdef __cplusplus
}
#endif

#endif /* Header */
