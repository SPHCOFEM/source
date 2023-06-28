/*!
  @par Revision history:
  - 18.01.2001, c
*/
#ifndef _IOUTILS_H_
#define _IOUTILS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

//#include "types.h"

#define IO_LevelMax 10

typedef struct
{
  int stdoutLog;
  int stderrLog;
  char stdoutLogFileName[256]; /* !!! */
  char stderrLogFileName[256]; /* !!! */

  int stdoutTermOutput;
  int stderrTermOutput;
  int supress;

  int level;
  char *tabbing[IO_LevelMax];

  char programName[64]; /* !!! */
} IOConf;

extern IOConf g_ioc;

int io_init( IOConf *obj );
int io_alloc( IOConf *obj );
int io_free( IOConf *obj );
int io_setProgName( IOConf *obj, char *path );
void errput( const char *what, ... );
void output( const char *what, ... );

//  void local_debput( char *fname, int line, char *buffer );
//  #ifdef DEBUG
//  #define debput( buffer ) (void) local_debput( __FILE__, __LINE__, buffer );
//  #elde
//  #define debput( buffer ) ((void) 0)
//  #endif

void io_up( const char *fnname );
void io_down( const char *fnname );
void io_suppress();
void io_allow();

#ifdef __cplusplus
}
#endif

#endif /* Header */
