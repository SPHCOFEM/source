/*!
  @par Revision history:
  - 18.01.2001, c
*/
#include <stdarg.h>
#include <string.h>
#include <errno.h>

#include "allocator.h"
#include "errutils.h"
#include "ioutils.h"

IOConf g_ioc;

#undef __FUNC__
#define __FUNC__ "io_init"
/*!
  @par Revision history:
  - 18.01.2001, c
*/
int io_init( IOConf *obj )
{
  memset( obj, 0, sizeof( IOConf ) );

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "io_alloc"
/*!
  @par Revision history:
  - 18.01.2001, c
  - 08.11.2001
*/
int io_alloc( IOConf *obj )
{
  int i, j;

  io_init( obj );

  for (i = 0; i < IO_LevelMax; i++) {
    obj->tabbing[i] = (char *) malloc( (2 * i + 1) * sizeof( char ) );
    for (j = 0; j < i; j++) {
      obj->tabbing[i][2*j] = ' ';
      obj->tabbing[i][2*j+1] = '.';
/*        if (j % 2) obj->tabbing[i][j] = '.'; */
/*        else obj->tabbing[i][j] = ' '; */
    }
    obj->tabbing[i][2*i] = '\0';
  }
  obj->level = 0;

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "io_free"
/*!
  @par Revision history:
  - 18.01.2001, c
*/
int io_free( IOConf *obj )
{
  int i;

  if (obj == 0) return( RET_OK );

  for (i = 0; i < IO_LevelMax; i++) {
    free( obj->tabbing[i] );
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "io_setProgName"
/*!
  @par Revision history:
  - 18.01.2001, c
*/
int io_setProgName( IOConf *obj, char *path )
{
  int i, l;

  if (obj == 0) {
    ERR_SetAndDie( ERR_NotAllocated );
  }

  l = strlen( path );
  i = l; 
  while (i >= 0) {
    if (path[i] == '/') break;
    i--;
  }
  strncpy( obj->programName, &path[i+1], l - i ); // There should be null.

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "errput"
/*!
  @par Revision history:
  - 18.01.2001, c
  - 25.01.2001
*/
void errput( const char *what, ... )
{
  va_list ap;
  char buf[256]; /* !!! */
  FILE *file;

  if ((!g_ioc.stderrTermOutput) && (!g_ioc.stderrLog)) {
    return;
  }
  va_start( ap, what );
  sprintf( buf, "%s: **ERROR** %s -> %s", g_ioc.programName,
	   g_ioc.tabbing[g_ioc.level], what );
  if (g_ioc.stderrTermOutput) {
    vfprintf( stdout, buf, ap );
  }
  if (g_ioc.stderrLog) {
    file = fopen( g_ioc.stderrLogFileName, "a" );
    vfprintf( file, buf, ap );
    fclose( file );
  }
  va_end( ap );
}

#undef __FUNC__
#define __FUNC__ "output"
/*!
  @par Revision history:
  - 18.01.2001, c
  - 15.11.2001
*/
void output( const char *what, ... )
{
  va_list ap;
  char buf[256]; /* !!! */
  FILE *file;

  if (g_ioc.supress) return;
  if ((!g_ioc.stdoutTermOutput) && (!g_ioc.stdoutLog)) {
    return;
  }
  va_start( ap, what );
  sprintf( buf, "%s:%s %s", g_ioc.programName,
	   g_ioc.tabbing[g_ioc.level], what );
  if (g_ioc.stdoutTermOutput) {
    vfprintf( stdout, buf, ap );
  }
  if (g_ioc.stdoutLog) {
    file = fopen( g_ioc.stdoutLogFileName, "a" );
    vfprintf( file, buf, ap );
    fclose( file );
  }
  va_end( ap );
}

#undef __FUNC__
#define __FUNC__ "io_up"
/*!
  @par Revision history:
  - 18.01.2001, c
  - 08.02.2001
*/
void io_up( const char *fnname )
{
  if (g_ioc.level >= 1) g_ioc.level--;
  else {
    errput( "too low io level!\n" );
  }

  if (fnname == 0) return;
  output( "Leaving %s().\n", fnname );
}

#undef __FUNC__
#define __FUNC__ "io_down"
/*!
  @par Revision history:
  - 18.01.2001, c
  - 08.02.2001
*/
void io_down( const char *fnname )
{
  if (fnname == 0) return;
  output( "Entering %s()...\n", fnname );

  if (g_ioc.level < (IO_LevelMax - 1)) g_ioc.level++;
  else {
    errput( "too much io levels!\n" );
  }
}

#undef __FUNC__
#define __FUNC__ "io_suppress"
/*!
  @par Revision history:
  - 15.11.2001, c
  - 16.11.2001
*/
void io_suppress()
{
  g_ioc.supress++;
}

#undef __FUNC__
#define __FUNC__ "io_allow"
/*!
  @par Revision history:
  - 15.11.2001, c
  - 16.11.2001
*/
void io_allow()
{
  g_ioc.supress--;
  if (g_ioc.supress < 0) {
    errput( "too low io level!\n" );
    g_ioc.supress = 0;
  }
}
