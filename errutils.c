#include "ioutils.h"
#include "errutils.h"

/*!
  @par Revision history:
  - 03.01.2001, c
  - 05.01.2001
  - 06.02.2001
  - 07.02.2001
  - 15.05.2001
  - 25.05.2001
  - 05.06.2001
  - 18.12.2001
  - 13.01.2002
*/
int g_error;

int g_nerr = 23;

char *g_errList[] = {
  "No error :)",
  "Corrupted file!",
  "End of file!",
  "Unknown keyword!",
  "Memory allocation failed!",
  "Structure not allocated!",
  "Array sizes do no match properly!",
  "Bad array limit!",
  "Cannot open file!",
  "Data read failed (data corrupted)!",
  "Wrong object (argument) type!",
  "Null pointer !",
  "Verification failed!",
  "Comparison error!",
  "Switch failed!",
  "Item not found!",
  "Numerical error! (e.g. division by zero)",
  "Python interpret error!",
  "Expat error!",
  "External library error!",
  "Corrupted memory!",
  "Parse error!",
  "Operation failed!"
};

/*!
  @short Print error message determined by 'g_error'.

  @par Revision history:
  - 03.01.2001, c
  - 04.01.2001
  - 18.01.2001
  - 25.01.2001
*/
void printError( const char *s, const char *sourcedir, const char *fname,
		 const char *funname, int line )
{
  if ((g_error >= 0) && (g_error < g_nerr)) { 
    if (s == 0) {
      errput( "%s\n", g_errList[g_error] );
      output( " - error in %s%s, %s(), line %d.\n",
	      sourcedir, fname, funname, line );
    } else {
      errput( "%s\n", g_errList[g_error] );
      output( " - error in %s%s, %s(), line %d:\n%s.\n",
	      sourcedir, fname, funname, line, s );
    }
  }
}
