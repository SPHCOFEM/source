#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <stdlib.h>
#include <errno.h>

#include "allocator.h"
#include "errutils.h"
#include "ioutils.h"

#define PI 3.141592653589793
#define EPS 1e-9
#define MAX 1e9
#define TMP "sphcofem.tmp"
#define OUTFILE "sphcofem.out"
#define TXTFILE "sphcofem.txt"
#define MSGFILE "sphcofem.msg"
#define SIGFILE "signal"
#define LOGNAME "sphcofem.log"

#define sqr(x) ((x)*(x))
#define norm(x,y,z) (sqrt(sqr(x)+sqr(y)+sqr(z)))
