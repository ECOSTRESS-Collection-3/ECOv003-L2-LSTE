// error.c
// See error.h for credits and documentation
#include "bool.h"
#include "error.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <linux/limits.h>

static char flogn[PATH_MAX];
FILE *flog = NULL;
static bool echo_stdout = false;

void set_log_filename(const char *logname)
{
    strncpy(flogn, logname, PATH_MAX);
    flogn[PATH_MAX-1] = '\0';
    if (flog) fclose(flog);
    flog = fopen(flogn, "w");
    if (!flog)
    {
        fprintf(stderr, "ERROR: cannot open output log file %s\n",
                 logname);
        exit(-1);
    }
}

void echo_log_to_stdout()
{
    echo_stdout = true;
}

void logmsg(Severity severe, const char *srcfile, int linenum, 
            int ernum, const char *message, ...)
{
    va_list valist;
    const char *severe_str;
    char msgout[1024];
    char msgfmt[1024];
    char tmpBff[PATH_MAX];
    size_t check_4_truncated = -1;
    
    // prepare va_list struct
    va_start(valist, message);

    // Set severity string
    if (severe == Error)
        severe_str = "ERROR";
    else if (severe == Warning)
        severe_str = "WARNING";
    else
        severe_str = "INFO";

    // prepare msgfmt
    check_4_truncated = snprintf(msgfmt, sizeof(msgfmt), "%s at %s:%d -- %s\n",
            severe_str, srcfile, linenum, message);
    if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(msgfmt))
    {
        ERREXIT(113, "Buffer is too small.", NONE);
    }
    // Apply var args to output string.
    check_4_truncated = vsnprintf(msgout, sizeof(msgout), msgfmt, valist);
    if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(msgout))
    {
        ERREXIT(113, "Buffer is too small.", NONE);
    }
    // Cleanup va_list memory
    va_end(valist);
    
    // Output message to log file.
    if (flog)
    {        
        memset(tmpBff, 0, sizeof(tmpBff));
        check_4_truncated = snprintf(tmpBff, sizeof(tmpBff), "%s", msgout);
        if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(tmpBff))
        {
            ERREXIT(113, "Buffer is too small.", NONE);
        }
        fprintf(flog, tmpBff);
    }
    if (echo_stdout)
    {        
        memset(tmpBff, 0, sizeof(tmpBff));
        check_4_truncated = snprintf(tmpBff, sizeof(tmpBff), "%s", msgout);
        if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(tmpBff))
        {
            ERREXIT(113, "Buffer is too small.", NONE);
        }
        fprintf(stdout, tmpBff);
    }

    // Exit if error
    if (severe == Error)
    {
        fflush(stdout);
        fflush(stderr);
        fflush(flog);
        if (ernum == 0)
            ernum = -1; // Zero is a "success" exit code.
        exit(ernum); 
    }
}
