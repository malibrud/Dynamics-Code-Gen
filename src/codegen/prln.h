#pragma once
#include <stdio.h>
#include <stdarg.h>

static FILE *prln_file;
static int prln_indent_level = 0;
const int INDENT = 4;

void prln( const char *format, ... ) {
    for ( int i = 0 ; i < prln_indent_level * INDENT ; i++ ) putc( ' ', prln_file );
    va_list args;
    va_start( args, format );

    vfprintf( prln_file, format, args );
    putc( '\n', prln_file );

    va_end( args );
}
