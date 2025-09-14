#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <assert.h>
#include <stdbool.h>
#include <intrin.h>

#define ERR_RET_VAL 1

void check( bool success, const char *format, ... ) {
    if ( success ) return;
    va_list args;
    va_start( args, format );
    vfprintf( stderr, format, args );
    va_end( args ); 
    fprintf( stderr, "\n");
    __debugbreak();
    exit( ERR_RET_VAL );
}