#pragma once

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>
#include "system.h"
#include "prln.h"

#define VAR_SIZE 64
typedef struct {
    char s[ VAR_SIZE ];
} Var;

Var _Var( const char *f, ... ) {
    Var v;
    va_list args;
    va_start( args, f );
    vsnprintf( v.s, VAR_SIZE, f, args );
    va_end( args );
    return v;
}

typedef struct {
    Var v[ 6 ];
} Svec_;

typedef struct {
    Var E[ 3 ][ 3 ];
    Var r[ 3 ];
} xForm;

Svec_ genS( int iBody, int jType ) {
    Svec_ S;
    switch( jType ) {
        case J_REVOLUTE:
            S.v[ 0 ] = _Var( "0" );
            S.v[ 1 ] = _Var( "0" );
            S.v[ 2 ] = _Var( "1" );
            S.v[ 3 ] = _Var( "0" );
            S.v[ 4 ] = _Var( "0" );
            S.v[ 5 ] = _Var( "0" );
            break;
        default:
            assert( false );
    }
    return S;
}

Svec_ genVJ( int iBody, Svec_ *S, Var *qd ) {
    Svec_ vJ;
    for ( int i = 0 ; i < 6 ; i++ ) {
        bool isElement0 = !strcmp( S->v[ i ].s, "0" );
        vJ.v[ i ] = isElement0 ? _Var( "0") : *qd;
    }
    return vJ;
}

void gen_fwd_dyn_crba( FILE *f, System *sys ) {
    prln_file = f;

    int NB = sys->nBod;
    prln( "void fwd_dyn_crba( double *q, double *qd, double *tau, double *out_qdd ) {" );
    prln_indent_level++;

    prln( "// The outward pass of the computation of C.  This is the first for loop of Table 5.1 in RBDA. ");
    prln( "" );
    for ( int i = 0 ; i < NB ; i++ ) {
        prln( "// Body %d", i );
        int jType = sys->jtype[ i ];
        Var q = _Var( "q%d", i );
        Var qd = _Var( "qd%d", i );
        Svec_ S = genS( i, jType );
        Svec_ vJ = genVJ( i, &S, &qd );

        prln( "double %s = q[ %d ];", q.s, i );
        prln( "double %s = qd[ %d ];", qd.s, i );
        prln( "double S%d[ 6 ] = { %s, %s, %s, %s, %s, %s };", 
            i,
            S.v[0], S.v[1], S.v[2], S.v[3], S.v[4], S.v[5]
        );
        prln( "double sin_%s = sin( %s );", q.s, q.s );
        prln( "double cos_%s = cos( %s );", q.s, q.s );

        sys->xFromParent[ i ].emit();
        sys->xFromBase[ i ].emit();

        prln( "" );
    }

    prln_indent_level--;
    prln( "}" );
}