#pragma once

#include <string>
#include "symmath.h"
#include "prln.h"

using namespace std;

struct Svec {
    string baseName;
    string v[ 6 ];

    Svec() : Svec( "" ) {}
    Svec( string &&name ) : Svec( name ) {}
    Svec( string &name )
        :baseName( name )
    {
        for ( int i = 0 ; i < 6 ; i++ ) {
            v[ i ] = "0";
        }
    }

    string& operator[]( int i ) {
        return v[ i ];
    }

    Svec operator*( const string& s ) {
        Svec p;
        for ( int i = 0 ; i < 6 ; i++ ) {
            p.v[ i ] = multFactors( v[ i ], s );
        }
        return p;
    }
    
    void emit() {
        for ( int i = 0 ; i < 6 ; i++ ) {
            string outStr = format( "double {}{} = {};", baseName, i, v[ i ] );
            prln( "%s", outStr.c_str() );
        }
    }
};