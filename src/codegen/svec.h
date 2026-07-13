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
    
    // Emits each element as a named temp, then replaces the element with the
    // temp's name (literals 0 and 1 stay as-is to preserve sparsity).
    void emit() {
        for ( int i = 0 ; i < 6 ; i++ ) {
            string outStr = format( "double {}{} = {};", baseName, i, v[ i ] );
            prln( "%s", outStr.c_str() );
            if ( v[ i ] != "0" && v[ i ] != "1" ) v[ i ] = format( "{}{}", baseName, i );
        }
    }
};

string svecDot( Svec &a, Svec &b ) {
    string sum = "0";
    for ( int i = 0 ; i < 6 ; i++ ) sum = addTerms( sum, multFactors( a.v[i], b.v[i] ) );
    return sum;
}

// out = a x b, for two 3-vectors of expression strings.
void cross3( const string a[3], const string b[3], string out[3] ) {
    out[0] = subTerms( multFactors( a[1], b[2] ), multFactors( a[2], b[1] ) );
    out[1] = subTerms( multFactors( a[2], b[0] ), multFactors( a[0], b[2] ) );
    out[2] = subTerms( multFactors( a[0], b[1] ), multFactors( a[1], b[0] ) );
}

// Motion cross product, v x m (RBDA eq 2.33): both operands are motion Svecs.
// Used for the v_i x vJ_i term in the acceleration recursion.
Svec svCrossMotion( string name, Svec &v, Svec &m ) {
    Svec s( name );
    const string *w   = &v.v[0];
    const string *rho = &v.v[3];
    const string *w2  = &m.v[0];
    const string *rho2 = &m.v[3];

    string wxw2[3], wxrho2[3], rhoxw2[3];
    cross3( w, w2, wxw2 );
    cross3( w, rho2, wxrho2 );
    cross3( rho, w2, rhoxw2 );

    for ( int i = 0 ; i < 3 ; i++ ) {
        s.v[i]   = wxw2[i];
        s.v[i+3] = addTerms( wxrho2[i], rhoxw2[i] );
    }
    return s;
}

// Force cross product, v x* f (RBDA eq 2.34): v is a motion Svec, f a force
// Svec. Used for the v_i x* (I_i * v_i) bias term in the force recursion.
Svec svCrossForce( string name, Svec &v, Svec &f ) {
    Svec s( name );
    const string *w   = &v.v[0];
    const string *rho = &v.v[3];
    const string *n   = &f.v[0];
    const string *f2  = &f.v[3];

    string wxn[3], rhoxf2[3], wxf2[3];
    cross3( w, n, wxn );
    cross3( rho, f2, rhoxf2 );
    cross3( w, f2, wxf2 );

    for ( int i = 0 ; i < 3 ; i++ ) {
        s.v[i]   = addTerms( wxn[i], rhoxf2[i] );
        s.v[i+3] = wxf2[i];
    }
    return s;
}

Svec svecAdd( string name, Svec &a, Svec &b ) {
    Svec s( name );
    for ( int i = 0 ; i < 6 ; i++ ) s.v[i] = addTerms( a.v[i], b.v[i] );
    return s;
}

Svec svecSub( string name, Svec &a, Svec &b ) {
    Svec s( name );
    for ( int i = 0 ; i < 6 ; i++ ) s.v[i] = subTerms( a.v[i], b.v[i] );
    return s;
}