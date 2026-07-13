#pragma once
#include <string>
#include <format>

using namespace std;

// Wrap a compound expression (one containing top-level + or -) in parentheses.
// Products and single terms never contain spaces; sums and differences always do.
string paren( const string &t ) {
    if ( t.find( ' ' ) != string::npos ) return format( "( {} )", t );
    return t;
}

string addTerms( const string &t1, const string &t2 ) {
    if ( t1 == "0" ) return t2;
    if ( t2 == "0" ) return t1;
    return format( "{} + {}", t1, t2 );
}

// Negates a term. If it's already a bare negation (e.g. "-1" or "-sin_q0",
// as opposed to a parenthesized sum like "( a + b )"), strips the leading '-'
// instead of stacking another one -- "--1" is not valid C (it parses as
// pre-decrement of a non-lvalue), so double negatives must cancel here.
string negate( const string &t ) {
    if ( t == "0" ) return "0";
    string p = paren( t );
    if ( p[ 0 ] == '-' ) return p.substr( 1 );
    return "-" + p;
}

string subTerms( const string &t1, const string &t2 ) {
    if ( t2 == "0" ) return t1;
    if ( t1 == "0" ) return negate( t2 );
    return format( "{} - {}", t1, paren( t2 ) );
}

string multFactors( const string &f1, const string &f2 ) {
    if ( f1 == "0" ) return "0";
    if ( f2 == "0" ) return "0";
    if ( f1 == "1" ) return f2;
    if ( f2 == "1" ) return f1;
    return format( "{}*{}", paren( f1 ), paren( f2 ) );
}
