#pragma once
#include <string>
#include <format>

using namespace std;

string addTerms( const string &t1, const string &t2 ) {
    if ( t1 == "0" ) return t2;
    if ( t2 == "0" ) return t1;
    return format( "{} + {}", t1, t2 );
}

string multFactors( const string &f1, const string &f2 ) {
    if ( f1 == "0" ) return "0";
    if ( f2 == "0" ) return "0";
    if ( f1 == "1" ) return f2;
    if ( f2 == "1" ) return f1;
    return format( "{}*{}", f1, f2 );
}
