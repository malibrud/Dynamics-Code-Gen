#pragma once
#include <inttypes.h>
#include <string>
#include <format>
#include "svec.h"
#include "prln.h"
#include "symmath.h"

using namespace std;
using namespace std::string_literals;

typedef int16_t Xtype;

#define XFORM_AXIS_R ( (Xtype)0x20 )
#define XFORM_AXIS_T ( (Xtype)0x40 )
#define XFORM_AXIS_P ( (Xtype)0x80 )
#define XFORM_AXIS_TR ( XFORM_AXIS_T | XFORM_AXIS_R )

#define XFORM_IDENT      ( (Xtype)0 )
#define XFORM_6DOF       ( (Xtype)1 )
#define XFORM_3DOF_ROT   ( (Xtype)2 )
#define XFORM_3DOF_TRANS ( (Xtype)3 )

// Transformations around the X, Y, and Z axes.
#define XFORM_RX ( XFORM_AXIS_R | (Xtype)4 )
#define XFORM_RY ( XFORM_AXIS_R | (Xtype)5 )
#define XFORM_RZ ( XFORM_AXIS_R | (Xtype)6 )

// Transformations along the X, Y, and Z axes.
#define XFORM_TX ( XFORM_AXIS_T | (Xtype)8 )
#define XFORM_TY ( XFORM_AXIS_T | (Xtype)9 )
#define XFORM_TZ ( XFORM_AXIS_T | (Xtype)10 )

// Transformations with translation along an axis then rotation around an axis
#define XFORM_TX_RX ( XFORM_AXIS_TR | (Xtype)12 )
#define XFORM_TX_RY ( XFORM_AXIS_TR | (Xtype)13 )
#define XFORM_TX_RZ ( XFORM_AXIS_TR | (Xtype)14 )

#define XFORM_TY_RX ( XFORM_AXIS_TR | (Xtype)16 )
#define XFORM_TY_RY ( XFORM_AXIS_TR | (Xtype)17 )
#define XFORM_TY_RZ ( XFORM_AXIS_TR | (Xtype)18 )

#define XFORM_TZ_RX ( XFORM_AXIS_TR | (Xtype)21 )
#define XFORM_TZ_RY ( XFORM_AXIS_TR | (Xtype)22 )
#define XFORM_TZ_RZ ( XFORM_AXIS_TR | (Xtype)23 )

// Transformations in the planes XY, YZ, ZX
#define XFORM_PXY ( XFORM_AXIS_P | (Xtype)25 )
#define XFORM_PYZ ( XFORM_AXIS_P | (Xtype)26 )
#define XFORM_PZX ( XFORM_AXIS_P | (Xtype)27 )

// Transformations defined by D-H parameters
#define XFORM_DH_LINK  (Xtype)29
#define XFORM_DH_JOINT (Xtype)30

static std::string blank = ""s;
struct Xform {
    Xtype type;
    std::string baseName;
    std::string E[ 3 ][ 3 ];
    std::string r[ 3 ];

    Xform(): Xform( XFORM_IDENT, "", blank ) {}
    Xform( Xtype type, std::string name, std::string q )
        :baseName( name ) 
    {
        using namespace std::literals::string_literals;
        using namespace std;

        string cq = "cos_"s + q;
        string sq = "sin_"s + q;
        switch( type ) {
        case XFORM_IDENT:
            E[0][0] = "1"; E[0][1] = "0"; E[0][2] = "0";
            E[1][0] = "0"; E[1][1] = "1"; E[1][2] = "0";
            E[2][0] = "0"; E[2][1] = "0"; E[2][2] = "1";

            r[0] = "0";  r[1] = "0";  r[2] = "0";
            break;
        
        case XFORM_RX:
            E[0][0] = "1"; E[0][1] = "0"; E[0][2] = "0";
            E[1][0] = "0"; E[1][1] = cq;  E[1][2] = "-"s + sq;
            E[2][0] = "0"; E[2][1] = sq;  E[2][2] = cq;

            r[0] = "0";  r[1] = "0";  r[2] = "0";
            break;

        case XFORM_RY:
            E[0][0] = cq;        E[0][1] = "0"; E[0][2] = sq;
            E[1][0] = "0";       E[1][1] = "1"; E[1][2] = "0";
            E[2][0] = "-"s + sq; E[2][1] = "0"; E[2][2] = cq;

            r[0] = "0";  r[1] = "0";  r[2] = "0";
            break;

        case XFORM_RZ:
            E[0][0] = cq;  E[0][1] = "-"s + sq; E[0][2] = "0";
            E[1][0] = sq;  E[1][1] = cq;        E[1][2] = "0";
            E[2][0] = "0"; E[2][1] = "0";       E[2][2] = "1";

            r[0] = "0";  r[1] = "0";  r[2] = "0";
            break;

        // Translations along one axis by a symbolic amount q -- these are
        // the prismatic joint transforms.
        case XFORM_TX:
            E[0][0] = "1"; E[0][1] = "0"; E[0][2] = "0";
            E[1][0] = "0"; E[1][1] = "1"; E[1][2] = "0";
            E[2][0] = "0"; E[2][1] = "0"; E[2][2] = "1";

            r[0] = q;  r[1] = "0";  r[2] = "0";
            break;

        case XFORM_TY:
            E[0][0] = "1"; E[0][1] = "0"; E[0][2] = "0";
            E[1][0] = "0"; E[1][1] = "1"; E[1][2] = "0";
            E[2][0] = "0"; E[2][1] = "0"; E[2][2] = "1";

            r[0] = "0";  r[1] = q;  r[2] = "0";
            break;

        case XFORM_TZ:
            E[0][0] = "1"; E[0][1] = "0"; E[0][2] = "0";
            E[1][0] = "0"; E[1][1] = "1"; E[1][2] = "0";
            E[2][0] = "0"; E[2][1] = "0"; E[2][2] = "1";

            r[0] = "0";  r[1] = "0";  r[2] = q;
            break;

        default:
            break;
        }
    }

    // General fixed translation (no rotation): E = 1, r = (x, y, z), with
    // numeric components baked in at generation time. numStr keeps exact
    // zeros as "0" so the data's own sparsity folds through the symbolic
    // math just like structural sparsity does.
    Xform( std::string name, double x, double y, double z )
        :type( XFORM_3DOF_TRANS )
        ,baseName( name )
    {
        E[0][0] = "1"; E[0][1] = "0"; E[0][2] = "0";
        E[1][0] = "0"; E[1][1] = "1"; E[1][2] = "0";
        E[2][0] = "0"; E[2][1] = "0"; E[2][2] = "1";

        r[0] = numStr( x );  r[1] = numStr( y );  r[2] = numStr( z );
    }

    // Emits each entry as a named temp, then replaces the entry with the temp's
    // name so later compositions reference the emitted variable instead of
    // re-expanding the expression.  Literals 0 and 1 are kept as-is so that
    // sparsity keeps propagating through subsequent symbolic math.
    void emit() {
        // Generate E
        for ( int r = 0 ; r < 3 ; r++ )
        for ( int c = 0 ; c < 3 ; c++ )
        {
            prln( "double %s_E%d%d = %s;", baseName.c_str(), r, c, E[r][c].c_str() );
            if ( E[r][c] != "0" && E[r][c] != "1" ) E[r][c] = format( "{}_E{}{}", baseName, r, c );
        }
        for ( int i = 0 ; i < 3 ; i++ ) {
            prln( "double %s_r%d = %s;", baseName.c_str(), i, r[i].c_str() );
            if ( r[i] != "0" && r[i] != "1" ) r[i] = format( "{}_r{}", baseName, i );
        }
    }
};

Xform compose( std::string base, Xform &X_C_B, Xform &X_B_A ) {
    Xform X;
    X.baseName = base;

    // Generate E
    for ( int r = 0 ; r < 3 ; r++ )
    for ( int c = 0 ; c < 3 ; c++ )
    {
        string v = "0"s;
        for ( int i = 0 ; i < 3 ; i++ )
        {
            auto term = multFactors( X_C_B.E[r][i], X_B_A.E[i][c] );
            v = addTerms( v, term );
        }
        X.E[r][c] = v;
    }

    // Generate r = r_B_A + E_B_A^T * r_C_B
    for ( int r = 0 ; r < 3 ; r++ )
    {
        string v = X_B_A.r[r];
        for ( int i = 0 ; i < 3 ; i++ )
        {
            auto term = multFactors( X_B_A.E[i][r], X_C_B.r[i] );
            v = addTerms( v, term );
        }
        X.r[ r ] = v;
    }

    return X;
}

// Inverts a transform: if X maps A to B via (E, r), the inverse maps B to A
// via (E^T, -E*r). Derivation: w_A = E^T w_B and v_A = E^T v_B + r x (E^T w_B),
// which is the (E', r') form with E' = E^T, r' = -E*r since E (r x u) = (E r) x (E u)
// for a rotation E. Call this on an already-emitted transform so the E
// entries stay materialized names and only r' introduces new expressions.
Xform invertXform( string name, Xform &X ) {
    Xform inv;
    inv.baseName = name;

    for ( int r = 0 ; r < 3 ; r++ )
    for ( int c = 0 ; c < 3 ; c++ )
        inv.E[r][c] = X.E[c][r];

    for ( int k = 0 ; k < 3 ; k++ ) {
        string sum = "0";
        for ( int i = 0 ; i < 3 ; i++ )
            sum = addTerms( sum, multFactors( X.E[k][i], X.r[i] ) );
        inv.r[k] = negate( sum );
    }
    return inv;
}

// Transforms a motion vector from frame A to frame B, where X is the A-to-B transform.
// vB_angular = E * w
// vB_linear  = E * ( v - r x w )
Svec fwdXformMotion( string name, Xform &X, Svec &vA ) {
    Svec s( name );
    const string *w = &vA.v[ 0 ];
    const string *v = &vA.v[ 3 ];

    // Angular part: E * w
    for ( int r = 0 ; r < 3 ; r++ )
    {
        string sum = "0";
        for ( int i = 0 ; i < 3 ; i++ )
        {
            sum = addTerms( sum, multFactors( X.E[r][i], w[i] ) );
        }
        s.v[ r ] = sum;
    }

    // t = v - r x w
    string t[ 3 ];
    t[0] = subTerms( v[0], subTerms( multFactors( X.r[1], w[2] ), multFactors( X.r[2], w[1] ) ) );
    t[1] = subTerms( v[1], subTerms( multFactors( X.r[2], w[0] ), multFactors( X.r[0], w[2] ) ) );
    t[2] = subTerms( v[2], subTerms( multFactors( X.r[0], w[1] ), multFactors( X.r[1], w[0] ) ) );

    // Linear part: E * t
    for ( int r = 0 ; r < 3 ; r++ )
    {
        string sum = "0";
        for ( int i = 0 ; i < 3 ; i++ )
        {
            sum = addTerms( sum, multFactors( X.E[r][i], t[i] ) );
        }
        s.v[ r + 3 ] = sum;
    }

    return s;
}