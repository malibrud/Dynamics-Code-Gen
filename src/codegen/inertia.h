#pragma once
#include <string>
#include <format>
#include <array>
#include "vec3.h"
#include "xform.h"
#include "svec.h"
#include "symmath.h"
#include "prln.h"

using namespace std;

// Physical description of a single rigid body, expressed about its own
// body-fixed frame origin. These are ordinary numbers, known at
// code-generation time -- only the joint variables q, qd are symbolic/runtime.
typedef struct {
    double xx, yy, zz, xy, yz, zx; // about the body's own CG, off-diagonal signs absorbed
} Inertia;

typedef struct {
    double m;
    double cg[3]; // CG location relative to the body's own frame origin
    Inertia I;    // rotational inertia about the CG
} Body;

// A spatial (6x6, symmetric) rigid-body inertia in the compact 10-parameter
// form (RBDA sec 2.13), with each entry a symbolic expression string. Field
// layout mirrors src/numeric/inertia.h's SpatialInertia so the two stay easy
// to compare.
struct SpatialInertia {
    string baseName;
    string m;
    string sx, sy, sz;                   // first moment of mass, h = m*cg
    string Ixx, Iyy, Izz, Ixy, Iyz, Izx;  // rotational inertia about the frame origin

    SpatialInertia() : SpatialInertia( ""s ) {}
    SpatialInertia( string name )
        :baseName( name )
        ,m("0"), sx("0"), sy("0"), sz("0")
        ,Ixx("0"), Iyy("0"), Izz("0"), Ixy("0"), Iyz("0"), Izx("0")
    {}

    // Emits each field as a named temp, then replaces the field with the
    // temp's name (0/1 literals are left alone so sparsity keeps propagating).
    void emit() {
        auto emitField = [&]( const char *suffix, string &field ) {
            prln( "double %s_%s = %s;", baseName.c_str(), suffix, field.c_str() );
            if ( field != "0" && field != "1" ) field = format( "{}_{}", baseName, suffix );
        };
        emitField( "m",   m );
        emitField( "sx",  sx );
        emitField( "sy",  sy );
        emitField( "sz",  sz );
        emitField( "Ixx", Ixx );
        emitField( "Iyy", Iyy );
        emitField( "Izz", Izz );
        emitField( "Ixy", Ixy );
        emitField( "Iyz", Iyz );
        emitField( "Izx", Izx );
    }
};

// Builds the spatial inertia of a body about its own frame origin, given its
// mass, rotational inertia about its CG, and the CG offset -- i.e. the
// parallel axis theorem. All inputs are plain numbers (known at generation
// time), so the arithmetic happens now in C++ and the result is baked in as
// numeric-literal expression strings; only later composition with joint
// transforms introduces q.
SpatialInertia bodySpatialInertia( string name, Body &B ) {
    SpatialInertia I( name );
    double x = B.cg[0], y = B.cg[1], z = B.cg[2];
    double m = B.m;

    I.m   = format( "{}", m );
    I.sx  = format( "{}", m*x );
    I.sy  = format( "{}", m*y );
    I.sz  = format( "{}", m*z );
    I.Ixx = format( "{}", B.I.xx + m*( z*z + y*y ) );
    I.Iyy = format( "{}", B.I.yy + m*( x*x + z*z ) );
    I.Izz = format( "{}", B.I.zz + m*( y*y + x*x ) );
    I.Ixy = format( "{}", B.I.xy - m*x*y );
    I.Iyz = format( "{}", B.I.yz - m*y*z );
    I.Izx = format( "{}", B.I.zx - m*z*x );

    return I;
}

SpatialInertia addInertia( string name, SpatialInertia &a, SpatialInertia &b ) {
    SpatialInertia s( name );
    s.m   = addTerms( a.m,   b.m );
    s.sx  = addTerms( a.sx,  b.sx );
    s.sy  = addTerms( a.sy,  b.sy );
    s.sz  = addTerms( a.sz,  b.sz );
    s.Ixx = addTerms( a.Ixx, b.Ixx );
    s.Iyy = addTerms( a.Iyy, b.Iyy );
    s.Izz = addTerms( a.Izz, b.Izz );
    s.Ixy = addTerms( a.Ixy, b.Ixy );
    s.Iyz = addTerms( a.Iyz, b.Iyz );
    s.Izx = addTerms( a.Izx, b.Izx );
    return s;
}

// Expands the compact 10-parameter form into the full symmetric 6x6 block
// matrix
//     [  Ibar    h x  ]
//     [ (h x)^T  m I  ]
// (RBDA eq 2.63; h x is the skew matrix of h = (sx,sy,sz), and the lower-left
// block is its transpose, i.e. -(h x)), so it can be pushed through generic
// matrix algebra. Sanity anchor: for pure rotation w = z about the origin the
// linear momentum must be m*(w x c) = (-sy, sx, 0) -- that's the lower-left
// block's third column. The lower-left originally had the sign flipped
// (copied from src/numeric's si_multBySvec, which shares the same bug); the
// dense validator in src/validate caught it.
array<array<string,6>,6> spatialInertiaTo6x6( SpatialInertia &I ) {
    array<array<string,6>,6> M;
    string nsx = subTerms( "0", I.sx );
    string nsy = subTerms( "0", I.sy );
    string nsz = subTerms( "0", I.sz );

    M[0] = { I.Ixx, I.Ixy, I.Izx, "0",  nsz,  I.sy };
    M[1] = { I.Ixy, I.Iyy, I.Iyz, I.sz, "0",  nsx  };
    M[2] = { I.Izx, I.Iyz, I.Izz, nsy,  I.sx, "0"  };
    M[3] = { "0",   I.sz,  nsy,   I.m,  "0",  "0"  };
    M[4] = { nsz,   "0",   I.sx,  "0",  I.m,  "0"  };
    M[5] = { I.sy,  nsx,   "0",   "0",  "0",  I.m  };
    return M;
}

// Packs a (mathematically symmetric) full 6x6 block matrix back into the
// compact 10-parameter form, reading only the upper triangle.
SpatialInertia spatialInertiaFrom6x6( string name, array<array<string,6>,6> &M ) {
    SpatialInertia I( name );
    I.Ixx = M[0][0]; I.Iyy = M[1][1]; I.Izz = M[2][2];
    I.Ixy = M[0][1]; I.Iyz = M[1][2]; I.Izx = M[0][2];
    I.sx  = M[2][4]; I.sy  = M[0][5]; I.sz  = M[1][3];
    I.m   = M[3][3];
    return I;
}

// F = I * v (RBDA's "m2f" operation), computed generically from the full 6x6
// form. For our single-DOF joints v = S is hyper-sparse ([0,0,1,0,0,0] for a
// z-axis revolute), so multFactors/addTerms collapse this down to picking out
// one column of M. (Matches src/numeric's bd_m2f; note src/numeric's
// si_multBySvec had the linear-part signs flipped -- see spatialInertiaTo6x6.)
Svec spatialInertiaMulVec( string name, SpatialInertia &I, Svec &v ) {
    auto M = spatialInertiaTo6x6( I );
    Svec out( name );
    for ( int r = 0 ; r < 6 ; r++ ) {
        string sum = "0";
        for ( int c = 0 ; c < 6 ; c++ ) sum = addTerms( sum, multFactors( M[r][c], v.v[c] ) );
        out.v[r] = sum;
    }
    return out;
}

// Expands a motion transform X into its full 6x6 matrix form by pushing each
// of the 6 standard basis vectors through the already-validated
// fwdXformMotion -- this guarantees the matrix form agrees with the vector
// form used everywhere else, instead of re-deriving the block structure
// [[E,0],[-E*rx,E]] by hand.
array<array<string,6>,6> xformTo6x6( Xform &X ) {
    array<array<string,6>,6> M;
    for ( int col = 0 ; col < 6 ; col++ ) {
        Svec e( ""s );
        e.v[ col ] = "1";
        Svec c = fwdXformMotion( "tmp"s, X, e );
        for ( int row = 0 ; row < 6 ; row++ ) M[row][col] = c.v[row];
    }
    return M;
}

array<array<string,6>,6> smat6Transpose( array<array<string,6>,6> &A ) {
    array<array<string,6>,6> T;
    for ( int r = 0 ; r < 6 ; r++ )
    for ( int c = 0 ; c < 6 ; c++ )
        T[r][c] = A[c][r];
    return T;
}

array<array<string,6>,6> smat6Mul( array<array<string,6>,6> &A, array<array<string,6>,6> &B ) {
    array<array<string,6>,6> C;
    for ( int r = 0 ; r < 6 ; r++ )
    for ( int c = 0 ; c < 6 ; c++ ) {
        string sum = "0";
        for ( int k = 0 ; k < 6 ; k++ ) sum = addTerms( sum, multFactors( A[r][k], B[k][c] ) );
        C[r][c] = sum;
    }
    return C;
}

Svec smat6MulVec( array<array<string,6>,6> &A, Svec &v ) {
    Svec out( ""s );
    for ( int r = 0 ; r < 6 ; r++ ) {
        string sum = "0";
        for ( int c = 0 ; c < 6 ; c++ ) sum = addTerms( sum, multFactors( A[r][c], v.v[c] ) );
        out.v[r] = sum;
    }
    return out;
}

// Congruence transform of a composite spatial inertia across a joint, i.e.
// the CRBA recursion  Ic_parent += X^T * Ic_child * X  (RBDA Table 6.2/6.4).
// This is exactly what spatial_v2's HandC.m does with plain 6x6 matrices --
// IC{p} = IC{p} + Xup{i}'*IC{i}*Xup{i} -- just done symbolically here. X is
// child_X_parent (our usual xFromParent[i] convention, matching Xup there).
SpatialInertia transformInertia( string name, Xform &X, SpatialInertia &Ic ) {
    auto X6  = xformTo6x6( X );
    auto Xt6 = smat6Transpose( X6 );
    auto M6  = spatialInertiaTo6x6( Ic );
    auto N   = smat6Mul( M6, X6 );   // Ic * X
    auto R   = smat6Mul( Xt6, N );   // X^T * Ic * X
    return spatialInertiaFrom6x6( name, R );
}
