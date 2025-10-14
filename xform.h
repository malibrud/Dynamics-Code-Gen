#pragma once
#include <math.h>
#include <assert.h>
#include <inttypes.h>
#include <stdbool.h>
#include "arena.h"

Arena *x_arena = NULL;
void x_init( Arena *a ) {
    assert( a );
    x_arena = a;
}

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
#define XFORM_PYZ ( XFORM_AXIS_P | (Xtype)25 )
#define XFORM_PZX ( XFORM_AXIS_P | (Xtype)27 )

// Transformations defined by D-H parameters
#define XFORM_DH_LINK  (Xtype)29
#define XFORM_DH_JOINT (Xtype)30

typedef struct {
    Xtype type;
} Xform;

typedef struct {
    Xform X;
    double a;
    double sina;
    double cosa;
} Xform_AxisRot;

typedef struct {
    Xform X;
    double d;
} Xform_AxisTrans;

typedef struct {
    Xform X;
    double c1; // Coordinate 1 ...  if this is   x, y, or z
    double c2; // Coordinate 2 ...  then this is y, z, or x
    double a;
    double cosa;
    double sina;
    double c1cosa;
    double c1sina;
    double c2cosa;
    double c2sina;
} Xform_AxisPlane;

typedef struct {
    Xform X;
    double d;
    double a;
    double sina;
    double cosa;
    double dsina;
    double dcosa;
} Xform_AxisTransRot;

// Cannonical transforms
Xform x_identity = { XFORM_IDENT };

void x_update( Xform *X, double *q ) {
    switch ( X->type )
    {
    case XFORM_IDENT:
        return;
    case XFORM_RX:
    case XFORM_RY:
    case XFORM_RZ:
        Xform_AxisRot *Xar = (Xform_AxisRot*)X;
        Xar->a =  *q;
        Xar->sina = sin( *q );
        Xar->cosa = cos( *q );
        return;
    default:
        assert(false);
        break;
    }
}

Xform *x_newAxisTrans( Xtype type, double d ) {
    assert( type & XFORM_AXIS_T );

    Xform_AxisTrans *X;
    X = ar_alloc_aligned( x_arena, sizeof( Xform_AxisTrans ), ALIGNMENT );
    X->X.type = type;
    X->d = d;
    return (Xform*)X;
}

Xform *x_newAxisRot( Xtype type, double a ) {
    assert( type & XFORM_AXIS_R );

    Xform_AxisRot *X;
    X = ar_alloc_aligned( x_arena, sizeof( Xform_AxisRot ), ALIGNMENT );
    X->X.type = type;
    X->a = a;
    X->cosa = cos( a );
    X->sina = sin( a );
    return (Xform*)X;
}

Xform *x_newAxisTransRot( Xtype type, double d, double a ) {
    assert( type & XFORM_AXIS_T );
    assert( type & XFORM_AXIS_R );

    Xform_AxisTransRot *X;
    X = ar_alloc_aligned( x_arena, sizeof( Xform_AxisTransRot ), ALIGNMENT );
    X->X.type = type;
    X->a = a;
    X->d = d;
    X->cosa = cos( a );
    X->sina = sin( a );
    X->dcosa = d*X->cosa;
    X->dsina = d*X->sina;
    return (Xform*)X;
}

Xform *x_newAxisPlane( Xtype type, double c1, double c2, double a ) {
    assert( type & XFORM_AXIS_P );

    Xform_AxisPlane *X;
    X = ar_alloc_aligned( x_arena, sizeof( Xform_AxisPlane ), ALIGNMENT );
    X->X.type = type;
    X->a = a;
    X->c1 = c1;
    X->c2 = c2;
    double c = cos( a );
    double s = sin( a );
    X->cosa = c;
    X->sina = s;
    X->c1cosa = c1*c;
    X->c1sina = c1*s;
    X->c2cosa = c2*c;
    X->c2sina = c2*s;
    return (Xform*)X;
}

Xform *x_newTx( double d ) {
    return x_newAxisTrans( XFORM_TX, d );
}

Xform *x_newRz( double a ) {
    return x_newAxisRot( XFORM_RZ, a );
}

Xform *x_newTxRz( double d, double a ) {
    return x_newAxisTransRot( XFORM_TX_RZ, d, a );
}

Xform *x_newPxy( double x, double y, double a ) {
    return x_newAxisPlane( XFORM_PXY, x, y, a );
}

// Xform *x_new( Xtype type ) {
//     Xform *X;
//     switch ( type )
//     {
//     case XFORM_IDENT:
//         X = ar_alloc_aligned( x_arena, sizeof( Xform ), sizeof( void* ) );
//         X->type = type;
//         return X;
//     
//     default:
//         break;
//     }
// 
// }

void x_compose_tx_rz_to_txrz( Xform *_Xba, Xform *_Xcb, Xform *_Xca ) {
    assert( _Xba->type == XFORM_TX );
    assert( _Xcb->type == XFORM_RZ );
    assert( _Xca->type == XFORM_TX_RZ );

    Xform_AxisTrans    *Xba = (Xform_AxisTrans*)_Xba;
    Xform_AxisRot      *Xcb = (Xform_AxisRot*)_Xcb;
    Xform_AxisTransRot *Xca = (Xform_AxisTransRot*)_Xca;

    Xca->d     = Xba->d;
    Xca->a     = Xcb->a;
    Xca->sina  = Xcb->sina;
    Xca->cosa  = Xcb->cosa;
    Xca->dsina = Xba->d * Xcb->sina;
    Xca->dcosa = Xba->d * Xcb->cosa;
}

void x_compose_txrz_txrz_to_pxy( Xform *_Xba, Xform *_Xcb, Xform *_Xca ) {
    assert( _Xba->type == XFORM_TX_RZ );
    assert( _Xcb->type == XFORM_TX_RZ );
    assert( _Xca->type == XFORM_PXY );

    Xform_AxisTransRot *Xba = (Xform_AxisTransRot*)_Xba;
    Xform_AxisTransRot *Xcb = (Xform_AxisTransRot*)_Xcb;
    Xform_AxisPlane *Xca = (Xform_AxisPlane*)_Xca;

    double d1 = Xba->d;
    double d2 = Xcb->d;
    double a1 = Xba->a;
    double a2 = Xcb->a;
    double cos1 = Xba->cosa;
    double cos2 = Xcb->cosa;
    double sin1 = Xba->sina;
    double sin2 = Xcb->sina;

    Xca->a     = a1 + a2;
    Xca->c1    = d1 + d2 * cos1;
    Xca->c2    =  0 + d2 * sin1;
    Xca->cosa  =  cos1 * cos2 + sin1 * sin2;
    Xca->sina  =  cos1 * sin2 + sin1 * cos2;
    Xca->c1cosa  =  Xca->c1 * Xca->cosa;
    Xca->c1sina  =  Xca->c1 * Xca->sina;
    Xca->c2cosa  =  Xca->c2 * Xca->cosa;
    Xca->c2sina  =  Xca->c2 * Xca->sina;
}

void x_compose_pxy_pxy_to_pxy( Xform *_Xba, Xform *_Xcb, Xform *_Xca ) {
    assert( _Xba->type == XFORM_PXY );
    assert( _Xcb->type == XFORM_PXY );
    assert( _Xca->type == XFORM_PXY );

    Xform_AxisPlane *Xba = (Xform_AxisPlane*)_Xba;
    Xform_AxisPlane *Xcb = (Xform_AxisPlane*)_Xcb;
    Xform_AxisPlane *Xca = (Xform_AxisPlane*)_Xca;

    double x1 = Xba->c1;
    double y1 = Xba->c2;
    double x2 = Xcb->c1;
    double y2 = Xcb->c2;
    double a1 = Xba->a;
    double a2 = Xcb->a;
    double cos1 = Xba->cosa;
    double cos2 = Xcb->cosa;
    double sin1 = Xba->sina;
    double sin2 = Xcb->sina;

    Xca->a     = a1 + a2;
    Xca->c1    = x1 + x2*Xba->cosa - y2*Xba->sina;
    Xca->c2    = y1 + x2*Xba->sina + y2*Xba->cosa;
    Xca->cosa  =  cos1 * cos2 + sin1 * sin2;
    Xca->sina  =  cos1 * sin2 + sin1 * cos2;
    Xca->c1cosa  =  Xca->c1 * Xca->cosa;
    Xca->c1sina  =  Xca->c1 * Xca->sina;
    Xca->c2cosa  =  Xca->c2 * Xca->cosa;
    Xca->c2sina  =  Xca->c2 * Xca->sina;
}

void x_copy( Xform *Xsrc, Xform *Xdst ) {
    if ( Xsrc->type == XFORM_TX_RZ && Xdst->type == XFORM_PXY ) {
        Xform_AxisTransRot *src = (Xform_AxisTransRot*)Xsrc;
        Xform_AxisPlane    *dst = (Xform_AxisPlane*)Xdst;
        dst->a      = src->a;
        dst->cosa   = src->cosa;
        dst->sina   = src->sina;
        dst->c1     = src->d;
        dst->c2     = 0;
        dst->c1cosa = src->dcosa;
        dst->c1sina = src->dsina;
        dst->c2cosa = 0;
        dst->c2sina = 0;
    } else {
        assert( false );
    }
}

// Composition of transforms
//   Xca = Xcb * Xba
// See RBDA Sec 2.8, pp 20. 
void x_compose( Xform *Xba, Xform *Xcb, Xform *Xca ) {
    if ( Xba->type == XFORM_TX && Xcb->type == XFORM_RZ ) {
        assert( Xca->type == XFORM_TX_RZ );
        x_compose_tx_rz_to_txrz( Xba, Xcb, Xca );
        return;
    } else
    if ( Xba->type == XFORM_TX_RZ && Xcb->type == XFORM_TX_RZ ) {
        assert( Xca->type == XFORM_PXY );
        x_compose_txrz_txrz_to_pxy( Xba, Xcb, Xca );
        return;
    } else
    if ( Xba->type == XFORM_PXY && Xcb->type == XFORM_PXY ) {
        assert( Xca->type == XFORM_PXY );
        x_compose_pxy_pxy_to_pxy( Xba, Xcb, Xca );
        return;
    } else
    if ( Xba->type == XFORM_IDENT ) {
        x_copy( Xcb, Xca );
        return;
    } else {
        assert( false );
    }
}

void x_mapply_TxRz( Xform *xIn, Svec *vIn, Svec *vOut ) {
    assert( xIn->type == XFORM_TX_RZ );
    Xform_AxisTransRot *X = (Xform_AxisTransRot*)xIn;
    double s   = X->sina;
    double c   = X->cosa;
    double rx  = X->d;
    double rxs = X->dsina;
    double rxc = X->dcosa;

    double *vi = vIn->v;
    double *vo = vOut->v;

    vo[ 0 ] =  c*vi[0] + s*vi[1] +  0;
    vo[ 1 ] = -s*vi[0] + c*vi[1] +  0;
    vo[ 2 ] =    0     +   0     + vi[2];
    
    vo[ 3 ] = 0 +  0       + rxs*vi[2] + c*vi[3] + s*vi[4] + 0;
    vo[ 4 ] = 0 +  0       + rxc*vi[2] - s*vi[3] + c*vi[4] + 0;
    vo[ 5 ] = 0 - rx*vi[1] +     0     +   0     +    0    + vi[5];
}

void x_fapply_TxRz( Xform *xIn, Svec *vIn, Svec *vOut ) {
    assert( xIn->type == XFORM_TX_RZ );
    Xform_AxisTransRot *X = (Xform_AxisTransRot*)xIn;
    double s   = X->sina;
    double c   = X->cosa;
    double rx  = X->d;
    double rxs = X->dsina;
    double rxc = X->dcosa;

    double *vi = vIn->v;
    double *vo = vOut->v;

    vo[ 0 ] =  c*vi[0] + s*vi[1] +  0;
    vo[ 1 ] = -s*vi[0] + c*vi[1] +  0;
    vo[ 2 ] =    0     +   0     + vi[2];
    
    vo[ 3 ] = 0 +  0       + rxs*vi[2] + c*vi[3] + s*vi[4] + 0;
    vo[ 4 ] = 0 +  0       + rxc*vi[2] - s*vi[3] + c*vi[4] + 0;
    vo[ 5 ] = 0 - rx*vi[1] +     0     +   0     +    0    + vi[5];
}

// 
void x_finvapply_TxRz( Xform *xIn, Svec *vIn, Svec *vOut ) {
    assert( xIn->type == XFORM_TX_RZ );
    Xform_AxisTransRot *X = (Xform_AxisTransRot*)xIn;
    double s   = X->sina;
    double c   = X->cosa;
    double rx  = X->d;
    double rxs = X->dsina;
    double rxc = X->dcosa;

    double *vi = vIn->v;
    double *vo = vOut->v;

    vo[ 0 ] = c*vi[0] - s*vi[1] +  0    +    0      +    0      + 0;
    vo[ 1 ] = s*vi[0] + c*vi[1] +  0    +    0      +    0      - rx;
    vo[ 2 ] =   0     +   0     + vi[2] + rxs*vi[3] + rxc*vi[4] + 0;
    
    vo[ 3 ] =   0     +   0     +  0    +  c*vi[3]  -  s*vi[4]  + 0;
    vo[ 4 ] =   0     +   0     +  0    +  s*vi[3]  +  c*vi[4]  + 0;
    vo[ 5 ] =   0     +   0     +  0    +    0      +     0     + vi[5];
}

// Apply the transformation to a motion vector, i.e. velocity or acceleration
// See RBDA eq (2.24) pp 22
void x_mapply( Xform *X, Svec *vMult, Svec *vResult ) {
    switch ( X->type )
    {
    case XFORM_TX_RZ:
        x_mapply_TxRz( X, vMult, vResult );
        break;
    
    default:
        assert( false );
        break;
    }
}

// Apply the transformation to a force vector, i.e. force or momentum
// See RBDA eq (2.25) pp 22
void x_fapply( Xform *X, Svec *vMult, Svec *vResult ) {
    switch ( X->type )
    {
    case XFORM_TX_RZ:
        x_fapply_TxRz( X, vMult, vResult );
        break;
    
    default:
        assert( false );
        break;
    }
}

// Apply the inverse transformation to a force vector, i.e. force or momentum
// See RBDA eq (2.27) pp 22
void x_finvapply( Xform *X, Svec *vMult, Svec *vResult ) {
    switch ( X->type )
    {
    case XFORM_TX_RZ:
        x_finvapply_TxRz( X, vMult, vResult );
        break;
    
    default:
        assert( false );
        break;
    }
}