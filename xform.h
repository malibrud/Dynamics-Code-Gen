#pragma once
#include <math.h>
#include <assert.h>
#include <inttypes.h>
#include <stdbool.h>

#define XFORM_AXIS_R 0x20
#define XFORM_AXIS_T 0x40
#define XFORM_AXIS_TR 0x40

#define XFORM_IDENT 0
#define XFORM_6DOF 1
#define XFORM_3DOF_ROT 2
#define XFORM_3DOF_TRANS 3

// Transformations around the X, Y, and Z axes.
#define XFORM_RX ( XFORM_AXIS_R | 4 )
#define XFORM_RY ( XFORM_AXIS_R | 5 )
#define XFORM_RZ ( XFORM_AXIS_R | 6 )

// Transformations along the X, Y, and Z axes.
#define XFORM_TX ( XFORM_AXIS_T | 8 )
#define XFORM_TY ( XFORM_AXIS_T | 9 )
#define XFORM_TZ ( XFORM_AXIS_T | 10 )

// Transformations with translation along an axis then rotation around an axis
#define XFORM_TX_RX ( XFORM_AXIS_TR | 12 )
#define XFORM_TX_RY ( XFORM_AXIS_TR | 13 )
#define XFORM_TX_RZ ( XFORM_AXIS_TR | 14 )

#define XFORM_TY_RX ( XFORM_AXIS_TR | 16 )
#define XFORM_TY_RY ( XFORM_AXIS_TR | 17 )
#define XFORM_TY_RZ ( XFORM_AXIS_TR | 18 )

#define XFORM_TZ_RX ( XFORM_AXIS_TR | 21 )
#define XFORM_TZ_RY ( XFORM_AXIS_TR | 22 )
#define XFORM_TZ_RZ ( XFORM_AXIS_TR | 23 )

// Transformations in the planes XY, YZ, ZX
#define XFORM_PXY 12
#define XFORM_PYZ 13
#define XFORM_PZX 14

// Transformations defined by D-H parameters
#define XFORM_DH_LINK  16
#define XFORM_DH_JOINT 17

typedef struct {
    int8_t type;
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
    double d;
} Xform_AxisPlaneXY;

typedef struct {
    Xform X;
    double d;
    double a;
    double dsina;
    double dcosa;
} Xform_XtransZrot;

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

void x_compose_tx_rz( Xform *_Xba, Xform *_Xcb, Xform *_Xca ) {
    assert( _Xba->type == XFORM_TX );
    assert( _Xcb->type == XFORM_RZ );
    assert( _Xca->type == XFORM_TX_RZ );

    Xform_AxisTrans  *Xba = (Xform_AxisTrans*)_Xba;
    Xform_AxisRot    *Xcb = (Xform_AxisRot*)_Xcb;
    Xform_XtransZrot *Xca = (Xform_XtransZrot*)_Xca;

    Xca->d = Xba->d;
    Xca->a = Xcb->a;
    Xca->dsina = Xba->d * Xcb->sina;
    Xca->dcosa = Xba->d * Xcb->cosa;
}

// Composition of transforms
//   Xca = Xcb * Xba
// See (Featherstone, 2008) Sec 2.8, pp 20. 
void x_compose( Xform *Xba, Xform *Xcb, Xform *Xca ) {
    if ( Xba->type == XFORM_TX && Xcb->type == XFORM_RZ ) {
        assert( Xca->type == XFORM_TX_RZ );
        x_compose_tx_rz( Xba, Xcb, Xca );
        return;
    }
}