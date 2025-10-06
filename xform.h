#pragma once
#include <math.h>
#include <assert.h>
#include <inttypes.h>
#include <stdbool.h>

#define XFORM_AXIS_R 0x20
#define XFORM_AXIS_T 0x40
#define XFORM_AXIS_TR ( XFORM_AXIS_T | XFORM_AXIS_R )

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
    double rx;
    double ry;
    double a;
    double cosa;
    double sina;
    double rxcosa;
    double rxsina;
    double rycosa;
    double rysina;
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
    Xca->rx    = d1 + d2 * cos1;
    Xca->ry    =  0 + d2 * sin1;
    Xca->cosa  =  cos1 * cos2 + sin1 * sin2;
    Xca->sina  =  cos1 * sin2 + sin1 * cos2;
    Xca->rxcosa  =  Xca->rx * Xca->cosa;
    Xca->rxsina  =  Xca->rx * Xca->sina;
    Xca->rycosa  =  Xca->ry * Xca->cosa;
    Xca->rysina  =  Xca->ry * Xca->sina;
}

void x_compose_pxy_pxy_to_pxy( Xform *_Xba, Xform *_Xcb, Xform *_Xca ) {
    assert( _Xba->type == XFORM_PXY );
    assert( _Xcb->type == XFORM_PXY );
    assert( _Xca->type == XFORM_PXY );

    Xform_AxisPlane *Xba = (Xform_AxisPlane*)_Xba;
    Xform_AxisPlane *Xcb = (Xform_AxisPlane*)_Xcb;
    Xform_AxisPlane *Xca = (Xform_AxisPlane*)_Xca;

    double x1 = Xba->rx;
    double y1 = Xba->ry;
    double x2 = Xcb->rx;
    double y2 = Xcb->ry;
    double a1 = Xba->a;
    double a2 = Xcb->a;
    double cos1 = Xba->cosa;
    double cos2 = Xcb->cosa;
    double sin1 = Xba->sina;
    double sin2 = Xcb->sina;

    Xca->a     = a1 + a2;
    Xca->rx    = x1 + x2*Xba->cosa - y2*Xba->sina;
    Xca->ry    = y1 + x2*Xba->sina + y2*Xba->cosa;
    Xca->cosa  =  cos1 * cos2 + sin1 * sin2;
    Xca->sina  =  cos1 * sin2 + sin1 * cos2;
    Xca->rxcosa  =  Xca->rx * Xca->cosa;
    Xca->rxsina  =  Xca->rx * Xca->sina;
    Xca->rycosa  =  Xca->ry * Xca->cosa;
    Xca->rysina  =  Xca->ry * Xca->sina;
}

// Composition of transforms
//   Xca = Xcb * Xba
// See (Featherstone, 2008) Sec 2.8, pp 20. 
void x_compose( Xform *Xba, Xform *Xcb, Xform *Xca ) {
    if ( Xba->type == XFORM_TX && Xcb->type == XFORM_RZ ) {
        assert( Xca->type == XFORM_TX_RZ );
        x_compose_tx_rz( Xba, Xcb, Xca );
        return;
    } else
    if ( Xba->type == XFORM_TX_RZ && Xcb->type == XFORM_TX_RZ ) {
        assert( Xca->type == XFORM_PXY );
        x_compose_txrz_txrz_to_pxy( Xba, Xcb, Xca );
        return;
    }
}