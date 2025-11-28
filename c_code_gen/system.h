#pragma once
#include <vector>
#include <inttypes.h>
#include <assert.h>
#include "joint.h"
#include "xform.h"

using namespace std;

struct System{
    int nBod;
    vector<int16_t> parent;
    vector<int16_t> children;
    vector<uint8_t> jtype;

    vector<Xform> xTree;
    vector<Xform> xJoint;
    vector<Xform> xFromParent;
    vector<Xform> xFromBase;

    System(): System(0) {}
    System( int N )
        :nBod( N )
        ,parent( N )
        ,children( N )
        ,jtype( N )
        ,xTree( N )
        ,xJoint( N )
        ,xFromParent( N )
        ,xFromBase( N )
    {
    }
};
