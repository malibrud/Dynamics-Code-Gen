#pragma once
#include <vector>
#include <inttypes.h>
#include <assert.h>
#include "joint.h"
#include "xform.h"
#include "svec.h"

using namespace std;

struct System{
    int nBod;
    vector<int16_t> parent;
    vector<int16_t> children;
    vector<uint8_t> jtype;

    vector<Svec> S;
    vector<Svec> vJ;
    vector<Svec> v;
    vector<Svec> a;
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
        ,S( N )
        ,vJ( N )
        ,v( N )
        ,a( N )
        ,xTree( N )
        ,xJoint( N )
        ,xFromParent( N )
        ,xFromBase( N )
    {
    }
};
