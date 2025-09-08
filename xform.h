#pragma once
#include <inttypes.h>

#define XFORM_IDENT 0

#define XFORM_GEN 1

#define XFORM_ROT 2

#define XFORM_ROT_X 3
#define XFORM_ROT_Y 4
#define XFORM_ROT_Z 5

#define XFORM_TRANS_X 6
#define XFORM_TRANS_Y 7
#define XFORM_TRANS_Z 8

#define XFORM_DH_LINK 9
#define XFORM_DH_JOINT 10

typedef struct {
    int8_t type;


} Xform;