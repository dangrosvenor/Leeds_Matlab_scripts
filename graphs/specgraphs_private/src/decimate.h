/* Copyright 1993-2002 The MathWorks, Inc.  */

/*
  Source file for reducep MEX file
*/

/* $Revision: 1.4 $ */

#include "quadrics.h"

extern Vertex *decimate_last_v0;
extern Vertex *decimate_last_v1;

extern bool decimate_quadric(Vertex *v, Mat4& Q);
extern real decimate_min_error();
extern real decimate_max_error(Model& m);
extern real decimate_error(Vertex *);
extern void decimate_contract(Model& m);
extern void decimate_init(Model& m, real limit);
extern void decimate_cleanup();
