/*
 *  VTK_Out.h
 *  Read UGRID
 *
 *  Created by Matthew O'Connell on 6/2/10.
 *  Copyright 2010  . All rights reserved.
 *
 */


#ifndef VTK_Grid_Out_h
#define VTK_Grid_Out_h
int VTK_Grid_Out(char* name, int nnodes, double *xyz, int *nelem, int **c2n);

int VTK_Grid_Out_with_solution(char* name, int nnodes, double *xyz, int *nelem, int **c2n, double *solution, int nsolution_variables/*, int **face_tag*/);

#endif