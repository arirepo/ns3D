/*
 *  VTK_Out.c
 *  Read UGRID
 *
 *  Created by Matthew O'Connell on 6/2/10.
 *  Copyright 2010  . All rights reserved.
 *
 */

#include "VTK_Out.h"
#include "generic.h"

int VTK_Grid_Out_with_solution(char* name, int nnodes, double *xyz, int *nelem, int **c2n, double *solution, int nsolution_variables/*, int **face_tag*/){
    
    
    int vtk_etype[7] = {-1,5,9,10,14,13,12};
    
    printf("\nPARAVIEW VTK OUTPUT");
    FILE *fp;
    char filename[32];
    char suffix[5] = ".vtk";
    sprintf(filename,"%s",name);
    strncat(filename,suffix,4);
    //strcat(filename,".vtk");
    printf("\nVTK Filename = <%s>\n",filename);
    // Open file for write
    if ((fp = fopen(filename,"w")) == NULL)
    {
        printf("\nError opening file <%s>.",filename);
        return 0;
    }
    
    //Print Header
    fprintf(fp, "# vtk DataFile Version 3.0");
    fprintf(fp, "\n%s",name);
    fprintf(fp, "\nASCII");
    fprintf(fp, "\nDATASET UNSTRUCTURED_GRID");
    
    //Print Points
    fprintf(fp, "\nPOINTS %d double",nnodes); // mod
    int i;
    for (i = 1; i <= nnodes; i++) {
        fprintf(fp, "\n%.15e %.15e %.15e",xyz[3*i + 0], xyz[3*i + 1], xyz[3*i+2]);
    }
    
    //Print Cells
    int total_cells = 0;
    int etype, mnode;
    int size = 0;
    
    for (etype = Triangle; etype <= Hex; etype++) {
        mnode = ElemTypeToMnode(etype);
        total_cells += nelem[etype];
        size += nelem[etype]*mnode;
    }
    size += total_cells;
    fprintf(fp, "\n\nCELLS %d %d", total_cells, size);
    
    int j;
    for (etype = Triangle; etype <= Hex; etype++) {
        mnode = ElemTypeToMnode(etype);
        for (i = 1; i <= nelem[etype]; i++) {
            fprintf(fp, "\n%d ",mnode);
            for (j = 0; j < mnode; j++) {
                fprintf(fp, "%d ",c2n[etype][mnode*i + j]-1);
            }
        }
    }
    
    fprintf(fp, "\n");
    
    fprintf(fp, "\nCELL_TYPES %d",total_cells);
    
    for (etype = Triangle ; etype <= Hex; etype++) {
        for (i = 0; i < nelem[etype]; i++) {
            fprintf(fp, "\n%d",vtk_etype[etype]);
        }
    }
    
    fprintf(fp, "\n\nPOINT_DATA %d",nnodes);
    fprintf(fp, "\nSCALARS density double 1"); // mod
    fprintf(fp, "\nLOOKUP_TABLE default");
    for (i = 1; i < nnodes+1; i++) {
        fprintf(fp, "\n%.15e",solution[5*i + 0]);
    }
 
    fprintf(fp, "\nSCALARS x-momentum double 1"); // mod
    fprintf(fp, "\nLOOKUP_TABLE default");
    for (i = 1; i < nnodes+1; i++) {
        fprintf(fp, "\n%.15e",solution[5*i + 1]);
    }  
    
    fprintf(fp, "\nSCALARS y-momentum double 1"); // mod
    fprintf(fp, "\nLOOKUP_TABLE default");
    for (i = 1; i < nnodes+1; i++) {
        fprintf(fp, "\n%.15e",solution[5*i + 2]);
    }  
    
    fprintf(fp, "\nSCALARS z-momentum double 1"); // mod
    fprintf(fp, "\nLOOKUP_TABLE default");
    for (i = 1; i < nnodes+1; i++) {
        fprintf(fp, "\n%.15e",solution[5*i + 3]);
    }  
    
    fprintf(fp, "\nSCALARS energytotal_density double 1"); // mod
    fprintf(fp, "\nLOOKUP_TABLE default");
    for (i = 1; i < nnodes+1; i++) {
        fprintf(fp, "\n%.15e",solution[5*i + 4]);
    }  
    
    int nelements = 0;
    for (etype = Triangle; etype <= Hex; etype++) {
        nelements += nelem[etype];
    }
    
    /*
    fprintf(fp, "\nCELL DATA %d",nelements + 3 + 1);
    fprintf(fp, "\nSCALARS surfacetag FLOAT");
    fprintf(fp, "\nLOOKUP_TABLE default");
    int e;
    for (etype = Triangle; etype <= Quad; etype++) {
        for (e = 1; e <= nelem[etype]; e++) {
            fprintf(fp, "\n%d",face_tag[etype][e]);
        }
    }
    for (; etype <= Hex; etype++) {
        for (e = 1; e <= nelem[etype]; etype++) {
            fprintf(fp, "\n0");
        }
    }
    */
    
    fprintf(fp, "\n");
    fclose(fp);
    return 1;
}

int VTK_Grid_Out(char* name, int nnodes, double *xyz, int *nelem, int **c2n){
    
    
    int vtk_etype[7] = {-1,5,9,10,14,13,12};
    
    printf("\nPARAVIEW VTK OUTPUT");
    FILE *fp;
    char filename[32];
    filename[0]='\0';
    strcat(filename, name);
    strcat(filename,".vtk");
    printf("\nVTK Filename = <%s>\n",filename);
    // Open file for write
    if ((fp = fopen(filename,"w")) == NULL)
    {
        printf("\nError opening file <%s>.",filename);
        return 0;
    }
    
    //Print Header
    fprintf(fp, "# vtk DataFile Version 3.0");
    fprintf(fp, "\n%s",name);
    fprintf(fp, "\nASCII");
    fprintf(fp, "\nDATASET UNSTRUCTURED_GRID");
    
    //Print Points
    fprintf(fp, "\nPOINTS %d double",nnodes);
    int i;
    for (i = 1; i <= nnodes; i++) {
        fprintf(fp, "\n%.15e %.15e %.15e",xyz[3*i + 0], xyz[3*i + 1], xyz[3*i+2]);
    }

    //Print Cells
    int total_cells = 0;
    int etype, mnode;
    int size = 0;
    
    for (etype = Triangle; etype <= Hex; etype++) {
        mnode = ElemTypeToMnode(etype);
        total_cells += nelem[etype];
        size += nelem[etype]*mnode;
    }
    size += total_cells;
    fprintf(fp, "\n\nCELLS %d %d", total_cells, size);
    
    int j;
    for (etype = Triangle; etype <= Hex; etype++) {
        mnode = ElemTypeToMnode(etype);
        for (i = 1; i <= nelem[etype]; i++) {
            fprintf(fp, "\n%d ",mnode);
            for (j = 0; j < mnode; j++) {
                fprintf(fp, "%d ",c2n[etype][mnode*i + j]-1);
            }
        }
    }
    
    
    fprintf(fp, "\n");
    
    fprintf(fp, "\nCELL_TYPES %d",total_cells);

    for (etype = Triangle ; etype <= Hex; etype++) {
        for (i = 0; i < nelem[etype]; i++) {
            fprintf(fp, "\n%d",vtk_etype[etype]);
        }
    }

     
    fprintf(fp, "\n");

    fclose(fp);
    return 1;
}