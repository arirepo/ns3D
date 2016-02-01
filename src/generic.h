#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#define NUMBER_OF_ELEM_TYPES 7
typedef int ElemType;
#define Invalid_Element_Type -1
#define Edge     0
#define Triangle 1
#define Quad     2
#define Tet      3
#define Pyramid  4
#define Prism    5
#define Hex      6

#define MAX_ELEM_FACES 6
#define MAX_ELEM_EDGES 12
#define MAX_MNODE 8

#define Invalid_Element -999

#define _ASSERT(x) {assert(x);}

#define UNSIO_FLOWTYPE_UNKNOWN        -1
#define UNSIO_FLOWTYPE_INCOMPRESSIBLE 0
#define UNSIO_FLOWTYPE_COMPRESSIBLE   1
#define UNSIO_FLOWTYPE_VARMACHRUP     2

// Some function prototypes

int ElemTypeToMnode(ElemType e);
void TranslateElementWinding(int *nodes, 
                             int translation[NUMBER_OF_ELEM_TYPES][MAX_MNODE],
                             ElemType etype, int to_other_format);
int DetectNumberOfSurfaceTags(const int nelem[],
                              const int* const ifactag[],
                              int **taglist, int *maxtag_);
int IsInList(const int *list, const int nlist, const int item, int *loc);

typedef struct {
  double mach_number;
  double reynolds_number;
  double angle_of_attack;
  double elapsed_time;
} FlowParameters;

int Read_asciiugrid_Grid(FILE *fp,
			 int *nnodes,
			 int nelem[],
			 double**xyz,
			 int*c2n[],
			 int*factag[]
			 );

int Write_Fieldview_Grid(const char *filename,
			 FILE *fp,
			 const int nnodes,
			 const int nelem[],
			 const double *xyz,
			 const int * const c2n[],
			 const int * const factag[],
			 const int nID,
			 const int * IDs,
			 const char* const *names,
			 const int nvariables,
			 const int flowtype_,
			 const FlowParameters *flp,
			 const double *qvars,
			 const int nbvariables,
			 const double * const bvariables[]
			 );
///////////////////////////////////////////////////////////////////////////////
///MY CODE AND WRITING
///THE FOLLOWING DECLERATIONS ARE NOT IN THE ORIGINAL CODE BUT ADDED BY ME
///////////////////////////////////////////////////////////////////////////////
int Local_Elem_To_Global_Elem(int ielem,
                              int nelem[],
                              int etype);

int Gen_Map_Elem__Sur_Point(int *nnodes,
			 int nelem[],
			 int*c2n[],
			 int **nesp,
                         int **esp
			 );
void Global_Elem_To_Local_Elem(int *ielem,
                              int nelem[],
                              int gelem,
                              int *etype);
//Vector linked-list wrapper for holding integer arrays of variable size
//it can hold vector of same type or an array of integer. The size of enteties
//is specified by length
struct int_vec
{
    int length;
    int *int_data;
    struct int_vec* nexts;
};

//Vector linked-list wrapper for holding double arrays of variable size
//it can hold vector of same type or an array of doubles. The size of enteties
//is specified by length
struct double_vec
{
    int length;
    double *db_data;
    struct double_vec* nexts;
};

int Create_Element_Node_Conn_Table(struct int_vec **elem_node_conn);

int Gen_Map_Point__Sur_Point_LookUp(int *nnodes,
			 int nelem[],
			 int*c2n[],
			 int **nesp,
                         int **esp,
			 int your_point, int **surr_points, struct int_vec **enct);

int Gen_Map_Point__Sur_Point(int *nnodes,
			 int nelem[],
			 int*c2n[],
			 int **nesp,
                         int **esp,
			 struct int_vec **psp);

int Gen_Map_Point__Sur_Point_CRS(int *nnodes,
			 int nelem[],
			 int *c2n[],
			 int **nesp,
                         int **esp,
			 int **npsp,
                         int **psp);
int Gen_Edg_Point(struct int_vec **psp,int *nedge, int **e2n);

int Gen_Map_Elem__Sur_Elem(int nelem[],
			 int *c2n[],
			 int **nesp,
                         int **esp,
			 struct int_vec **ese);

int Create_Edges_Surr_Point(int *nnodes, int *nedge, int **e2n, struct int_vec **edsp);

int Create_Surr_Face_El_Node(struct int_vec **surr_face_el_node);

int Create_Points_In_Face_El(struct int_vec **pnts_face_el);

int obtain_common_ints(int *arr1, int *N1, int *arr2, int *N2, int *result);

int Create_Node_Centered_Grid(int *nnodes,int nelem[],
			 int **nesp,
                         int **esp,
			 int *c2n[],
                         const double *xyz,
                         int *nedge,
                         int **e2n,
                         struct double_vec **node_centered, //output
                         struct double_vec **bn_struct);     //output

int Node_Based_To_Array_Convertor(int *nnodes, int *nedge,                              // input variables
                         int **e2n, int nelem[], int *c2n[],                            // input variables
                         struct int_vec **edsp, int*factag[],                           // input variables
                         struct double_vec **node_centered,                             // input variables
                         struct double_vec **bn_struct,                                 // input variables
                         int **bn_edge_to_node_map, double **Vol,                       // output variables
                         double **Sx_i, double **Sy_i, double **Sz_i, double **S_i,     // output variables
                         double **Sx_b, double **Sy_b, double **Sz_b, double **S_b,     // output variables
                         int *bn_size, int **tag_b);                                    // output variables

#define BC_Wall 1
#define BC_FreeFace 2


//A userdefined sign function, which returns 1. for x>0. and 0. for x=.0 and -1. for x<0.
//Note that we can use copysign(x,y) in math.h but the latter is C99 and is not supported
//in all platforms
//#define sign2(x) (double)(( x > 0. ) - ( x < 0. ))
