#include "generic.h"

/////////////TEST SEGMENT OF THE ORIGINAL CODE/////////////////////////////////
///INSTRUCTIONS: MOVE THE FOLLOWING CODES TO THE CORRESPONDING HEADER AND
///IMPLEMENTATION FILES AFTER VALIDATION OF EACH FUNCTION.
///////////////////////////////////////////////////////////////////////////////


//this function initializes the conservative variables
int Init_Flow_Field(int *nnodes, int *bn_size, double *gamma, double **qinf, double *R, double *Cv, double *xyz, int **bn_edge_to_node_map, int **tag_b, double **Q1i, double **Q2i, double **Q3i, double **Q4i, double **Q5i, // input/output variables (conservatives vars for interior nodes)
                         double **Q1b, double **Q2b, double **Q3b, double **Q4b, double **Q5b) // input/output variables (conservatives vars for boundary nodes)
{
    //introducing local variables
    int i, node1;

    //declaring and initializing free stream variables
    double rho_inf = (*qinf)[1];
    double u_inf = (*qinf)[2];
    double v_inf = (*qinf)[3];
    double w_inf = (*qinf)[4];
    double p_inf = (*qinf)[5];
    double T_inf = p_inf / ((*R) * rho_inf);
    double et_inf = (*Cv)*T_inf + .5 *(pow(u_inf,2.) + pow(v_inf,2.) + pow(w_inf,2.));

    //initializing all nodes except ghost nodes
    for(i = 1; i<= (*nnodes); i++)
    {

        //Updating final conservative variables
        (*Q1i)[i] = rho_inf;
        (*Q2i)[i] = rho_inf * u_inf;
        (*Q3i)[i] = rho_inf * v_inf;
        (*Q4i)[i] = rho_inf * w_inf;
        (*Q5i)[i] = rho_inf * et_inf;
    }

    //correcting bcs and updating ghost nodes
    for(i = 1; i<= (*bn_size); i++)
    {
        //finding the boundary node(real) on that ghost edge
        node1 = (*bn_edge_to_node_map)[i];

        if((*tag_b)[i] == 3) //viscous wall
        {
        //Updating final conservative variables
        (*Q1i)[node1] = rho_inf;
        (*Q2i)[node1] = 0.;
        (*Q3i)[node1] = 0.;
        (*Q4i)[node1] = 0.;
        (*Q5i)[node1] = rho_inf * (*Cv)*T_inf;

        }

        (*Q1b)[i] = (*Q1i)[node1];
        (*Q2b)[i] = (*Q2i)[node1];
        (*Q3b)[i] = (*Q3i)[node1];
        (*Q4b)[i] = (*Q4i)[node1];
        (*Q5b)[i] = (*Q5i)[node1];

    }

//completed successfully!
return 0;
}

//this function initializes the interior conservative variables based on boundary tags, if not on boundary the node gets 0 other wise the node gets boundary element tag.
int Init_Tags(int *nnodes, int *bn_size, int **tag_b, int **bn_edge_to_node_map, double **Q1i, double **Q2i, double **Q3i, double **Q4i, double **Q5i) // input/output variables (conservatives vars for interior nodes)

{
    //introducing local variables
    int i;
    int inode = 0;

    //initializing all nodes except ghost nodes
    for(i = 1; i<= (*nnodes); i++)
    {
        (*Q1i)[i] = 0.;
        (*Q2i)[i] = 0.;
        (*Q3i)[i] = 0.;//.5 * (*Q1i)[i];
        (*Q4i)[i] = 0.;//.5 * (*Q1i)[i];
        (*Q5i)[i] = 0.;
    }

    for(i = 1; i<= (*bn_size); i++)
    {
        inode = (*bn_edge_to_node_map)[i];

        (*Q1i)[inode] = (*tag_b)[i];
        (*Q2i)[inode] = 0.;
        (*Q3i)[inode] = 0.;//.5 * (*Q1i)[i];
        (*Q4i)[inode] = 0.;//.5 * (*Q1i)[i];
        (*Q5i)[inode] = 0.;

    }

//completed successfully!
return 0;
}

//this function calculates the matrix product of two squire matrices A(N*N) and B(N*N) and puts the result 
//into C(N*N)

int Sq_Matrix_Product(double A[][6], double B[][6], double C[][6], int N)
{
    double temp;
    int i,j,k;
    
    for(i = 1; i <= N; i++)
        for(j = 1; j <= N; j++)
        {
            temp = 0.;
            for(k = 1;k <= N; k++)
                temp += (A[i][k] * B[k][j]);
            
            C[i][j] = temp;
        }
    
    //completed successfully
    return 0;
    
}

//This function calculates the residiuals
int Calc_Residuals(int *nnodes, int *nedge,                                                                                                      // input variables
                         int **e2n,                                                                                                              // input variables
                         int **bn_edge_to_node_map,int *bn_size,                                                                                 // input variables
                         double **Sx_i, double **Sy_i, double **Sz_i, double **S_i,                                                              // input variables
                         double **Sx_b, double **Sy_b, double **Sz_b, double **S_b,                                                              // input variables
                         double **Q1i, double **Q2i, double **Q3i, double **Q4i, double **Q5i,                                                   // input variables (conservatives vars for interior nodes)
                         double **Q1b, double **Q2b, double **Q3b, double **Q4b, double **Q5b,                                                   // input variables (conservatives vars for boundary nodes)
                         double *Cv, double *Cp, double *R, double *gamma,                                                                       // input vars (fluid characteristics)
                         double **Wx_e, double **Wy_e, double **Wz_e,                                      //weight for gradient prestored for each edge
                         double *xyz, double **uu1, double **uu2, double **uu3, double **uu4, double **uu5, double **grad_u1, double **grad_u2, double **grad_u3, double **grad_u4, double **grad_u5,  //preallocated gradient-related arrays
                         double **res1, double **res2, double **res3, double **res4, double **res5,double **lambda0i_S,                          // output var, residuals at all nodes, lambda * S
                         double **lambda1_ghost, double **lambda2_ghost, double **lambda3_ghost, double **lambda4_ghost, double **lambda5_ghost) // output var, lambdai at ghost edges
{
//decelaring local variables
int i, j, iedge, node1, node2; //indexers

//combinations used in eigen-value based computations
double phi =0., Thetan=0.;

//Primitive
double rho1, u1, v1, w1, et1, p1, ht1, c1, T1;
double rho2, u2, v2, w2, et2, p2, ht2, c2, T2;

//viscous terms
double mu = 1.8208e-5*40.0, kappa = 0.0257/40.0; //physical properties !!!! should be revised as a input variable to function
double txx,tyy,tzz,txy,txz,tyz,qx,qy,qz; //average viscous terms per edge
double u1_x,v1_x,w1_x,T1_x;
double u1_y,v1_y,w1_y,T1_y;
double u1_z,v1_z,w1_z,T1_z;

double u2_x,v2_x,w2_x,T2_x;
double u2_y,v2_y,w2_y,T2_y;
double u2_z,v2_z,w2_z,T2_z;

double uu_x,vv_x,ww_x,TT_x;
double uu_y,vv_y,ww_y,TT_y;
double uu_z,vv_z,ww_z,TT_z;

double vis_Minv[6][6];
double vis_f[6], vis_g[6], vis_h[6];

//hat variables (Roe vars)
double rho_hat, u_hat, v_hat, w_hat, ht_hat, c_hat;

//Allocating matrices
double M[6][6], Minv[6][6], P[6][6], Pinv[6][6], T[6][6], Tinv[6][6], A[6][6], Lambda[6][6], TempMatrix[6][6];
double f1[6], g1[6], h1[6]; //flux vectors evaluated at point 1 across the edge
double f2[6], g2[6], h2[6]; //flux vectors evaluated at point 2 across the edge
double F_ave[6], half_AQ[6]; //temp vectors to evaluate phi from formula

//// Important ////////////
//resetting residuals and lambda0i_S to zero
for(i = 1; i<= (*nnodes); i++)
{
    (*res1)[i] = 0.0;
    (*res2)[i] = 0.0;
    (*res3)[i] = 0.0;
    (*res4)[i] = 0.0;
    (*res5)[i] = 0.0;
    (*lambda0i_S)[i] = 0.0;
}
//resetting eigen-values at ghost edges
for(i = 1; i<= (*bn_size); i++)
{
    (*lambda1_ghost)[i] = 0.0;
    (*lambda2_ghost)[i] = 0.0;
    (*lambda3_ghost)[i] = 0.0;
    (*lambda4_ghost)[i] = 0.0;
    (*lambda5_ghost)[i] = 0.0;
}

///////////////////////////
//introducing average flux vector
double phi1[6];
double phi2[6];

//coordinate vars used for computing deltax,y,z for gradients.
double x1,y1,z1;
double x2,y2,z2;

//calling gradient calculater function to obtain grad Qi and store in grad_ui
    Calc_Grad_u(nnodes,                 //input - number of node
                         nedge,             //input- number of edges
                          e2n,              //input- edge to node map
                         Wx_e, Wy_e, Wz_e, //input- weight function used to compute gradients (edge-based)
                         Q1i,Q2i, Q3i, Q4i, Q5i, //input u vectors
                         grad_u1, grad_u2, grad_u3, grad_u4, grad_u5);   //output grad(u) vectors


//Precomputing vars uL,R to be used in main edge loop.
for(iedge = 1; iedge <= (*nedge); iedge++)
{

    //Resolving point1 and point2
    node1 = (*e2n)[2*iedge+ 0];
    node2 = (*e2n)[2*iedge+ 1];

    //computing coordinates point1 and point2
    x1 = xyz[3*node1 + 0]; //getting the x-coor of point1
    y1 = xyz[3*node1 + 1]; //getting the y-coor of point1
    z1 = xyz[3*node1 + 2]; //getting the z-coor of point1

    x2 = xyz[3*node2 + 0]; //getting the x-coor of point2
    y2 = xyz[3*node2 + 1]; //getting the y-coor of point2
    z2 = xyz[3*node2 + 2]; //getting the z-coor of point2


    (*uu1)[2*iedge + 0] = (*Q1i)[node1] + (*grad_u1)[3*node1 + 0] *0.5* (x2-x1) + (*grad_u1)[3*node1 + 1] *0.5* (y2-y1) + (*grad_u1)[3*node1 + 2] *0.5* (z2-z1);  //u1L
    (*uu2)[2*iedge + 0] = (*Q2i)[node1] + (*grad_u2)[3*node1 + 0] *0.5* (x2-x1) + (*grad_u2)[3*node1 + 1] *0.5* (y2-y1) + (*grad_u2)[3*node1 + 2] *0.5* (z2-z1);  //u2L
    (*uu3)[2*iedge + 0] = (*Q3i)[node1] + (*grad_u3)[3*node1 + 0] *0.5* (x2-x1) + (*grad_u3)[3*node1 + 1] *0.5* (y2-y1) + (*grad_u3)[3*node1 + 2] *0.5* (z2-z1);  //u3L
    (*uu4)[2*iedge + 0] = (*Q4i)[node1] + (*grad_u4)[3*node1 + 0] *0.5* (x2-x1) + (*grad_u4)[3*node1 + 1] *0.5* (y2-y1) + (*grad_u4)[3*node1 + 2] *0.5* (z2-z1);  //u4L
    (*uu5)[2*iedge + 0] = (*Q5i)[node1] + (*grad_u5)[3*node1 + 0] *0.5* (x2-x1) + (*grad_u5)[3*node1 + 1] *0.5* (y2-y1) + (*grad_u5)[3*node1 + 2] *0.5* (z2-z1);  //u5L

    (*uu1)[2*iedge + 1] = (*Q1i)[node2] + (*grad_u1)[3*node2 + 0] *0.5* (-(x2-x1)) + (*grad_u1)[3*node2 + 1] *0.5* (-(y2-y1)) + (*grad_u1)[3*node2 + 2] *0.5* (-(z2-z1));  //u1R
    (*uu2)[2*iedge + 1] = (*Q2i)[node2] + (*grad_u2)[3*node2 + 0] *0.5* (-(x2-x1)) + (*grad_u2)[3*node2 + 1] *0.5* (-(y2-y1)) + (*grad_u2)[3*node2 + 2] *0.5* (-(z2-z1));  //u2R
    (*uu3)[2*iedge + 1] = (*Q3i)[node2] + (*grad_u3)[3*node2 + 0] *0.5* (-(x2-x1)) + (*grad_u3)[3*node2 + 1] *0.5* (-(y2-y1)) + (*grad_u3)[3*node2 + 2] *0.5* (-(z2-z1));  //u3R
    (*uu4)[2*iedge + 1] = (*Q4i)[node2] + (*grad_u4)[3*node2 + 0] *0.5* (-(x2-x1)) + (*grad_u4)[3*node2 + 1] *0.5* (-(y2-y1)) + (*grad_u4)[3*node2 + 2] *0.5* (-(z2-z1));  //u4R
    (*uu5)[2*iedge + 1] = (*Q5i)[node2] + (*grad_u5)[3*node2 + 0] *0.5* (-(x2-x1)) + (*grad_u5)[3*node2 + 1] *0.5* (-(y2-y1)) + (*grad_u5)[3*node2 + 2] *0.5* (-(z2-z1));  //u5R

} //end of uL,R precomputation loop

//Looping over all interior edges (except ghost edges)
for(iedge = 1; iedge <= (*nedge); iedge++)
{
    //resetting all temp matrices
    for(i = 0; i<=5; i++)
        for(j = 0; j<=5; j++)
        {
            M[i][j] = 0.;
            Minv[i][j] = 0.;
            P[i][j] = 0.;
            Pinv[i][j] = 0.;
            T[i][j] = 0.;
            Tinv[i][j] = 0.;
            A[i][j] = 0.;
            Lambda[i][j] = 0.;
            TempMatrix[i][j] = 0.;
            vis_Minv[i][j] = 0.;
        }


    //Resolving point1 and point2
    node1 = (*e2n)[2*iedge+ 0];
    node2 = (*e2n)[2*iedge+ 1];

    //Computing primitive variables from conservative vars for point 1 [Q] -> [q]
    rho1 = (*uu1)[2*iedge + 0];
    u1 = (*uu2)[2*iedge + 0] / (*uu1)[2*iedge + 0];
    v1 = (*uu3)[2*iedge + 0] / (*uu1)[2*iedge + 0];
    w1 = (*uu4)[2*iedge + 0] / (*uu1)[2*iedge + 0];
    et1 = (*uu5)[2*iedge + 0] / (*uu1)[2*iedge + 0];
    T1 = 1./(*Cv)*(et1 - 1./2.*(pow(u1,2.) + pow(v1,2.) + pow(w1,2.))); //def total energy
    p1 = rho1 * (*R) * T1; //Perfect gas assumption
    ht1 = (*Cp)*T1 + 1./2. * (pow(u1,2.) + pow(v1,2.) + pow(w1,2.));
    c1 = sqrt((*gamma) * p1/rho1);

    //Computing primitive variables from conservative vars for point 2 [Q] -> [q]
    rho2 = (*uu1)[2*iedge + 1];
    u2 = (*uu2)[2*iedge + 1] / (*uu1)[2*iedge + 1];
    v2 = (*uu3)[2*iedge + 1] / (*uu1)[2*iedge + 1];
    w2 = (*uu4)[2*iedge + 1] / (*uu1)[2*iedge + 1];
    et2 = (*uu5)[2*iedge + 1] / (*uu1)[2*iedge + 1];
    T2 = 1./(*Cv)*(et2 - 1./2.*(pow(u2,2.) + pow(v2,2.) + pow(w2,2.))); //def total energy
    p2 = rho2 * (*R) * T2; //Perfect gas assumption
    ht2 = (*Cp)*T2 + 1./2. * (pow(u2,2.) + pow(v2,2.) + pow(w2,2.));
    c2 = sqrt((*gamma) * p2/rho2);

    //evaluating flux vectors at point 1
    f1[1] = rho1 * u1;
    f1[2] = rho1 * pow(u1,2.) + p1;
    f1[3] = rho1 * u1 * v1;
    f1[4] = rho1 * u1 * w1;
    f1[5] = rho1 * ht1 * u1;

    g1[1] = rho1 * v1;
    g1[2] = rho1 * u1 * v1;
    g1[3] = rho1 * pow(v1,2.) + p1;
    g1[4] = rho1 * v1 * w1;
    g1[5] = rho1 * ht1 * v1;

    h1[1] = rho1 * w1;
    h1[2] = rho1 * u1 * w1;
    h1[3] = rho1 * w1 * v1;
    h1[4] = rho1 * pow(w1,2.) + p1;
    h1[5] = rho1 * ht1 * w1;

    //evaluating flux vectors at point 2
    f2[1] = rho2 * u2;
    f2[2] = rho2 * pow(u2,2.) + p2;
    f2[3] = rho2 * u2 * v2;
    f2[4] = rho2 * u2 * w2;
    f2[5] = rho2 * ht2 * u2;

    g2[1] = rho2 * v2;
    g2[2] = rho2 * u2 * v2;
    g2[3] = rho2 * pow(v2,2.) + p2;
    g2[4] = rho2 * v2 * w2;
    g2[5] = rho2 * ht2 * v2;

    h2[1] = rho2 * w2;
    h2[2] = rho2 * u2 * w2;
    h2[3] = rho2 * w2 * v2;
    h2[4] = rho2 * pow(w2,2.) + p2;
    h2[5] = rho2 * ht2 * w2;

    
    //Computing hat variables (Roe variables or averages)
    rho_hat = sqrt(rho1*rho2);
    u_hat = (sqrt(rho1) *u1 + sqrt(rho2) *u2)  / (sqrt(rho1)+sqrt(rho2));
    v_hat = (sqrt(rho1) *v1 + sqrt(rho2) *v2)  / (sqrt(rho1)+sqrt(rho2));
    w_hat = (sqrt(rho1) *w1 + sqrt(rho2) *w2)  / (sqrt(rho1)+sqrt(rho2));
    ht_hat = (sqrt(rho1) *ht1 + sqrt(rho2) *ht2)  / (sqrt(rho1)+sqrt(rho2));
    c_hat = sqrt(((*gamma) - 1) * (ht_hat - 1./2. *(pow(u_hat,2.) + pow(v_hat,2.) + pow(w_hat,2.))));

    //Evaluating P, Pinv at Roe averages
    P[1][1] = (*Sx_i)[iedge];
    P[1][2] = (*Sy_i)[iedge];
    P[1][3] = (*Sz_i)[iedge];
    P[1][4] = rho_hat/c_hat;
    P[1][5] = rho_hat/c_hat;
    P[2][1] = 0.;
    P[2][2] = -(*Sz_i)[iedge];
    P[2][3] = (*Sy_i)[iedge];
    P[2][4] = (*Sx_i)[iedge];
    P[2][5] = -(*Sx_i)[iedge];
    P[3][1] = (*Sz_i)[iedge];
    P[3][2] = 0.;
    P[3][3] = -(*Sx_i)[iedge];
    P[3][4] = (*Sy_i)[iedge];
    P[3][5] = -(*Sy_i)[iedge];
    P[4][1] = -(*Sy_i)[iedge];
    P[4][2] = (*Sx_i)[iedge];
    P[4][3] = 0.;
    P[4][4] = (*Sz_i)[iedge];
    P[4][5] = -(*Sz_i)[iedge];
    P[5][1] = 0.;
    P[5][2] = 0.;
    P[5][3] = 0.;
    P[5][4] = rho_hat*c_hat;
    P[5][5] = rho_hat*c_hat;
    //filling inverse of P[]
    Pinv[1][1] = (*Sx_i)[iedge];
    Pinv[1][2] = 0.;
    Pinv[1][3] = (*Sz_i)[iedge];
    Pinv[1][4] = -(*Sy_i)[iedge];
    Pinv[1][5] = -(*Sx_i)[iedge]/pow(c_hat,2.);
    Pinv[2][1] = (*Sy_i)[iedge];
    Pinv[2][2] = -(*Sz_i)[iedge];
    Pinv[2][3] = 0.;
    Pinv[2][4] = (*Sx_i)[iedge];
    Pinv[2][5] = -(*Sy_i)[iedge]/pow(c_hat,2.);
    Pinv[3][1] = (*Sz_i)[iedge];
    Pinv[3][2] = (*Sy_i)[iedge];
    Pinv[3][3] = -(*Sx_i)[iedge];
    Pinv[3][4] = 0.;
    Pinv[3][5] = -(*Sz_i)[iedge]/pow(c_hat,2.);
    Pinv[4][1] = 0.;
    Pinv[4][2] = (*Sx_i)[iedge]/2.0;
    Pinv[4][3] = (*Sy_i)[iedge]/2.0;
    Pinv[4][4] = (*Sz_i)[iedge]/2.0;
    Pinv[4][5] = 1./(2. * rho_hat * c_hat);
    Pinv[5][1] = 0.;
    Pinv[5][2] = -(*Sx_i)[iedge]/2.0;
    Pinv[5][3] = -(*Sy_i)[iedge]/2.0;
    Pinv[5][4] = -(*Sz_i)[iedge]/2.0;
    Pinv[5][5] = 1./(2.*rho_hat*c_hat);

    //Evaluating M, Minv at Roe averages
    phi = ((*gamma) -1.)/2.*(pow(u_hat,2.) + pow(v_hat,2.) + pow(w_hat,2.));
    
    M[1][1] = 1.0;
    M[1][2] = 0.0;
    M[1][3] = 0.0;
    M[1][4] = 0.0;
    M[1][5] = 0.0;
    M[2][1] = u_hat;
    M[2][2] = rho_hat;
    M[2][3] = 0.;
    M[2][4] = 0.;
    M[2][5] = 0.;
    M[3][1] = v_hat;
    M[3][2] = 0.;
    M[3][3] = rho_hat;
    M[3][4] = 0.;
    M[3][5] = 0.;
    M[4][1] = w_hat;
    M[4][2] = 0.;
    M[4][3] = 0.;
    M[4][4] = rho_hat;
    M[4][5] = 0.;
    M[5][1] = phi/((*gamma)-1.);
    M[5][2] = rho_hat * u_hat;
    M[5][3] = rho_hat * v_hat;
    M[5][4] = rho_hat * w_hat;
    M[5][5] = 1./((*gamma)-1.);

    Minv[1][1] = 1.0;
    Minv[1][2] = 0.0;
    Minv[1][3] = 0.0;
    Minv[1][4] = 0.0;
    Minv[1][5] = 0.0;
    Minv[2][1] = -u_hat/rho_hat;
    Minv[2][2] = 1./rho_hat;
    Minv[2][3] = 0.;
    Minv[2][4] = 0.;
    Minv[2][5] = 0.;
    Minv[3][1] = -v_hat/rho_hat;
    Minv[3][2] = 0.;
    Minv[3][3] = 1./rho_hat;
    Minv[3][4] = 0.;
    Minv[3][5] = 0.;
    Minv[4][1] = -w_hat/rho_hat;
    Minv[4][2] = 0.;
    Minv[4][3] = 0.;
    Minv[4][4] = 1./rho_hat;
    Minv[4][5] = 0.;
    Minv[5][1] = phi;
    Minv[5][2] = -u_hat*((*gamma)-1.);
    Minv[5][3] = -v_hat*((*gamma)-1.);
    Minv[5][4] = -w_hat*((*gamma)-1.);
    Minv[5][5] = ((*gamma)-1.);

    //calculating T Tinv
    Sq_Matrix_Product(M, P, T,5);
    Sq_Matrix_Product(Pinv, Minv, Tinv,5);
    
    //calculatig Lambda matrix
    Thetan = u_hat*(*Sx_i)[iedge] + v_hat*(*Sy_i)[iedge] + w_hat*(*Sz_i)[iedge];

    Lambda[1][1] = fabs(Thetan);
    Lambda[2][2] = fabs(Thetan);
    Lambda[3][3] = fabs(Thetan);
    Lambda[4][4] = fabs(Thetan + c_hat);
    Lambda[5][5] = fabs(Thetan - c_hat);

    //calculating A, we need a temporary matrix TempMatrix
    //step 1: calculating Lambda T
    Sq_Matrix_Product(Lambda,Tinv,TempMatrix,5);
    //step 2: calculatingA =  T^-1 TempMatrix
    Sq_Matrix_Product(T,TempMatrix,A,5);

    //////////////////////////////////////////////////
    //calculating viscous terms
    //point 1
    //step 1: compute vis_Minv matrix at node 1 to find derivative of primitive using already-computed derivative of conservatives
    vis_Minv[1][1] = 1.0; vis_Minv[1][2] = 0.0; vis_Minv[1][3] = 0.0; vis_Minv[1][4] = 0.0; vis_Minv[1][5] = 0.0; //first row
    vis_Minv[2][1] = -u1/rho1; vis_Minv[2][2] = 1.0/rho1; vis_Minv[2][3] = 0.0; vis_Minv[2][4] = 0.0; vis_Minv[2][5] = 0.0; //second row
    vis_Minv[3][1] = -v1/rho1; vis_Minv[3][2] = 0.0; vis_Minv[3][3] = 1.0/rho1; vis_Minv[3][4] = 0.0; vis_Minv[3][5] = 0.0; //third row
    vis_Minv[4][1] = -w1/rho1; vis_Minv[4][2] = 0.; vis_Minv[4][3] = 0.; vis_Minv[4][4] = 1./rho1; vis_Minv[4][5] = 0.0; //fourth row
    vis_Minv[5][1] = -T1/rho1 + 1./(2.*rho1*(*Cv))*(pow(u1,2.) + pow(v1,2.) + pow(w1,2.)); vis_Minv[5][2] = -u1/(rho1*(*Cv)); vis_Minv[5][3] = -v1/(rho1*(*Cv)); vis_Minv[5][4] = -w1/(rho1*(*Cv)); vis_Minv[5][5] = 1./(rho1*(*Cv)); //fifth row
    ///step 2: compute primitive gradients for node 1
    u1_x = vis_Minv[2][1] * (*grad_u1)[3*node1 + 0] + vis_Minv[2][2] * (*grad_u2)[3*node1 + 0] + vis_Minv[2][3] * (*grad_u3)[3*node1 + 0] + vis_Minv[2][4] * (*grad_u4)[3*node1 + 0] +  vis_Minv[2][5] * (*grad_u5)[3*node1 + 0];
    u1_y = vis_Minv[2][1] * (*grad_u1)[3*node1 + 1] + vis_Minv[2][2] * (*grad_u2)[3*node1 + 1] + vis_Minv[2][3] * (*grad_u3)[3*node1 + 1] + vis_Minv[2][4] * (*grad_u4)[3*node1 + 1] +  vis_Minv[2][5] * (*grad_u5)[3*node1 + 1];
    u1_z = vis_Minv[2][1] * (*grad_u1)[3*node1 + 2] + vis_Minv[2][2] * (*grad_u2)[3*node1 + 2] + vis_Minv[2][3] * (*grad_u3)[3*node1 + 2] + vis_Minv[2][4] * (*grad_u4)[3*node1 + 2] +  vis_Minv[2][5] * (*grad_u5)[3*node1 + 2];

    v1_x = vis_Minv[3][1] * (*grad_u1)[3*node1 + 0] + vis_Minv[3][2] * (*grad_u2)[3*node1 + 0] + vis_Minv[3][3] * (*grad_u3)[3*node1 + 0] + vis_Minv[3][4] * (*grad_u4)[3*node1 + 0] +  vis_Minv[3][5] * (*grad_u5)[3*node1 + 0];
    v1_y = vis_Minv[3][1] * (*grad_u1)[3*node1 + 1] + vis_Minv[3][2] * (*grad_u2)[3*node1 + 1] + vis_Minv[3][3] * (*grad_u3)[3*node1 + 1] + vis_Minv[3][4] * (*grad_u4)[3*node1 + 1] +  vis_Minv[3][5] * (*grad_u5)[3*node1 + 1];
    v1_z = vis_Minv[3][1] * (*grad_u1)[3*node1 + 2] + vis_Minv[3][2] * (*grad_u2)[3*node1 + 2] + vis_Minv[3][3] * (*grad_u3)[3*node1 + 2] + vis_Minv[3][4] * (*grad_u4)[3*node1 + 2] +  vis_Minv[3][5] * (*grad_u5)[3*node1 + 2];

    w1_x = vis_Minv[4][1] * (*grad_u1)[3*node1 + 0] + vis_Minv[4][2] * (*grad_u2)[3*node1 + 0] + vis_Minv[4][3] * (*grad_u3)[3*node1 + 0] + vis_Minv[4][4] * (*grad_u4)[3*node1 + 0] +  vis_Minv[4][5] * (*grad_u5)[3*node1 + 0];
    w1_y = vis_Minv[4][1] * (*grad_u1)[3*node1 + 1] + vis_Minv[4][2] * (*grad_u2)[3*node1 + 1] + vis_Minv[4][3] * (*grad_u3)[3*node1 + 1] + vis_Minv[4][4] * (*grad_u4)[3*node1 + 1] +  vis_Minv[4][5] * (*grad_u5)[3*node1 + 1];
    w1_z = vis_Minv[4][1] * (*grad_u1)[3*node1 + 2] + vis_Minv[4][2] * (*grad_u2)[3*node1 + 2] + vis_Minv[4][3] * (*grad_u3)[3*node1 + 2] + vis_Minv[4][4] * (*grad_u4)[3*node1 + 2] +  vis_Minv[4][5] * (*grad_u5)[3*node1 + 2];

    T1_x = vis_Minv[5][1] * (*grad_u1)[3*node1 + 0] + vis_Minv[5][2] * (*grad_u2)[3*node1 + 0] + vis_Minv[5][3] * (*grad_u3)[3*node1 + 0] + vis_Minv[5][4] * (*grad_u4)[3*node1 + 0] +  vis_Minv[5][5] * (*grad_u5)[3*node1 + 0];
    T1_y = vis_Minv[5][1] * (*grad_u1)[3*node1 + 1] + vis_Minv[5][2] * (*grad_u2)[3*node1 + 1] + vis_Minv[5][3] * (*grad_u3)[3*node1 + 1] + vis_Minv[5][4] * (*grad_u4)[3*node1 + 1] +  vis_Minv[5][5] * (*grad_u5)[3*node1 + 1];
    T1_z = vis_Minv[5][1] * (*grad_u1)[3*node1 + 2] + vis_Minv[5][2] * (*grad_u2)[3*node1 + 2] + vis_Minv[5][3] * (*grad_u3)[3*node1 + 2] + vis_Minv[5][4] * (*grad_u4)[3*node1 + 2] +  vis_Minv[5][5] * (*grad_u5)[3*node1 + 2];

    //point 2
    //step 1: compute vis_Minv matrix at node 2 to find derivative of primitive using already-computed derivative of conservatives
    vis_Minv[1][1] = 1.0; vis_Minv[1][2] = 0.0; vis_Minv[1][3] = 0.0; vis_Minv[1][4] = 0.0; vis_Minv[1][5] = 0.0; //first row
    vis_Minv[2][1] = -u2/rho2; vis_Minv[2][2] = 1.0/rho2; vis_Minv[2][3] = 0.0; vis_Minv[2][4] = 0.0; vis_Minv[2][5] = 0.0; //second row
    vis_Minv[3][1] = -v2/rho2; vis_Minv[3][2] = 0.0; vis_Minv[3][3] = 1.0/rho2; vis_Minv[3][4] = 0.0; vis_Minv[3][5] = 0.0; //third row
    vis_Minv[4][1] = -w2/rho2; vis_Minv[4][2] = 0.; vis_Minv[4][3] = 0.; vis_Minv[4][4] = 1./rho2; vis_Minv[4][5] = 0.0; //fourth row
    vis_Minv[5][1] = -T2/rho2 + 1./(2.*rho2*(*Cv))*(pow(u2,2.) + pow(v2,2.) + pow(w2,2.)); vis_Minv[5][2] = -u2/(rho2*(*Cv)); vis_Minv[5][3] = -v2/(rho2*(*Cv)); vis_Minv[5][4] = -w2/(rho2*(*Cv)); vis_Minv[5][5] = 1./(rho2*(*Cv)); //fifth row
    ///step 2: compute primitive gradients for node 2
    u2_x = vis_Minv[2][1] * (*grad_u1)[3*node2 + 0] + vis_Minv[2][2] * (*grad_u2)[3*node2 + 0] + vis_Minv[2][3] * (*grad_u3)[3*node2 + 0] + vis_Minv[2][4] * (*grad_u4)[3*node2 + 0] +  vis_Minv[2][5] * (*grad_u5)[3*node2 + 0];
    u2_y = vis_Minv[2][1] * (*grad_u1)[3*node2 + 1] + vis_Minv[2][2] * (*grad_u2)[3*node2 + 1] + vis_Minv[2][3] * (*grad_u3)[3*node2 + 1] + vis_Minv[2][4] * (*grad_u4)[3*node2 + 1] +  vis_Minv[2][5] * (*grad_u5)[3*node2 + 1];
    u2_z = vis_Minv[2][1] * (*grad_u1)[3*node2 + 2] + vis_Minv[2][2] * (*grad_u2)[3*node2 + 2] + vis_Minv[2][3] * (*grad_u3)[3*node2 + 2] + vis_Minv[2][4] * (*grad_u4)[3*node2 + 2] +  vis_Minv[2][5] * (*grad_u5)[3*node2 + 2];

    v2_x = vis_Minv[3][1] * (*grad_u1)[3*node2 + 0] + vis_Minv[3][2] * (*grad_u2)[3*node2 + 0] + vis_Minv[3][3] * (*grad_u3)[3*node2 + 0] + vis_Minv[3][4] * (*grad_u4)[3*node2 + 0] +  vis_Minv[3][5] * (*grad_u5)[3*node2 + 0];
    v2_y = vis_Minv[3][1] * (*grad_u1)[3*node2 + 1] + vis_Minv[3][2] * (*grad_u2)[3*node2 + 1] + vis_Minv[3][3] * (*grad_u3)[3*node2 + 1] + vis_Minv[3][4] * (*grad_u4)[3*node2 + 1] +  vis_Minv[3][5] * (*grad_u5)[3*node2 + 1];
    v2_z = vis_Minv[3][1] * (*grad_u1)[3*node2 + 2] + vis_Minv[3][2] * (*grad_u2)[3*node2 + 2] + vis_Minv[3][3] * (*grad_u3)[3*node2 + 2] + vis_Minv[3][4] * (*grad_u4)[3*node2 + 2] +  vis_Minv[3][5] * (*grad_u5)[3*node2 + 2];

    w2_x = vis_Minv[4][1] * (*grad_u1)[3*node2 + 0] + vis_Minv[4][2] * (*grad_u2)[3*node2 + 0] + vis_Minv[4][3] * (*grad_u3)[3*node2 + 0] + vis_Minv[4][4] * (*grad_u4)[3*node2 + 0] +  vis_Minv[4][5] * (*grad_u5)[3*node2 + 0];
    w2_y = vis_Minv[4][1] * (*grad_u1)[3*node2 + 1] + vis_Minv[4][2] * (*grad_u2)[3*node2 + 1] + vis_Minv[4][3] * (*grad_u3)[3*node2 + 1] + vis_Minv[4][4] * (*grad_u4)[3*node2 + 1] +  vis_Minv[4][5] * (*grad_u5)[3*node2 + 1];
    w2_z = vis_Minv[4][1] * (*grad_u1)[3*node2 + 2] + vis_Minv[4][2] * (*grad_u2)[3*node2 + 2] + vis_Minv[4][3] * (*grad_u3)[3*node2 + 2] + vis_Minv[4][4] * (*grad_u4)[3*node2 + 2] +  vis_Minv[4][5] * (*grad_u5)[3*node2 + 2];

    T2_x = vis_Minv[5][1] * (*grad_u1)[3*node2 + 0] + vis_Minv[5][2] * (*grad_u2)[3*node2 + 0] + vis_Minv[5][3] * (*grad_u3)[3*node2 + 0] + vis_Minv[5][4] * (*grad_u4)[3*node2 + 0] +  vis_Minv[5][5] * (*grad_u5)[3*node2 + 0];
    T2_y = vis_Minv[5][1] * (*grad_u1)[3*node2 + 1] + vis_Minv[5][2] * (*grad_u2)[3*node2 + 1] + vis_Minv[5][3] * (*grad_u3)[3*node2 + 1] + vis_Minv[5][4] * (*grad_u4)[3*node2 + 1] +  vis_Minv[5][5] * (*grad_u5)[3*node2 + 1];
    T2_z = vis_Minv[5][1] * (*grad_u1)[3*node2 + 2] + vis_Minv[5][2] * (*grad_u2)[3*node2 + 2] + vis_Minv[5][3] * (*grad_u3)[3*node2 + 2] + vis_Minv[5][4] * (*grad_u4)[3*node2 + 2] +  vis_Minv[5][5] * (*grad_u5)[3*node2 + 2];

    //computing average primitive variable gradients at the center of the edge usind gradients at point1 and 2.
    uu_x = (u1_x + u2_x)/2.0; vv_x = (v1_x + v2_x)/2.0; ww_x = (w1_x + w2_x)/2.0; TT_x = (T1_x + T2_x)/2.0;
    uu_y = (u1_y + u2_y)/2.0; vv_y = (v1_y + v2_y)/2.0; ww_y = (w1_y + w2_y)/2.0; TT_y = (T1_y + T2_y)/2.0;
    uu_z = (u1_z + u2_z)/2.0; vv_z = (v1_z + v2_z)/2.0; ww_z = (w1_z + w2_z)/2.0; TT_z = (T1_z + T2_z)/2.0;

    //computing average viscouse fluxes at the center of the edges
    txx = 2./3. * mu * (2.* uu_x - vv_y - ww_z);
    tyy = 2./3. * mu * (2.* vv_y - uu_x - ww_z);
    tzz = 2./3. * mu * (2.* ww_z - uu_x - vv_y);
    txy = mu * (uu_y + vv_x);
    txz = mu * (ww_x + uu_z);
    tyz = mu * (vv_z + ww_y);
    qx = -kappa * TT_x;
    qy = -kappa * TT_y;
    qz = -kappa * TT_z;

    //Putting in viscous flux vector
    vis_f[1] = 0.; vis_f[2] = -txx; vis_f[3] = -txy; vis_f[4] = -txz; vis_f[5] = -0.5*(u1+u2)*txx -0.5*(v1+v2)*txy -0.5*(w1+w2)*txz + qx;
    vis_g[1] = 0.; vis_g[2] = -txy; vis_g[3] = -tyy; vis_g[4] = -tyz; vis_g[5] = -0.5*(u1+u2)*txy -0.5*(v1+v2)*tyy -0.5*(w1+w2)*tyz + qy;
    vis_h[1] = 0.; vis_h[2] = -txz; vis_h[3] = -tyz; vis_h[4] = -tzz; vis_h[5] = -0.5*(u1+u2)*txz -0.5*(v1+v2)*tyz -0.5*(w1+w2)*tzz + qz;

    //Updating total fluxes
    for (i=1; i<=5; i++)
    {
        f1[i] += vis_f[i];
        f2[i] += vis_f[i];

        g1[i] += vis_g[i];
        g2[i] += vis_g[i];

        h1[i] += vis_h[i];
        h2[i] += vis_h[i];
    }
    //end of calculating viscous fluxes
    /////////////////////////////////////

    //Computing normal flux vector averages over onde 1 and 2
    F_ave[1] = 1./2.*((f1[1]+f2[1])*(*Sx_i)[iedge] + (g1[1]+g2[1])*(*Sy_i)[iedge] + (h1[1]+h2[1])*(*Sz_i)[iedge]);
    F_ave[2] = 1./2.*((f1[2]+f2[2])*(*Sx_i)[iedge] + (g1[2]+g2[2])*(*Sy_i)[iedge] + (h1[2]+h2[2])*(*Sz_i)[iedge]);
    F_ave[3] = 1./2.*((f1[3]+f2[3])*(*Sx_i)[iedge] + (g1[3]+g2[3])*(*Sy_i)[iedge] + (h1[3]+h2[3])*(*Sz_i)[iedge]);
    F_ave[4] = 1./2.*((f1[4]+f2[4])*(*Sx_i)[iedge] + (g1[4]+g2[4])*(*Sy_i)[iedge] + (h1[4]+h2[4])*(*Sz_i)[iedge]);
    F_ave[5] = 1./2.*((f1[5]+f2[5])*(*Sx_i)[iedge] + (g1[5]+g2[5])*(*Sy_i)[iedge] + (h1[5]+h2[5])*(*Sz_i)[iedge]);

    //Computing 1/2 * A (Q2 - Q1)
    half_AQ[1] = 1./2.* A[1][1] * ((*uu1)[2*iedge + 1] - (*uu1)[2*iedge + 0]) + 1./2.* A[1][2] * ((*uu2)[2*iedge + 1] - (*uu2)[2*iedge + 0]) + 1./2.* A[1][3] * ((*uu3)[2*iedge + 1] - (*uu3)[2*iedge + 0]) + 1./2.* A[1][4] * ((*uu4)[2*iedge + 1] - (*uu4)[2*iedge + 0]) + 1./2.* A[1][5] * ((*uu5)[2*iedge + 1] - (*uu5)[2*iedge + 0]);
    half_AQ[2] = 1./2.* A[2][1] * ((*uu1)[2*iedge + 1] - (*uu1)[2*iedge + 0]) + 1./2.* A[2][2] * ((*uu2)[2*iedge + 1] - (*uu2)[2*iedge + 0]) + 1./2.* A[2][3] * ((*uu3)[2*iedge + 1] - (*uu3)[2*iedge + 0]) + 1./2.* A[2][4] * ((*uu4)[2*iedge + 1] - (*uu4)[2*iedge + 0]) + 1./2.* A[2][5] * ((*uu5)[2*iedge + 1] - (*uu5)[2*iedge + 0]);
    half_AQ[3] = 1./2.* A[3][1] * ((*uu1)[2*iedge + 1] - (*uu1)[2*iedge + 0]) + 1./2.* A[3][2] * ((*uu2)[2*iedge + 1] - (*uu2)[2*iedge + 0]) + 1./2.* A[3][3] * ((*uu3)[2*iedge + 1] - (*uu3)[2*iedge + 0]) + 1./2.* A[3][4] * ((*uu4)[2*iedge + 1] - (*uu4)[2*iedge + 0]) + 1./2.* A[3][5] * ((*uu5)[2*iedge + 1] - (*uu5)[2*iedge + 0]);
    half_AQ[4] = 1./2.* A[4][1] * ((*uu1)[2*iedge + 1] - (*uu1)[2*iedge + 0]) + 1./2.* A[4][2] * ((*uu2)[2*iedge + 1] - (*uu2)[2*iedge + 0]) + 1./2.* A[4][3] * ((*uu3)[2*iedge + 1] - (*uu3)[2*iedge + 0]) + 1./2.* A[4][4] * ((*uu4)[2*iedge + 1] - (*uu4)[2*iedge + 0]) + 1./2.* A[4][5] * ((*uu5)[2*iedge + 1] - (*uu5)[2*iedge + 0]);
    half_AQ[5] = 1./2.* A[5][1] * ((*uu1)[2*iedge + 1] - (*uu1)[2*iedge + 0]) + 1./2.* A[5][2] * ((*uu2)[2*iedge + 1] - (*uu2)[2*iedge + 0]) + 1./2.* A[5][3] * ((*uu3)[2*iedge + 1] - (*uu3)[2*iedge + 0]) + 1./2.* A[5][4] * ((*uu4)[2*iedge + 1] - (*uu4)[2*iedge + 0]) + 1./2.* A[5][5] * ((*uu5)[2*iedge + 1] - (*uu5)[2*iedge + 0]);

    //computing flux at node 1
    phi1[1] = F_ave[1] - half_AQ[1];
    phi1[2] = F_ave[2] - half_AQ[2];
    phi1[3] = F_ave[3] - half_AQ[3];
    phi1[4] = F_ave[4] - half_AQ[4];
    phi1[5] = F_ave[5] - half_AQ[5];

    //distributing flux to node 2
    phi2[1] =  phi1[1];
    phi2[2] =  phi1[2];
    phi2[3] =  phi1[3];
    phi2[4] =  phi1[4];
    phi2[5] =  phi1[5];

    //accumulating residuals to the corresponding nodes
    //for point 1 we have
    (*res1)[node1] += phi1[1] * (*S_i)[iedge];
    (*res2)[node1] += phi1[2] * (*S_i)[iedge];
    (*res3)[node1] += phi1[3] * (*S_i)[iedge];
    (*res4)[node1] += phi1[4] * (*S_i)[iedge];
    (*res5)[node1] += phi1[5] * (*S_i)[iedge];

    //for point 2 we obtain
    (*res1)[node2] += phi2[1] * (-(*S_i)[iedge]);
    (*res2)[node2] += phi2[2] * (-(*S_i)[iedge]);
    (*res3)[node2] += phi2[3] * (-(*S_i)[iedge]);
    (*res4)[node2] += phi2[4] * (-(*S_i)[iedge]);
    (*res5)[node2] += phi2[5] * (-(*S_i)[iedge]);

    //calculating lambda0i_S and contributing it to two nodes forming interior edge
    (*lambda0i_S)[node1] += ((fabs(Lambda[4][4]) >  fabs(Lambda[5][5]))?fabs(Lambda[4][4]):fabs(Lambda[5][5])) * (*S_i)[iedge];
    (*lambda0i_S)[node2] += ((fabs(Lambda[4][4]) >  fabs(Lambda[5][5]))?fabs(Lambda[4][4]):fabs(Lambda[5][5])) * (*S_i)[iedge];
    
} //end of loop over interior edges

//Looping over all ghost edges to complete flux distribution
for(iedge = 1; iedge <= (*bn_size); iedge++)
{
    //resetting all temp matrices
    for(i = 0; i<=5; i++)
        for(j = 0; j<=5; j++)
        {
            M[i][j] = 0.;
            Minv[i][j] = 0.;
            P[i][j] = 0.;
            Pinv[i][j] = 0.;
            T[i][j] = 0.;
            Tinv[i][j] = 0.;
            A[i][j] = 0.;
            Lambda[i][j] = 0.;
            TempMatrix[i][j] = 0.;
            vis_Minv[i][j] = 0.;
        }

    //Resolving point1 and point2
    node1 = (*bn_edge_to_node_map)[iedge];
    //NOTE: for node2, since [Qb] and all boundary data are based on the sequential number of ghost edges
    //that we define, we don't need to look it up.


    //Computing primitive variables from conservative vars for point 1 [Q] -> [q]
    rho1 = (*Q1i)[node1];
    u1 = (*Q2i)[node1] / (*Q1i)[node1];
    v1 = (*Q3i)[node1] / (*Q1i)[node1];
    w1 = (*Q4i)[node1] / (*Q1i)[node1];
    et1 = (*Q5i)[node1] / (*Q1i)[node1];
    T1 = 1./(*Cv)*(et1 - 1./2.*(pow(u1,2.) + pow(v1,2.) + pow(w1,2.))); //def total energy
    p1 = rho1 * (*R) * T1; //Perfect gas assumption
    ht1 = (*Cp)*T1 + 1./2. * (pow(u1,2.) + pow(v1,2.) + pow(w1,2.));
    c1 = sqrt((*gamma) * p1/rho1);

    //Computing primitive variables from conservative vars for point 2 [Q] -> [q]
    rho2 = (*Q1b)[iedge];
    u2 = (*Q2b)[iedge] / (*Q1b)[iedge];
    v2 = (*Q3b)[iedge] / (*Q1b)[iedge];
    w2 = (*Q4b)[iedge] / (*Q1b)[iedge];
    et2 = (*Q5b)[iedge] / (*Q1b)[iedge];
    T2 = 1./(*Cv)*(et2 - 1./2.*(pow(u2,2.) + pow(v2,2.) + pow(w2,2.))); //def total energy
    p2 = rho2 * (*R) * T2; //Perfect gas assumption
    ht2 = (*Cp)*T2 + 1./2. * (pow(u2,2.) + pow(v2,2.) + pow(w2,2.));
    c2 = sqrt((*gamma) * p2/rho2);

    //evaluating flux vectors at point 1
    f1[1] = rho1 * u1;
    f1[2] = rho1 * pow(u1,2.) + p1;
    f1[3] = rho1 * u1 * v1;
    f1[4] = rho1 * u1 * w1;
    f1[5] = rho1 * ht1 * u1;

    g1[1] = rho1 * v1;
    g1[2] = rho1 * u1 * v1;
    g1[3] = rho1 * pow(v1,2.) + p1;
    g1[4] = rho1 * v1 * w1;
    g1[5] = rho1 * ht1 * v1;

    h1[1] = rho1 * w1;
    h1[2] = rho1 * u1 * w1;
    h1[3] = rho1 * w1 * v1;
    h1[4] = rho1 * pow(w1,2.) + p1;
    h1[5] = rho1 * ht1 * w1;

    //evaluating flux vectors at point 2
    f2[1] = rho2 * u2;
    f2[2] = rho2 * pow(u2,2.) + p2;
    f2[3] = rho2 * u2 * v2;
    f2[4] = rho2 * u2 * w2;
    f2[5] = rho2 * ht2 * u2;

    g2[1] = rho2 * v2;
    g2[2] = rho2 * u2 * v2;
    g2[3] = rho2 * pow(v2,2.) + p2;
    g2[4] = rho2 * v2 * w2;
    g2[5] = rho2 * ht2 * v2;

    h2[1] = rho2 * w2;
    h2[2] = rho2 * u2 * w2;
    h2[3] = rho2 * w2 * v2;
    h2[4] = rho2 * pow(w2,2.) + p2;
    h2[5] = rho2 * ht2 * w2;

    //computing hat variables (Roe vars)
    rho_hat = sqrt(rho1*rho2);
    u_hat = (sqrt(rho1) *u1 + sqrt(rho2) *u2)  / (sqrt(rho1)+sqrt(rho2));
    v_hat = (sqrt(rho1) *v1 + sqrt(rho2) *v2)  / (sqrt(rho1)+sqrt(rho2));
    w_hat = (sqrt(rho1) *w1 + sqrt(rho2) *w2)  / (sqrt(rho1)+sqrt(rho2));
    ht_hat = (sqrt(rho1) *ht1 + sqrt(rho2) *ht2)  / (sqrt(rho1)+sqrt(rho2));
    c_hat = sqrt(((*gamma) - 1) * (ht_hat - 1./2. *(pow(u_hat,2.) + pow(v_hat,2.) + pow(w_hat,2.))));

    //Evaluating P, Pinv at Roe averages
    P[1][1] = (*Sx_b)[iedge];
    P[1][2] = (*Sy_b)[iedge];
    P[1][3] = (*Sz_b)[iedge];
    P[1][4] = rho_hat/c_hat;
    P[1][5] = rho_hat/c_hat;
    P[2][1] = 0.;
    P[2][2] = -(*Sz_b)[iedge];
    P[2][3] = (*Sy_b)[iedge];
    P[2][4] = (*Sx_b)[iedge];
    P[2][5] = -(*Sx_b)[iedge];
    P[3][1] = (*Sz_b)[iedge];
    P[3][2] = 0.;
    P[3][3] = -(*Sx_b)[iedge];
    P[3][4] = (*Sy_b)[iedge];
    P[3][5] = -(*Sy_b)[iedge];
    P[4][1] = -(*Sy_b)[iedge];
    P[4][2] = (*Sx_b)[iedge];
    P[4][3] = 0.;
    P[4][4] = (*Sz_b)[iedge];
    P[4][5] = -(*Sz_b)[iedge];
    P[5][1] = 0.;
    P[5][2] = 0.;
    P[5][3] = 0.;
    P[5][4] = rho_hat*c_hat;
    P[5][5] = rho_hat*c_hat;

    Pinv[1][1] = (*Sx_b)[iedge];
    Pinv[1][2] = 0.;
    Pinv[1][3] = (*Sz_b)[iedge];
    Pinv[1][4] = -(*Sy_b)[iedge];
    Pinv[1][5] = -(*Sx_b)[iedge]/pow(c_hat,2.);
    Pinv[2][1] = (*Sy_b)[iedge];
    Pinv[2][2] = -(*Sz_b)[iedge];
    Pinv[2][3] = 0.;
    Pinv[2][4] = (*Sx_b)[iedge];
    Pinv[2][5] = -(*Sy_b)[iedge]/pow(c_hat,2.);
    Pinv[3][1] = (*Sz_b)[iedge];
    Pinv[3][2] = (*Sy_b)[iedge];
    Pinv[3][3] = -(*Sx_b)[iedge];
    Pinv[3][4] = 0.;
    Pinv[3][5] = -(*Sz_b)[iedge]/pow(c_hat,2.);
    Pinv[4][1] = 0.;
    Pinv[4][2] = (*Sx_b)[iedge]/2.0;
    Pinv[4][3] = (*Sy_b)[iedge]/2.0;
    Pinv[4][4] = (*Sz_b)[iedge]/2.0;
    Pinv[4][5] = 1./(2. * rho_hat * c_hat);
    Pinv[5][1] = 0.;
    Pinv[5][2] = -(*Sx_b)[iedge]/2.0;
    Pinv[5][3] = -(*Sy_b)[iedge]/2.0;
    Pinv[5][4] = -(*Sz_b)[iedge]/2.0;
    Pinv[5][5] = 1./(2.*rho_hat*c_hat);

    //Evaluating M, Minv at Roe averages
    phi = ((*gamma) -1.)/2.*(pow(u_hat,2.) + pow(v_hat,2.) + pow(w_hat,2.));

    M[1][1] = 1.0;
    M[1][2] = 0.0;
    M[1][3] = 0.0;
    M[1][4] = 0.0;
    M[1][5] = 0.0;
    M[2][1] = u_hat;
    M[2][2] = rho_hat;
    M[2][3] = 0.;
    M[2][4] = 0.;
    M[2][5] = 0.;
    M[3][1] = v_hat;
    M[3][2] = 0.;
    M[3][3] = rho_hat;
    M[3][4] = 0.;
    M[3][5] = 0.;
    M[4][1] = w_hat;
    M[4][2] = 0.;
    M[4][3] = 0.;
    M[4][4] = rho_hat;
    M[4][5] = 0.;
    M[5][1] = phi/((*gamma)-1.);
    M[5][2] = rho_hat * u_hat;
    M[5][3] = rho_hat * v_hat;
    M[5][4] = rho_hat * w_hat;
    M[5][5] = 1./((*gamma)-1.);

    Minv[1][1] = 1.0;
    Minv[1][2] = 0.0;
    Minv[1][3] = 0.0;
    Minv[1][4] = 0.0;
    Minv[1][5] = 0.0;
    Minv[2][1] = -u_hat/rho_hat;
    Minv[2][2] = 1./rho_hat;
    Minv[2][3] = 0.;
    Minv[2][4] = 0.;
    Minv[2][5] = 0.;
    Minv[3][1] = -v_hat/rho_hat;
    Minv[3][2] = 0.;
    Minv[3][3] = 1./rho_hat;
    Minv[3][4] = 0.;
    Minv[3][5] = 0.;
    Minv[4][1] = -w_hat/rho_hat;
    Minv[4][2] = 0.;
    Minv[4][3] = 0.;
    Minv[4][4] = 1./rho_hat;
    Minv[4][5] = 0.;
    Minv[5][1] = phi;
    Minv[5][2] = -u_hat*((*gamma)-1.);
    Minv[5][3] = -v_hat*((*gamma)-1.);
    Minv[5][4] = -w_hat*((*gamma)-1.);
    Minv[5][5] = ((*gamma)-1.);

    //calculating T Tinv
    Sq_Matrix_Product(M, P, T,5);
    Sq_Matrix_Product(Pinv, Minv, Tinv,5);

    //calculatig Lambda matrix
    Thetan = u_hat*(*Sx_b)[iedge] + v_hat*(*Sy_b)[iedge] + w_hat*(*Sz_b)[iedge];

    //Storing the value of eigen-values for ghost edges based on the index of ghost edges
    (*lambda1_ghost)[iedge] = Thetan;
    (*lambda2_ghost)[iedge] = Thetan;
    (*lambda3_ghost)[iedge] = Thetan;
    (*lambda4_ghost)[iedge] = Thetan + c_hat;
    (*lambda5_ghost)[iedge] = Thetan - c_hat;


    //Calculating absolute value of eigen-values required for A^{tilde} matrix
    Lambda[1][1] = fabs(Thetan);
    Lambda[2][2] = fabs(Thetan);
    Lambda[3][3] = fabs(Thetan);
    Lambda[4][4] = fabs(Thetan + c_hat);
    Lambda[5][5] = fabs(Thetan - c_hat);

    //calculating A, we need a temporary matrix TempMatrix
    //step 1: calculating Lambda T
    Sq_Matrix_Product(Lambda,Tinv,TempMatrix,5);
    //step 2: calculatingA =  T^-1 TempMatrix
    Sq_Matrix_Product(T,TempMatrix,A,5);

        //////////////////////////////////////////////////
    //calculating viscous terms
    //point 1
    //step 1: compute vis_Minv matrix at node 1 to find derivative of primitive using already-computed derivative of conservatives
    vis_Minv[1][1] = 1.0; vis_Minv[1][2] = 0.0; vis_Minv[1][3] = 0.0; vis_Minv[1][4] = 0.0; vis_Minv[1][5] = 0.0; //first row
    vis_Minv[2][1] = -u1/rho1; vis_Minv[2][2] = 1.0/rho1; vis_Minv[2][3] = 0.0; vis_Minv[2][4] = 0.0; vis_Minv[2][5] = 0.0; //second row
    vis_Minv[3][1] = -v1/rho1; vis_Minv[3][2] = 0.0; vis_Minv[3][3] = 1.0/rho1; vis_Minv[3][4] = 0.0; vis_Minv[3][5] = 0.0; //third row
    vis_Minv[4][1] = -w1/rho1; vis_Minv[4][2] = 0.; vis_Minv[4][3] = 0.; vis_Minv[4][4] = 1./rho1; vis_Minv[4][5] = 0.0; //fourth row
    vis_Minv[5][1] = -T1/rho1 + 1./(2.*rho1*(*Cv))*(pow(u1,2.) + pow(v1,2.) + pow(w1,2.)); vis_Minv[5][2] = -u1/(rho1*(*Cv)); vis_Minv[5][3] = -v1/(rho1*(*Cv)); vis_Minv[5][4] = -w1/(rho1*(*Cv)); vis_Minv[5][5] = 1./(rho1*(*Cv)); //fifth row
    ///step 2: compute primitive gradients for node 1
    u1_x = vis_Minv[2][1] * (*grad_u1)[3*node1 + 0] + vis_Minv[2][2] * (*grad_u2)[3*node1 + 0] + vis_Minv[2][3] * (*grad_u3)[3*node1 + 0] + vis_Minv[2][4] * (*grad_u4)[3*node1 + 0] +  vis_Minv[2][5] * (*grad_u5)[3*node1 + 0];
    u1_y = vis_Minv[2][1] * (*grad_u1)[3*node1 + 1] + vis_Minv[2][2] * (*grad_u2)[3*node1 + 1] + vis_Minv[2][3] * (*grad_u3)[3*node1 + 1] + vis_Minv[2][4] * (*grad_u4)[3*node1 + 1] +  vis_Minv[2][5] * (*grad_u5)[3*node1 + 1];
    u1_z = vis_Minv[2][1] * (*grad_u1)[3*node1 + 2] + vis_Minv[2][2] * (*grad_u2)[3*node1 + 2] + vis_Minv[2][3] * (*grad_u3)[3*node1 + 2] + vis_Minv[2][4] * (*grad_u4)[3*node1 + 2] +  vis_Minv[2][5] * (*grad_u5)[3*node1 + 2];

    v1_x = vis_Minv[3][1] * (*grad_u1)[3*node1 + 0] + vis_Minv[3][2] * (*grad_u2)[3*node1 + 0] + vis_Minv[3][3] * (*grad_u3)[3*node1 + 0] + vis_Minv[3][4] * (*grad_u4)[3*node1 + 0] +  vis_Minv[3][5] * (*grad_u5)[3*node1 + 0];
    v1_y = vis_Minv[3][1] * (*grad_u1)[3*node1 + 1] + vis_Minv[3][2] * (*grad_u2)[3*node1 + 1] + vis_Minv[3][3] * (*grad_u3)[3*node1 + 1] + vis_Minv[3][4] * (*grad_u4)[3*node1 + 1] +  vis_Minv[3][5] * (*grad_u5)[3*node1 + 1];
    v1_z = vis_Minv[3][1] * (*grad_u1)[3*node1 + 2] + vis_Minv[3][2] * (*grad_u2)[3*node1 + 2] + vis_Minv[3][3] * (*grad_u3)[3*node1 + 2] + vis_Minv[3][4] * (*grad_u4)[3*node1 + 2] +  vis_Minv[3][5] * (*grad_u5)[3*node1 + 2];

    w1_x = vis_Minv[4][1] * (*grad_u1)[3*node1 + 0] + vis_Minv[4][2] * (*grad_u2)[3*node1 + 0] + vis_Minv[4][3] * (*grad_u3)[3*node1 + 0] + vis_Minv[4][4] * (*grad_u4)[3*node1 + 0] +  vis_Minv[4][5] * (*grad_u5)[3*node1 + 0];
    w1_y = vis_Minv[4][1] * (*grad_u1)[3*node1 + 1] + vis_Minv[4][2] * (*grad_u2)[3*node1 + 1] + vis_Minv[4][3] * (*grad_u3)[3*node1 + 1] + vis_Minv[4][4] * (*grad_u4)[3*node1 + 1] +  vis_Minv[4][5] * (*grad_u5)[3*node1 + 1];
    w1_z = vis_Minv[4][1] * (*grad_u1)[3*node1 + 2] + vis_Minv[4][2] * (*grad_u2)[3*node1 + 2] + vis_Minv[4][3] * (*grad_u3)[3*node1 + 2] + vis_Minv[4][4] * (*grad_u4)[3*node1 + 2] +  vis_Minv[4][5] * (*grad_u5)[3*node1 + 2];

    T1_x = vis_Minv[5][1] * (*grad_u1)[3*node1 + 0] + vis_Minv[5][2] * (*grad_u2)[3*node1 + 0] + vis_Minv[5][3] * (*grad_u3)[3*node1 + 0] + vis_Minv[5][4] * (*grad_u4)[3*node1 + 0] +  vis_Minv[5][5] * (*grad_u5)[3*node1 + 0];
    T1_y = vis_Minv[5][1] * (*grad_u1)[3*node1 + 1] + vis_Minv[5][2] * (*grad_u2)[3*node1 + 1] + vis_Minv[5][3] * (*grad_u3)[3*node1 + 1] + vis_Minv[5][4] * (*grad_u4)[3*node1 + 1] +  vis_Minv[5][5] * (*grad_u5)[3*node1 + 1];
    T1_z = vis_Minv[5][1] * (*grad_u1)[3*node1 + 2] + vis_Minv[5][2] * (*grad_u2)[3*node1 + 2] + vis_Minv[5][3] * (*grad_u3)[3*node1 + 2] + vis_Minv[5][4] * (*grad_u4)[3*node1 + 2] +  vis_Minv[5][5] * (*grad_u5)[3*node1 + 2];

    //point 2
    //step 1: compute vis_Minv matrix at node 2 to find derivative of primitive using already-computed derivative of conservatives
    ///step 2: compute primitive gradients for node 2 (since this is ghost node we don't have gradient there hence we use gradients at point 1. The accuracy is lower.
    u2_x = u1_x;
    u2_y = u1_y;
    u2_z = u1_z;

    v2_x = v1_x;
    v2_y = v1_y;
    v2_z = v1_z;

    w2_x = w1_x;
    w2_y = w1_y;
    w2_z = w1_z;

    T2_x = T1_x;
    T2_y = T1_y;
    T2_z = T1_z;

    //computing average primitive variable gradients at the center of the edge usind gradients at point1 and 2.
    uu_x = (u1_x + u2_x)/2.0; vv_x = (v1_x + v2_x)/2.0; ww_x = (w1_x + w2_x)/2.0; TT_x = (T1_x + T2_x)/2.0;
    uu_y = (u1_y + u2_y)/2.0; vv_y = (v1_y + v2_y)/2.0; ww_y = (w1_y + w2_y)/2.0; TT_y = (T1_y + T2_y)/2.0;
    uu_z = (u1_z + u2_z)/2.0; vv_z = (v1_z + v2_z)/2.0; ww_z = (w1_z + w2_z)/2.0; TT_z = (T1_z + T2_z)/2.0;

    //computing average viscouse fluxes at the center of the edges
    txx = 2./3. * mu * (2.* uu_x - vv_y - ww_z);
    tyy = 2./3. * mu * (2.* vv_y - uu_x - ww_z);
    tzz = 2./3. * mu * (2.* ww_z - uu_x - vv_y);
    txy = mu * (uu_y + vv_x);
    txz = mu * (ww_x + uu_z);
    tyz = mu * (vv_z + ww_y);
    qx = -kappa * TT_x;
    qy = -kappa * TT_y;
    qz = -kappa * TT_z;

    //Putting in viscous flux vector
    vis_f[1] = 0.; vis_f[2] = -txx; vis_f[3] = -txy; vis_f[4] = -txz; vis_f[5] = -0.5*(u1+u2)*txx -0.5*(v1+v2)*txy -0.5*(w1+w2)*txz + qx;
    vis_g[1] = 0.; vis_g[2] = -txy; vis_g[3] = -tyy; vis_g[4] = -tyz; vis_g[5] = -0.5*(u1+u2)*txy -0.5*(v1+v2)*tyy -0.5*(w1+w2)*tyz + qy;
    vis_h[1] = 0.; vis_h[2] = -txz; vis_h[3] = -tyz; vis_h[4] = -tzz; vis_h[5] = -0.5*(u1+u2)*txz -0.5*(v1+v2)*tyz -0.5*(w1+w2)*tzz + qz;

    //Updating total fluxes
    for (i=1; i<=5; i++)
    {
        f1[i] += vis_f[i];
        f2[i] += vis_f[i];

        g1[i] += vis_g[i];
        g2[i] += vis_g[i];

        h1[i] += vis_h[i];
        h2[i] += vis_h[i];
    }
    //end of calculating viscous fluxes
    /////////////////////////////////////

    //Computing normal flux vector averages over onde 1 and 2
    F_ave[1] = 1./2.*((f1[1]+f2[1])*(*Sx_b)[iedge] + (g1[1]+g2[1])*(*Sy_b)[iedge] + (h1[1]+h2[1])*(*Sz_b)[iedge]);
    F_ave[2] = 1./2.*((f1[2]+f2[2])*(*Sx_b)[iedge] + (g1[2]+g2[2])*(*Sy_b)[iedge] + (h1[2]+h2[2])*(*Sz_b)[iedge]);
    F_ave[3] = 1./2.*((f1[3]+f2[3])*(*Sx_b)[iedge] + (g1[3]+g2[3])*(*Sy_b)[iedge] + (h1[3]+h2[3])*(*Sz_b)[iedge]);
    F_ave[4] = 1./2.*((f1[4]+f2[4])*(*Sx_b)[iedge] + (g1[4]+g2[4])*(*Sy_b)[iedge] + (h1[4]+h2[4])*(*Sz_b)[iedge]);
    F_ave[5] = 1./2.*((f1[5]+f2[5])*(*Sx_b)[iedge] + (g1[5]+g2[5])*(*Sy_b)[iedge] + (h1[5]+h2[5])*(*Sz_b)[iedge]);

    //Computing 1/2 * A (Q2 - Q1)
    half_AQ[1] = 1./2.* A[1][1] * ((*Q1b)[iedge] - (*Q1i)[node1]) + 1./2.* A[1][2] * ((*Q2b)[iedge] - (*Q2i)[node1]) + 1./2.* A[1][3] * ((*Q3b)[iedge] - (*Q3i)[node1]) + 1./2.* A[1][4] * ((*Q4b)[iedge] - (*Q4i)[node1]) + 1./2.* A[1][5] * ((*Q5b)[iedge] - (*Q5i)[node1]);
    half_AQ[2] = 1./2.* A[2][1] * ((*Q1b)[iedge] - (*Q1i)[node1]) + 1./2.* A[2][2] * ((*Q2b)[iedge] - (*Q2i)[node1]) + 1./2.* A[2][3] * ((*Q3b)[iedge] - (*Q3i)[node1]) + 1./2.* A[2][4] * ((*Q4b)[iedge] - (*Q4i)[node1]) + 1./2.* A[2][5] * ((*Q5b)[iedge] - (*Q5i)[node1]);
    half_AQ[3] = 1./2.* A[3][1] * ((*Q1b)[iedge] - (*Q1i)[node1]) + 1./2.* A[3][2] * ((*Q2b)[iedge] - (*Q2i)[node1]) + 1./2.* A[3][3] * ((*Q3b)[iedge] - (*Q3i)[node1]) + 1./2.* A[3][4] * ((*Q4b)[iedge] - (*Q4i)[node1]) + 1./2.* A[3][5] * ((*Q5b)[iedge] - (*Q5i)[node1]);
    half_AQ[4] = 1./2.* A[4][1] * ((*Q1b)[iedge] - (*Q1i)[node1]) + 1./2.* A[4][2] * ((*Q2b)[iedge] - (*Q2i)[node1]) + 1./2.* A[4][3] * ((*Q3b)[iedge] - (*Q3i)[node1]) + 1./2.* A[4][4] * ((*Q4b)[iedge] - (*Q4i)[node1]) + 1./2.* A[4][5] * ((*Q5b)[iedge] - (*Q5i)[node1]);
    half_AQ[5] = 1./2.* A[5][1] * ((*Q1b)[iedge] - (*Q1i)[node1]) + 1./2.* A[5][2] * ((*Q2b)[iedge] - (*Q2i)[node1]) + 1./2.* A[5][3] * ((*Q3b)[iedge] - (*Q3i)[node1]) + 1./2.* A[5][4] * ((*Q4b)[iedge] - (*Q4i)[node1]) + 1./2.* A[5][5] * ((*Q5b)[iedge] - (*Q5i)[node1]);

    //computing flux at node 1
    phi1[1] = F_ave[1] - half_AQ[1];
    phi1[2] = F_ave[2] - half_AQ[2];
    phi1[3] = F_ave[3] - half_AQ[3];
    phi1[4] = F_ave[4] - half_AQ[4];
    phi1[5] = F_ave[5] - half_AQ[5];

    //distributing flux to node 2
    //Note : we don't need a flux distribution at ghost point because we don't update this point and we just
    //fix it via boundary conditions

    //accumulating residuals to the corresponding nodes
    //for point 1 we have
    (*res1)[node1] += phi1[1] * (*S_b)[iedge];
    (*res2)[node1] += phi1[2] * (*S_b)[iedge];
    (*res3)[node1] += phi1[3] * (*S_b)[iedge];
    (*res4)[node1] += phi1[4] * (*S_b)[iedge];
    (*res5)[node1] += phi1[5] * (*S_b)[iedge];

    //calculating lambda0i_S and contributing it to the first node forming fictious edge
    (*lambda0i_S)[node1] += ((fabs(Lambda[4][4]) >  fabs(Lambda[5][5]))?fabs(Lambda[4][4]):fabs(Lambda[5][5])) * (*S_b)[iedge];



} //end of loop for ghost edges

//competed successfully!
return 0;
}

//This function update conservative vars at ghost edge boundary points according to characteristic boundary conditions
int Impose_Boundary_Conditions(int **bn_edge_to_node_map,int *bn_size, int **tag_b,                                           // input vars - boundary spec
                               double **Sx_b, double **Sy_b, double **Sz_b,                                                   // input boundary normals
      double **lambda1_ghost, double **lambda2_ghost, double **lambda3_ghost, double **lambda4_ghost, double **lambda5_ghost, // input var, lambdai at ghost edges
                               double **Q1i, double **Q2i, double **Q3i, double **Q4i, double **Q5i,                          // input the value of conservative variables at interior nodes having ghost edges
                               double **qinf, double *Cv, double *Cp, double *R, double *gamma,                               // input primitive variables at infinity [rho_inf, u_inf, v_inf, w_inf, p_inf] and flow properties
                               double **Q1b, double **Q2b, double **Q3b, double **Q4b, double **Q5b)                          // input/output the value of conservative variables at ghost points
{
    //decelaring local variables
    int iedge = 0; //edge index, can be used as ghost edge
    int node1 = 0; //the boundary node (real) on each ghost edge.
    //Primitive
    double rho1, u1, v1, w1, et1, p1, ht1, c1, T1; //vars for real boundary point on the ghost edge
    double rho2, u2, v2, w2, et2, p2, ht2, c2, T2; //vars for ghost point on the ghost edge
    //refrence values, for each ghost edge using relation var0 = ((var(ghost) = var2) + (var(boundary node) = var1))/2.;
    double rho0 = 0., u0 = 0., v0 = 0., w0 = 0., c0 = 0.;
    // infinity conditions readable format
    double rho_inf = (*qinf)[1];
    double u_inf = (*qinf)[2];
    double v_inf = (*qinf)[3];
    double w_inf = (*qinf)[4];
    double p_inf = (*qinf)[5];

    //looping over all ghost edges to calculate update for Qib
    for(iedge = 1; iedge <= (*bn_size); iedge++) 
    {
        //finding the boundary node(real) on that ghost edge
        node1 = (*bn_edge_to_node_map)[iedge];

        //Computing primitive variables from conservative vars for point 1 [Q] -> [q]
        rho1 = (*Q1i)[node1];
        u1 = (*Q2i)[node1] / (*Q1i)[node1];
        v1 = (*Q3i)[node1] / (*Q1i)[node1];
        w1 = (*Q4i)[node1] / (*Q1i)[node1];
        et1 = (*Q5i)[node1] / (*Q1i)[node1];
        T1 = 1./(*Cv)*(et1 - 1./2.*(pow(u1,2.) + pow(v1,2.) + pow(w1,2.))); //def total energy
        p1 = rho1 * (*R) * T1; //Perfect gas assumption
        ht1 = (*Cp)*T1 + 1./2. * (pow(u1,2.) + pow(v1,2.) + pow(w1,2.));
        c1 = sqrt((*gamma) * p1/rho1);

        //Computing primitive variables from conservative vars for point 2 [Q] -> [q]
        rho2 = (*Q1b)[iedge];
        u2 = (*Q2b)[iedge] / (*Q1b)[iedge];
        v2 = (*Q3b)[iedge] / (*Q1b)[iedge];
        w2 = (*Q4b)[iedge] / (*Q1b)[iedge];
        et2 = (*Q5b)[iedge] / (*Q1b)[iedge];
        T2 = 1./(*Cv)*(et2 - 1./2.*(pow(u2,2.) + pow(v2,2.) + pow(w2,2.))); //def total energy
        p2 = rho2 * (*R) * T2; //Perfect gas assumption
        ht2 = (*Cp)*T2 + 1./2. * (pow(u2,2.) + pow(v2,2.) + pow(w2,2.));
        c2 = sqrt((*gamma) * p2/rho2);

        //calculating the refrences variables (arithmatic averages)
        rho0 = .5 * (rho1 + rho2);
        u0   = .5 * (u1 + u2);
        v0   = .5 * (v1 + v2);
        w0   = .5 * (w1 + w2);
        c0   = .5 * (c1 + c2);

        switch((*tag_b)[iedge])
        {
            case 1: //inviscid walls
        //Checking to see if the eigen values are correct for a wall condition
/*
                if( ((*lambda1_ghost)[iedge] == 0.) && ((*lambda2_ghost)[iedge] == 0.)
                        && ((*lambda3_ghost)[iedge] == 0.) && ((*lambda4_ghost)[iedge] > 0.) &&
                        ((*lambda5_ghost)[iedge] < 0.) ) //Impermeable Wall
                {
*/
                    p2 = p1+rho0*c0*((*Sx_b)[iedge] *u1 + (*Sy_b)[iedge] *v1 + (*Sz_b)[iedge] *w1);
                    rho2 = rho1 + (p2 - p1) / pow(c0,2.0);
                    u2 = u1 - (*Sx_b)[iedge] * (p2-p1) / (rho0*c0);
                    v2 = v1 - (*Sy_b)[iedge] * (p2-p1) / (rho0*c0);
                    w2 = w1 - (*Sz_b)[iedge] * (p2-p1) / (rho0*c0);

            break;
                    
            case 2: //inlet- outlet

                //Check to see what kind of free face we have. Subsonic in/out let, or supersonic in/out let.
               if( ((*lambda1_ghost)[iedge] > 0.) && ((*lambda2_ghost)[iedge] > 0.)
                        && ((*lambda3_ghost)[iedge] > 0.) && ((*lambda4_ghost)[iedge] > 0.)
                        && ((*lambda5_ghost)[iedge] < 0.)   ) //Subsonic Outflow
                {

                    p2 = p_inf;
                    rho2 = rho1 + (p2 - p1) / pow(c0,2.0);
                    u2 = u1 - (*Sx_b)[iedge] * (p2-p1) / (rho0*c0);
                    v2 = v1 - (*Sy_b)[iedge] * (p2-p1) / (rho0*c0);
                    w2 = w1 - (*Sz_b)[iedge] * (p2-p1) / (rho0*c0);


                }
                else if ( ((*lambda1_ghost)[iedge] < 0.) && ((*lambda2_ghost)[iedge] < 0.)
                        && ((*lambda3_ghost)[iedge] < 0.) && ((*lambda4_ghost)[iedge] > 0.)
                        && ((*lambda5_ghost)[iedge] < 0.)   ) //Subsonic Inflow
                {
                    p2 = .5* (p_inf+p1-rho0*c0*((*Sx_b)[iedge] *(u_inf-u1) + (*Sy_b)[iedge] *(v_inf-v1) + (*Sz_b)[iedge] *(w_inf-w1)) );
                    rho2 = rho_inf + (p2 - p_inf) / pow(c0,2.0);
                    u2 = u_inf + (*Sx_b)[iedge] * (p2-p_inf) / (rho0*c0);
                    v2 = v_inf + (*Sy_b)[iedge] * (p2-p_inf) / (rho0*c0);
                    w2 = w_inf + (*Sz_b)[iedge] * (p2-p_inf) / (rho0*c0);

                }
                else if( ((*lambda1_ghost)[iedge] < 0.) && ((*lambda2_ghost)[iedge] < 0.)
                        && ((*lambda3_ghost)[iedge] < 0.) && ((*lambda4_ghost)[iedge] < 0.) &&
                        ((*lambda5_ghost)[iedge] < 0.) ) //Supersonic Inflow
                {
                    p2 = p_inf;
                    rho2 = rho_inf;
                    u2 = u_inf;
                    v2 = v_inf;
                    w2 = w_inf;
                }
                else if( ((*lambda1_ghost)[iedge] > 0.) && ((*lambda2_ghost)[iedge] > 0.)
                        && ((*lambda3_ghost)[iedge] > 0.) && ((*lambda4_ghost)[iedge] > 0.) &&
                        ((*lambda5_ghost)[iedge] > 0.) ) //Supersonic Outflow
                {
                    p2 = p1;
                    rho2 = rho1;
                    u2 = u1;
                    v2 = v1;
                    w2 = w1;

                }
                else
                {
                    printf("\n\rUnrecognized free face boundary conditions. Strange eigen-value signs.\n\r Could not impose BCs");
                    exit(0);

                }

               break;
        }

        //Transfering primitive variables at point 2 (ghost point) to conservative ones and updating [Qbi]
        (*Q1b)[iedge] = rho2;
        (*Q2b)[iedge] = rho2*u2;
        (*Q3b)[iedge] = rho2*v2;
        (*Q4b)[iedge] = rho2*w2;
        //some calcs before obtaining final conservative variables
        T2 = p2 / ((*R) * rho2);
        et2 = (*Cv)*T2 + .5 *(pow(u2,2.) + pow(v2,2.) + pow(w2,2.));
        (*Q5b)[iedge] = rho2*et2;

    }//end of iedge loop . Yep all ghost edges are treated and the value of
    //conservative variables at ghost points Qb are updated

    double T_inf = p_inf/((*R)*rho_inf);

    //looping over all ghost edges to calculate update for Qib for viscous solid walls
    for(iedge = 1; iedge <= (*bn_size); iedge++)
    {
        //finding the boundary node(real) on that ghost edge
        node1 = (*bn_edge_to_node_map)[iedge];

        //Computing primitive variables from conservative vars for point 1 [Q] -> [q]
        rho1 = (*Q1i)[node1];
        u1 = (*Q2i)[node1] / (*Q1i)[node1];
        v1 = (*Q3i)[node1] / (*Q1i)[node1];
        w1 = (*Q4i)[node1] / (*Q1i)[node1];
        et1 = (*Q5i)[node1] / (*Q1i)[node1];
        T1 = 1./(*Cv)*(et1 - 1./2.*(pow(u1,2.) + pow(v1,2.) + pow(w1,2.))); //def total energy
        p1 = rho1 * (*R) * T1; //Perfect gas assumption
        ht1 = (*Cp)*T1 + 1./2. * (pow(u1,2.) + pow(v1,2.) + pow(w1,2.));
        c1 = sqrt((*gamma) * p1/rho1);

        //Computing primitive variables from conservative vars for point 2 [Q] -> [q]
        rho2 = (*Q1b)[iedge];
        u2 = (*Q2b)[iedge] / (*Q1b)[iedge];
        v2 = (*Q3b)[iedge] / (*Q1b)[iedge];
        w2 = (*Q4b)[iedge] / (*Q1b)[iedge];
        et2 = (*Q5b)[iedge] / (*Q1b)[iedge];
        T2 = 1./(*Cv)*(et2 - 1./2.*(pow(u2,2.) + pow(v2,2.) + pow(w2,2.))); //def total energy
        p2 = rho2 * (*R) * T2; //Perfect gas assumption
        ht2 = (*Cp)*T2 + 1./2. * (pow(u2,2.) + pow(v2,2.) + pow(w2,2.));
        c2 = sqrt((*gamma) * p2/rho2);

        //calculating the refrences variables (arithmatic averages)
        rho0 = .5 * (rho1 + rho2);
        u0   = .5 * (u1 + u2);
        v0   = .5 * (v1 + v2);
        w0   = .5 * (w1 + w2);
        c0   = .5 * (c1 + c2);

        //Checking to see if the eigen values are correct for a wall condition
        if( (*tag_b)[iedge] == 3)
        {
                    p2 = p1;//p1+rho0*c0*((*Sx_b)[iedge] *u1 + (*Sy_b)[iedge] *v1 + (*Sz_b)[iedge] *w1);
                    rho2 = rho1 + (p2 - p1) / pow(c0,2.0);
                    u2 = 0;
                    v2 = 0;
                    w2 = 0;

            //Transfering primitive variables at point 2 (ghost point) to conservative ones and updating [Qbi]
            (*Q1b)[iedge] = rho2;
            (*Q2b)[iedge] = rho2*u2;
            (*Q3b)[iedge] = rho2*v2;
            (*Q4b)[iedge] = rho2*w2;
            //some calcs before obtaining final conservative variables
            //T2 = p2 / ((*R) * rho2);
            et2 = (*Cv)*T_inf + .5 *(pow(u2,2.) + pow(v2,2.) + pow(w2,2.));
            (*Q5b)[iedge] = rho2*et2;
        }//end of iedge loop . Yep all ghost edges are treated and the value of
    }

        //conservative variables at ghost points Qb are updated

    //completed successfully!
    return 0;
}

double find_min(int *nnodes, double **lambda0i_S, double *CFL, double **Vol)
{
    int i =0;
    double min = (*CFL)*(*Vol)[1]/(*lambda0i_S)[1];

    for(i = 2; i<= (*nnodes); i++)
    {
        if(min > (*CFL)*(*Vol)[i]/(*lambda0i_S)[i])
            min = (*CFL)*(*Vol)[i]/(*lambda0i_S)[i];

    }
    return min;
}

// this function marchs the equations to steady-state solution using Jameson- Turkel- Backer RK scheme.
int Solve_RK5_Jameson(int *nnodes, int *nedge, int **e2n, int*c2n[], int nelem[],                                       // input vars - edge and node spec
                         int **bn_edge_to_node_map,int *bn_size,                                        // input vars - boundary spec
                         double **Sx_i, double **Sy_i, double **Sz_i, double **S_i, double **Vol,  // input vars- volume
                         double **Sx_b, double **Sy_b, double **Sz_b, double **S_b,                // and boundary and interior normals
                         double *Cv, double *Cp, double *R, double *gamma,                         //input vars (fluid characteristics)
                         double *CFL, int *N_ITR, short *init_flag, short *messages_on, short *time_accur_on, int *ITR_Anim,            //input vars - solver parameters
                         int **tag_b, double **qinf,                                               //input vars - boundary element tags and free stream conds
                         double *xyz,
                         double **Q1i, double **Q2i, double **Q3i, double **Q4i, double **Q5i,     // output variables (conservatives vars for interior nodes)
                         double **Q1b, double **Q2b, double **Q3b, double **Q4b, double **Q5b)     // output variables (conservatives vars for boundary nodes)           //input vars ()

{


    //Infinity variables
    double rho_inf = (*qinf)[1];
    double p_inf = (*qinf)[5];
    double T_inf = p_inf/((*R) * rho_inf);

    //declaring a pointer to the output file for residuals
    FILE *Res_out = NULL;
    // openning the output file for writing residuals
    Res_out = fopen("tmp_res.dat", "w");
    //output contour file name which change during the main loop to make a simulation
    char out_contours[20];
    int export_intervals = 1; //export contours counter starting from 1 (the first contour)
    double time_scale = 0.1*0.013/(*qinf)[2];

    //introducing local variables
    int i, inode, iedge, node1;
    double *res1, *res2, *res3, *res4, *res5;

    //RK-specific data structure
    double *Delta_t, deltat, time = 0.;
    double alpha_RK5[6];
    alpha_RK5[1] = 1.0/4.0;
    alpha_RK5[2] = 1.0/6.0;
    alpha_RK5[3] = 3.0/8.0;
    alpha_RK5[4] = 1.0/2.0;
    alpha_RK5[5] = 1.0;
    int RK_index = 0;
    double *Qn_1,*Qn_2,*Qn_3,*Qn_4,*Qn_5;
    //declaring gradient related arrays u and grad_u
    //These variables are the input of the Roe's averages
    double *Wx_e, *Wy_e, *Wz_e; //preallocated
    double *u1, *u2, *u3, *u4, *u5; //each 1 .. 2*(nedge+1) for uL and uR per edge
    double *grad_u1, *grad_u2, *grad_u3, *grad_u4, *grad_u5; //each 1 .. 3*(nnodes+1)


    double *lambda0i_S;
    double *lambda1_ghost, *lambda2_ghost, *lambda3_ghost, *lambda4_ghost, *lambda5_ghost;
    double L2norm[6] = {0., 0., 0., 0., 0., 0.};

    //allocating conservative vars (outputs) for the interior and boundary
    (*Q1i) = (double *)calloc((*nnodes)+1, sizeof(double));
    (*Q2i) = (double *)calloc((*nnodes)+1, sizeof(double));
    (*Q3i) = (double *)calloc((*nnodes)+1, sizeof(double));
    (*Q4i) = (double *)calloc((*nnodes)+1, sizeof(double));
    (*Q5i) = (double *)calloc((*nnodes)+1, sizeof(double));

    (*Q1b) = (double *)calloc((*bn_size)+1, sizeof(double));
    (*Q2b) = (double *)calloc((*bn_size)+1, sizeof(double));
    (*Q3b) = (double *)calloc((*bn_size)+1, sizeof(double));
    (*Q4b) = (double *)calloc((*bn_size)+1, sizeof(double));
    (*Q5b) = (double *)calloc((*bn_size)+1, sizeof(double));


    //allocating residuals
    res1 = (double *)calloc((*nnodes) + 1, sizeof( double));
    res2 = (double *)calloc((*nnodes) + 1, sizeof( double));
    res3 = (double *)calloc((*nnodes) + 1, sizeof( double));
    res4 = (double *)calloc((*nnodes) + 1, sizeof( double));
    res5 = (double *)calloc((*nnodes) + 1, sizeof( double));

    //Delta_t
    Delta_t = (double *)calloc((*nnodes) + 1, sizeof( double));

    //Qn
    Qn_1 = (double *)calloc((*nnodes) + 1, sizeof( double));
    Qn_2 = (double *)calloc((*nnodes) + 1, sizeof( double));
    Qn_3 = (double *)calloc((*nnodes) + 1, sizeof( double));
    Qn_4 = (double *)calloc((*nnodes) + 1, sizeof( double));
    Qn_5 = (double *)calloc((*nnodes) + 1, sizeof( double));

    //allocating gradient-related arrays q and grad_q
    u1 = (double *)calloc(2*(*nedge+1), sizeof(double));
    u2 = (double *)calloc(2*(*nedge+1), sizeof(double));
    u3 = (double *)calloc(2*(*nedge+1), sizeof(double));
    u4 = (double *)calloc(2*(*nedge+1), sizeof(double));
    u5 = (double *)calloc(2*(*nedge+1), sizeof(double));
    grad_u1 = (double *)calloc(3*(*nnodes+1), sizeof(double));
    grad_u2 = (double *)calloc(3*(*nnodes+1), sizeof(double));
    grad_u3 = (double *)calloc(3*(*nnodes+1), sizeof(double));
    grad_u4 = (double *)calloc(3*(*nnodes+1), sizeof(double));
    grad_u5 = (double *)calloc(3*(*nnodes+1), sizeof(double));


    //allocating an array for Sum(max(lambda_i)Si)
    lambda0i_S = (double *)calloc(*nnodes + 1, sizeof( double));

    //allocating arrays for eigen-values at ghost edges
    lambda1_ghost = (double *)calloc((*bn_size)+1, sizeof( double));
    lambda2_ghost = (double *)calloc((*bn_size)+1, sizeof( double));
    lambda3_ghost = (double *)calloc((*bn_size)+1, sizeof( double));
    lambda4_ghost = (double *)calloc((*bn_size)+1, sizeof( double));
    lambda5_ghost = (double *)calloc((*bn_size)+1, sizeof( double));


    //Before anything, the weights for gradient calculations are evaluated and stored (edge-based) for the rest of program.
    Calculate_W_e(nnodes,                         //input- number of nodes
                        xyz,                      //input- coordinates
                         nedge,                   //input- number of edges
                         e2n,                     //input- edge to node map
                         &Wx_e, &Wy_e, &Wz_e);    //output- weight function used to compute gradients (edge-based)

    //Initializing the flow field
    if(*init_flag)
        Init_Flow_Field(nnodes,bn_size, gamma, qinf, R, Cv, xyz, bn_edge_to_node_map, tag_b, Q1i, Q2i, Q3i,Q4i, Q5i,
                         Q1b, Q2b, Q3b, Q4b, Q5b);

   //solving the equations
   for(i = 1;i <= *N_ITR; i++) //do the same thing for N_ITR number of iterations
   {
       //resetting norms
       L2norm[1] = 0.;
       L2norm[2] = 0.;
       L2norm[3] = 0.;
       L2norm[4] = 0.;
       L2norm[5] = 0.;

       //PreRK operation: calculate the residuals at k1_i((Qn , tn)) NOTE : THIS STEP IS A RESIDUAL COMPUTATION AND ALSO IS A DELTA_T COMPUTATION FOR ALL FOUR SUBITERATION COMMING NEXT.
       Calc_Residuals(nnodes, nedge,                                                                                                      // input variables
                         e2n,                                                                                                              // input variables
                         bn_edge_to_node_map,bn_size,                                                                                 // input variables
                         Sx_i, Sy_i, Sz_i, S_i,                                                              // input variables
                         Sx_b, Sy_b, Sz_b, S_b,                                                              // input variables
                         Q1i, Q2i, Q3i, Q4i, Q5i,                                                   // input variables (conservatives vars for interior nodes)
                         Q1b, Q2b, Q3b, Q4b, Q5b,                                                   // input variables (conservatives vars for boundary nodes)
                         Cv, Cp, R, gamma,                                                                       // input vars (fluid characteristics)
                         &Wx_e, &Wy_e, &Wz_e,                                    //weight for gradient prestored for each edge
                         xyz, &u1, &u2, &u3, &u4, &u5, &grad_u1, &grad_u2, &grad_u3, &grad_u4, &grad_u5,  //preallocated gradient-related arrays
                         &res1, &res2, &res3, &res4, &res5,&lambda0i_S,                          // output var, residuals at all nodes, lambda * S
                         &lambda1_ghost, &lambda2_ghost, &lambda3_ghost, &lambda4_ghost, &lambda5_ghost); // output var, lambdai at ghost edges


       //impose viscous wall BCs
       for(iedge = 1; iedge <= (*bn_size); iedge++)
           if((*tag_b)[iedge] == 3) //viscous wall
           {
               //find the interior node due to that ghost edge
               node1 = (*bn_edge_to_node_map)[iedge];
               //setting viscous residual for that node
               res2[node1] = 0.; //u=0 du/dt = 0
               res3[node1] = 0.; //v=0 dv/dt = 0
               res4[node1] = 0.; //w=0 dw/dt = 0
               res5[node1] = res1[node1] * (*Cv) * T_inf;
           }

        
       //time-step calculations
       if(*time_accur_on) //time-accurate
       {
        //Calculating time step - for time accurate solutions
        deltat = find_min(nnodes, &lambda0i_S, CFL, Vol);
       //PreRK operation: Substep 1: Computing Delta_t based on the first iteration for entire RK scheme
        for(inode = 1; inode <= (*nnodes); inode++)
            Delta_t[inode] = deltat;
       }
       else //steady state
       {
       //PreRK operation: Substep 1: Computing Delta_t based on the first iteration for entire RK scheme
        for(inode = 1; inode <= (*nnodes); inode++)
                  Delta_t[inode] = (*CFL) * (*Vol)[inode]/lambda0i_S[inode];
       }

       //PreRK operation: Substep 2: Taking the backup of conservative variables at interior nodes
        for(inode = 1; inode <= (*nnodes); inode++)
        {
            Qn_1[inode] = (*Q1i)[inode];
            Qn_2[inode] = (*Q2i)[inode];
            Qn_3[inode] = (*Q3i)[inode];
            Qn_4[inode] = (*Q4i)[inode];
            Qn_5[inode] = (*Q5i)[inode];
        }

       for(RK_index = 1; RK_index <= 5; RK_index++ ) //RK5 loop
       {
           //computing update for subiterations
           for(inode = 1; inode <= (*nnodes); inode++)
           {
               (*Q1i)[inode] = Qn_1[inode] -alpha_RK5[RK_index] * Delta_t[inode] /(*Vol)[inode]* res1[inode];
               (*Q2i)[inode] = Qn_2[inode] -alpha_RK5[RK_index] * Delta_t[inode] /(*Vol)[inode]* res2[inode];
               (*Q3i)[inode] = Qn_3[inode] -alpha_RK5[RK_index] * Delta_t[inode] /(*Vol)[inode]* res3[inode];
               (*Q4i)[inode] = Qn_4[inode] -alpha_RK5[RK_index] * Delta_t[inode] /(*Vol)[inode]* res4[inode];
               (*Q5i)[inode] = Qn_5[inode] -alpha_RK5[RK_index] * Delta_t[inode] /(*Vol)[inode]* res5[inode];
           }
           if(RK_index <= 4)
           {
               //computing residuals
               Calc_Residuals(nnodes, nedge,                                                                                                      // input variables
                                 e2n,                                                                                                              // input variables
                                 bn_edge_to_node_map,bn_size,                                                                                 // input variables
                                 Sx_i, Sy_i, Sz_i, S_i,                                                              // input variables
                                 Sx_b, Sy_b, Sz_b, S_b,                                                              // input variables
                                 Q1i, Q2i, Q3i, Q4i, Q5i,                                                   // input variables (conservatives vars for interior nodes)
                                 Q1b, Q2b, Q3b, Q4b, Q5b,                                                   // input variables (conservatives vars for boundary nodes)
                                 Cv, Cp, R, gamma,                                                                       // input vars (fluid characteristics)
                                 &Wx_e, &Wy_e, &Wz_e,                                    //weight for gradient prestored for each edge
                                 xyz, &u1, &u2, &u3, &u4, &u5, &grad_u1, &grad_u2, &grad_u3, &grad_u4, &grad_u5,  //preallocated gradient-related arrays
                                 &res1, &res2, &res3, &res4, &res5,&lambda0i_S,                          // output var, residuals at all nodes, lambda * S
                                 &lambda1_ghost, &lambda2_ghost, &lambda3_ghost, &lambda4_ghost, &lambda5_ghost); // output var, lambdai at ghost edges
           }
                   //impose viscous wall BCs
                   for(iedge = 1; iedge <= (*bn_size); iedge++)
                       if((*tag_b)[iedge] == 3) //viscous wall
                       {
                           //find the interior node due to that ghost edge
                           node1 = (*bn_edge_to_node_map)[iedge];
                           //setting viscous residual for that node
                           res2[node1] = 0.; //u=0 du/dt = 0
                           res3[node1] = 0.; //v=0 dv/dt = 0
                           res4[node1] = 0.; //w=0 dw/dt = 0
                           res5[node1] = res1[node1] * (*Cv) * T_inf;
                       }

       } //end of RK5 loop


       //After RK operations: Applying boundary conditions
       Impose_Boundary_Conditions(bn_edge_to_node_map,bn_size, tag_b,                                           // input vars - boundary spec
                               Sx_b, Sy_b, Sz_b,                                                   // input boundary normals
                               &lambda1_ghost, &lambda2_ghost, &lambda3_ghost, &lambda4_ghost, &lambda5_ghost, // input var, lambdai at ghost edges
                               Q1i, Q2i, Q3i, Q4i, Q5i,                          // input the value of conservative variables at interior nodes having ghost edges
                               qinf, Cv, Cp, R, gamma,                               // input primitive variables at infinity [rho_inf, u_inf, v_inf, w_inf, p_inf] and flow properties
                               Q1b, Q2b, Q3b, Q4b, Q5b);                          // input/output the value of conservative variables at ghost points

       //(MINOR): calculating norms if message displaying is on
       if (*messages_on)
       {
           for(inode = 1; inode <= (*nnodes); inode++)
           {
               L2norm[1] += pow(((*Q1i)[inode] - Qn_1[inode]), 2.);
               L2norm[2] += pow(((*Q2i)[inode] - Qn_2[inode]), 2.);
               L2norm[3] += pow(((*Q3i)[inode] - Qn_3[inode]), 2.);
               L2norm[4] += pow(((*Q4i)[inode] - Qn_4[inode]), 2.);
               L2norm[5] += pow(((*Q5i)[inode] - Qn_5[inode]), 2.);
           }
            L2norm[1] = sqrt(L2norm[1]/(*nnodes));
            L2norm[2] = sqrt(L2norm[2]/(*nnodes));
            L2norm[3] = sqrt(L2norm[3]/(*nnodes));
            L2norm[4] = sqrt(L2norm[4]/(*nnodes));
            L2norm[5] = sqrt(L2norm[5]/(*nnodes));

            printf("Itr %d -> L2_res1 = %e, L2_res2 = %e, L2_res3 = %e, L2_res4 = %e, L2_res5 = %e, t/t_ref = %e,\n\r", i,  L2norm[1],  L2norm[2],  L2norm[3],  L2norm[4], L2norm[5], (time / (double)time_scale));
            //printf("Q1i = %.24e", (*Q1i)[10]);

       }

     //writing residuals to the output file
     fprintf(Res_out,"%d %e %e %e %e %e\n\r",i, L2norm[1],  L2norm[2],  L2norm[3],  L2norm[4], L2norm[5]); //writing the pressure versus time at one point on the inlet

     //resseting output filename
     out_contours[1] = '\0';
     if(*ITR_Anim) //animation maker based on iterations
         if((i % (*ITR_Anim)) == 0)
         {
             sprintf(out_contours, "%s%d","ITR",i);
             Export_Sol_To_VTK(out_contours, Q1i, Q2i, Q3i, Q4i, Q5i, nnodes, c2n, xyz, nelem);
         }

     //animation maker based on simulation time for time-accurate solutions
     if(*time_accur_on)
     {
         time += deltat;
         if((time / (double)time_scale) >= export_intervals)
         {
             sprintf(out_contours, "%s%d","vortex",export_intervals);
             export_intervals++;
             Export_Sol_To_VTK(out_contours, Q1i, Q2i, Q3i, Q4i, Q5i, nnodes, c2n, xyz, nelem);
         }
     }

   } //end of main iteration loop


    //closing the output file
    fclose(Res_out);

    //clean-up procedures
    free(res1);
    free(res2);
    free(res3);
    free(res4);
    free(res5);
    free(Wx_e); free(Wy_e); free(Wz_e);
    free(u1); free(u2); free(u3); free(u4); free(u5);
    free(grad_u1); free(grad_u2); free(grad_u3); free(grad_u4); free(grad_u5); //each 1 .. 3*(nnodes+1)
    free(lambda0i_S);

   //completed successfully.
    return 0;
}

//This function exports the solution to VTK format so we can view it with Visit.
int Export_Sol_To_VTK(char *filenom,double **Q1i, double **Q2i, double **Q3i, double **Q4i, double **Q5i,int *nnodes, int*c2n[], const double *xyz, int nelem[])
{
    //introducing local variables
    int i;
    //converting Q1 ... 5 to Q[]
    double *solution = (double *) calloc(5*(*nnodes + 1), sizeof(double));

    for (i = 1; i <= (*nnodes); i++)
        {
          solution[5*i + 0] = (*Q1i)[i];
          solution[5*i + 1] = (*Q2i)[i];
          solution[5*i + 2] = (*Q3i)[i];
          solution[5*i + 3] = (*Q4i)[i];
          solution[5*i + 4] = (*Q5i)[i];

        }


    VTK_Grid_Out_with_solution((char *)filenom, *nnodes, xyz, nelem, c2n, solution, 5);

    //clean-up
    free(solution);

    //completed successfully!
    return 0;
}


//This function calculates the weight function Wx_e, Wy_e, and Wz_e which will be used to compute gradients based on second-order least-square method.
//Each We is a 2*(nedge +1) array with the following addressing. For example for Wx_e we have
// Wx_e(2*i + 0) = weight for point 1 on the edge (i = 1 ... nedge)
// Wx_e(2*i + 1) = weight for point 2 on the edge (i = 1 ... nedge)
//Similar for other weights
int Calculate_W_e(int *nnodes,                //input- number of nodes
                         const double *xyz,      //input- coordinates
                         int *nedge,             //input- number of edges
                         int **e2n,              //input- edge to node map
                         double **Wx_e, double **Wy_e, double **Wz_e) //output- weight function used to compute gradients (edge-based)
{
    //introducing local variables
    //indices for node and edge loops
    int inode = 0, iedge = 0;

    //data structure for two points forming an edge
    int point1, point2;
    double x1,y1,z1;
    double x2,y2,z2;


    //local arrays definition
    double *deltaxe, *deltaye, *deltaze; //edge-based
    double *s11, *s12, *s13, *s22, *s23, *s33; //node-based
    double *r11, *r12, *r13, *r22, *r23, *r33, *r22_0, *r33_0; //node-based

    //allocating local and output variables
    s11 = (double *)calloc((*nnodes + 1), sizeof(double));
    s12 = (double *)calloc((*nnodes + 1), sizeof(double));
    s13 = (double *)calloc((*nnodes + 1), sizeof(double));
    s22 = (double *)calloc((*nnodes + 1), sizeof(double));
    s23 = (double *)calloc((*nnodes + 1), sizeof(double));
    s33 = (double *)calloc((*nnodes + 1), sizeof(double));

    //resetting sij arrays to zero
    for(inode = 1; inode <= (*nnodes); inode++)
    {
        s11[inode] = 0.0;
        s12[inode] = 0.0;
        s13[inode] = 0.0;
        s22[inode] = 0.0;
        s23[inode] = 0.0;
        s33[inode] = 0.0;
    }

    r11 = (double *)calloc((*nnodes + 1), sizeof(double));
    r12 = (double *)calloc((*nnodes + 1), sizeof(double));
    r13 = (double *)calloc((*nnodes + 1), sizeof(double));
    r22 = (double *)calloc((*nnodes + 1), sizeof(double));
    r23 = (double *)calloc((*nnodes + 1), sizeof(double));
    r33 = (double *)calloc((*nnodes + 1), sizeof(double));
    r22_0 = (double *)calloc((*nnodes + 1), sizeof(double));
    r33_0 = (double *)calloc((*nnodes + 1), sizeof(double));

    //allocating edge-based arrays
    deltaxe = (double *)calloc((*nedge + 1), sizeof(double));
    deltaye = (double *)calloc((*nedge + 1), sizeof(double));
    deltaze = (double *)calloc((*nedge + 1), sizeof(double));

    //allocating output weights
    (*Wx_e) = (double *)calloc(2*(*nedge + 1), sizeof(double));
    (*Wy_e) = (double *)calloc(2*(*nedge + 1), sizeof(double));
    (*Wz_e) = (double *)calloc(2*(*nedge + 1), sizeof(double));


    //STEP I: finding sij arrays
    for(iedge = 1; iedge <= (*nedge); iedge++)
    {

        //step 1: calculating deltae
        //NOTE: WE MADE THE ASSUMPTION THAT alpha_je = 1
        //find points forming that edge
        point1 = (*e2n)[2*iedge + 0];
        point2 = (*e2n)[2*iedge + 1];

        x1 = xyz[3*point1 + 0]; //getting the x-coor of point1
        y1 = xyz[3*point1 + 1]; //getting the y-coor of point1
        z1 = xyz[3*point1 + 2]; //getting the z-coor of point1

        x2 = xyz[3*point2 + 0]; //getting the x-coor of point2
        y2 = xyz[3*point2 + 1]; //getting the y-coor of point2
        z2 = xyz[3*point2 + 2]; //getting the z-coor of point2

        //computing deltae components
        deltaxe[iedge] = (x2- x1);
        deltaye[iedge] = (y2- y1);
        deltaze[iedge] = (z2- z1);

        //step 2: calculating sij arrays
        s11[point1] += (deltaxe[iedge] * deltaxe[iedge]);
        s11[point2] += (-deltaxe[iedge] * -deltaxe[iedge]);

        s12[point1] += (deltaxe[iedge] * deltaye[iedge]);
        s12[point2] += (-deltaxe[iedge] * -deltaye[iedge]);

        s13[point1] += (deltaxe[iedge] * deltaze[iedge]);
        s13[point2] += (-deltaxe[iedge] * -deltaze[iedge]);

        s22[point1] += (deltaye[iedge] * deltaye[iedge]);
        s22[point2] += (-deltaye[iedge] * -deltaye[iedge]);

        s23[point1] += (deltaye[iedge] * deltaze[iedge]);
        s23[point2] += (-deltaye[iedge] * -deltaze[iedge]);

        s33[point1] += (deltaze[iedge] * deltaze[iedge]);
        s33[point2] += (-deltaze[iedge] * -deltaze[iedge]);

    } //

    //STEP II: calculating rij (node-based)
    for(inode = 1; inode <= (*nnodes); inode++)
    {
        r11[inode] = s11[inode];
        r12[inode] = s12[inode];
        r13[inode] = s13[inode];

        r22[inode] = s22[inode] - r12[inode] * r12[inode] / r11[inode];
        r23[inode] = s23[inode] - r12[inode] * r13[inode] / r11[inode];

        r33[inode] = s33[inode] - r13[inode] * r13[inode] / r11[inode] - r23[inode] * r23[inode] / r22[inode];

        r22_0[inode] = r22[inode];

        r33_0[inode] = r33[inode];
    }

    //STEP III: calculating weight functions per edges (Wx_e,Wy_e,Wz_e)
    for(iedge = 1; iedge <= (*nedge); iedge++)
    {
        //step 1: calculating point indices forming that edge
        point1 = (*e2n)[2*iedge + 0];
        point2 = (*e2n)[2*iedge + 1];

        // calculating Wz_e for point 1 on that edge
        (*Wz_e)[2*iedge + 0] = 1.0/r33_0[point1] * (deltaze[iedge]- r13[point1]/ r11[point1] * deltaxe[iedge] - r23[point1]/ r22[point1] * (deltaye[iedge] - r12[point1]/ r11[point1] * deltaxe[iedge]) );

        // calculating Wz_e for point 2 on that edge
        (*Wz_e)[2*iedge + 1] = 1.0/r33_0[point2] * (-deltaze[iedge]- r13[point2]/ r11[point2] * (-deltaxe[iedge]) - r23[point2]/ r22[point2] * (-deltaye[iedge] - r12[point2]/ r11[point2] * (-deltaxe[iedge]) ) );

        // calculating Wy_e for point 1 on that edge
        (*Wy_e)[2*iedge + 0] = 1.0/r22_0[point1] * (deltaye[iedge]- r12[point1]/ r11[point1] * deltaxe[iedge] - r23[point1]* (*Wz_e)[2*iedge + 0] );

        // calculating Wy_e for point 2 on that edge
        (*Wy_e)[2*iedge + 1] = 1.0/r22_0[point2] * (-deltaye[iedge]- r12[point2]/ r11[point2] * (-deltaxe[iedge]) - r23[point2]* (*Wz_e)[2*iedge + 1] );

        // calculating Wx_e for point 1 on that edge
        (*Wx_e)[2*iedge + 0] = 1.0/r11[point1] * (deltaxe[iedge]- r12[point1] * (*Wy_e)[2*iedge + 0] - r13[point1]* (*Wz_e)[2*iedge + 0] );

        // calculating Wx_e for point 2 on that edge
        (*Wx_e)[2*iedge + 1] = 1.0/r11[point2] * (-deltaxe[iedge]- r12[point2] * (*Wy_e)[2*iedge + 1] - r13[point2]* (*Wz_e)[2*iedge + 1] );

    } //end We loop

    //clean-up
    free(deltaxe);
    free(deltaye);
    free(deltaze); //edge-based
    free(s11);
    free(s12);
    free(s13);
    free(s22);
    free(s23);
    free(s33); //node-based
    free(r11);
    free(r12);
    free(r13);
    free(r22); free(r23); free(r33); free(r22_0); free(r33_0); //node-based


    //completed successfully!
    return 0;
}

//This function calculates the gradients of five dummy vectors u1..5 over unstructured grid
//Each grad is a 3*(nnodes +1) array with the following addressing. For example for grad_u_1 we have
// grad_u_1(3*i + 0) = partialu1/partialx for point i
// grad_u_1(3*i + 1) = partialu1/partialy for point i
// grad_u_1(3*i + 2) = partialu1/partialz for point i
//all variable ARE PREALLOCATED TO IMPROVE PERFORMANCE DURING THE MAIN EDGE LOOP
int Calc_Grad_u(    int *nnodes,                 //input - number of node
                         int *nedge,             //input- number of edges
                         int **e2n,              //input- edge to node map
                         double **Wx_e, double **Wy_e, double **Wz_e, //input- weight function used to compute gradients (edge-based)
                         double **u_1, double **u_2, double **u_3, double **u_4, double **u_5, //input u vectors
                         double **grad_u_1, double **grad_u_2, double **grad_u_3, double **grad_u_4, double **grad_u_5)   //output grad(u) vectors
{
    //declaring local variables
    int inode, iedge, point1, point2;
    int i = 0;

    //resetting all output gradients
    for(inode = 1; inode <= (*nnodes); inode++)
        for(i = 0; i<3; i++)
        {
            (*grad_u_1)[3*inode + i] = 0.;
            (*grad_u_2)[3*inode + i] = 0.;
            (*grad_u_3)[3*inode + i] = 0.;
            (*grad_u_4)[3*inode + i] = 0.;
            (*grad_u_5)[3*inode + i] = 0.;
        }

    //looping over all edges in the domain
    for( iedge = 1; iedge <= (*nedge); iedge++)
    {
        //Finding two points of each edge.
        point1 = (*e2n)[2*iedge+ 0];
        point2 = (*e2n)[2*iedge+ 1];

        //calculating grad_u_1 and ditributing among two nodes 
        (*grad_u_1)[3*point1 + 0] += (*Wx_e)[2*iedge + 0] * ((*u_1)[point2] - (*u_1)[point1]);
        (*grad_u_1)[3*point1 + 1] += (*Wy_e)[2*iedge + 0] * ((*u_1)[point2] - (*u_1)[point1]);
        (*grad_u_1)[3*point1 + 2] += (*Wz_e)[2*iedge + 0] * ((*u_1)[point2] - (*u_1)[point1]);

        (*grad_u_1)[3*point2 + 0] += (*Wx_e)[2*iedge + 1] * ((*u_1)[point1] - (*u_1)[point2]);
        (*grad_u_1)[3*point2 + 1] += (*Wy_e)[2*iedge + 1] * ((*u_1)[point1] - (*u_1)[point2]);
        (*grad_u_1)[3*point2 + 2] += (*Wz_e)[2*iedge + 1] * ((*u_1)[point1] - (*u_1)[point2]);

        //calculating grad_u_2 and ditributing among two nodes
        (*grad_u_2)[3*point1 + 0] += (*Wx_e)[2*iedge + 0] * ((*u_2)[point2] - (*u_2)[point1]);
        (*grad_u_2)[3*point1 + 1] += (*Wy_e)[2*iedge + 0] * ((*u_2)[point2] - (*u_2)[point1]);
        (*grad_u_2)[3*point1 + 2] += (*Wz_e)[2*iedge + 0] * ((*u_2)[point2] - (*u_2)[point1]);

        (*grad_u_2)[3*point2 + 0] += (*Wx_e)[2*iedge + 1] * ((*u_2)[point1] - (*u_2)[point2]);
        (*grad_u_2)[3*point2 + 1] += (*Wy_e)[2*iedge + 1] * ((*u_2)[point1] - (*u_2)[point2]);
        (*grad_u_2)[3*point2 + 2] += (*Wz_e)[2*iedge + 1] * ((*u_2)[point1] - (*u_2)[point2]);

        //calculating grad_u_3 and ditributing among two nodes
        (*grad_u_3)[3*point1 + 0] += (*Wx_e)[2*iedge + 0] * ((*u_3)[point2] - (*u_3)[point1]);
        (*grad_u_3)[3*point1 + 1] += (*Wy_e)[2*iedge + 0] * ((*u_3)[point2] - (*u_3)[point1]);
        (*grad_u_3)[3*point1 + 2] += (*Wz_e)[2*iedge + 0] * ((*u_3)[point2] - (*u_3)[point1]);

        (*grad_u_3)[3*point2 + 0] += (*Wx_e)[2*iedge + 1] * ((*u_3)[point1] - (*u_3)[point2]);
        (*grad_u_3)[3*point2 + 1] += (*Wy_e)[2*iedge + 1] * ((*u_3)[point1] - (*u_3)[point2]);
        (*grad_u_3)[3*point2 + 2] += (*Wz_e)[2*iedge + 1] * ((*u_3)[point1] - (*u_3)[point2]);

        //calculating grad_u_4 and ditributing among two nodes
        (*grad_u_4)[3*point1 + 0] += (*Wx_e)[2*iedge + 0] * ((*u_4)[point2] - (*u_4)[point1]);
        (*grad_u_4)[3*point1 + 1] += (*Wy_e)[2*iedge + 0] * ((*u_4)[point2] - (*u_4)[point1]);
        (*grad_u_4)[3*point1 + 2] += (*Wz_e)[2*iedge + 0] * ((*u_4)[point2] - (*u_4)[point1]);

        (*grad_u_4)[3*point2 + 0] += (*Wx_e)[2*iedge + 1] * ((*u_4)[point1] - (*u_4)[point2]);
        (*grad_u_4)[3*point2 + 1] += (*Wy_e)[2*iedge + 1] * ((*u_4)[point1] - (*u_4)[point2]);
        (*grad_u_4)[3*point2 + 2] += (*Wz_e)[2*iedge + 1] * ((*u_4)[point1] - (*u_4)[point2]);

        //calculating grad_u_5 and ditributing among two nodes
        (*grad_u_5)[3*point1 + 0] += (*Wx_e)[2*iedge + 0] * ((*u_5)[point2] - (*u_5)[point1]);
        (*grad_u_5)[3*point1 + 1] += (*Wy_e)[2*iedge + 0] * ((*u_5)[point2] - (*u_5)[point1]);
        (*grad_u_5)[3*point1 + 2] += (*Wz_e)[2*iedge + 0] * ((*u_5)[point2] - (*u_5)[point1]);

        (*grad_u_5)[3*point2 + 0] += (*Wx_e)[2*iedge + 1] * ((*u_5)[point1] - (*u_5)[point2]);
        (*grad_u_5)[3*point2 + 1] += (*Wy_e)[2*iedge + 1] * ((*u_5)[point1] - (*u_5)[point2]);
        (*grad_u_5)[3*point2 + 2] += (*Wz_e)[2*iedge + 1] * ((*u_5)[point1] - (*u_5)[point2]);

    }



    //completed successfully!
    return 0;
}

//this function test the validity of the gradient calculation
int Grad_Benchmarks(int *nnodes,                 //input- number of nodes
                         const double *xyz,      //input- coordinates
                         int *nedge,             //input- number of edges
                         int **e2n)              //input- edge to node map
{
 
    //declaring local variables
    int inode, i;
    double xi,yi,zi;
    double *Wx_e, *Wy_e, *Wz_e; //preallocated
    double *u_1, *u_2, *u_3, *u_4, *u_5;
    double *grad_u_1, *grad_u_2, *grad_u_3, *grad_u_4, *grad_u_5;
    
    
    //Allocating local vars
    u_1 =  (double *)calloc((*nnodes + 1), sizeof(double));
    u_2 =  (double *)calloc((*nnodes + 1), sizeof(double));
    u_3 =  (double *)calloc((*nnodes + 1), sizeof(double));
    u_4 =  (double *)calloc((*nnodes + 1), sizeof(double));
    u_5 =  (double *)calloc((*nnodes + 1), sizeof(double));
    grad_u_1 =  (double *)calloc(3*(*nnodes + 1), sizeof(double));
    grad_u_2 =  (double *)calloc(3*(*nnodes + 1), sizeof(double));
    grad_u_3 =  (double *)calloc(3*(*nnodes + 1), sizeof(double));
    grad_u_4 =  (double *)calloc(3*(*nnodes + 1), sizeof(double));
    grad_u_5 =  (double *)calloc(3*(*nnodes + 1), sizeof(double));

    
    //initializing u vectors
    for(inode = 1; inode <= (*nnodes); inode++)
    {
        xi = xyz[3*inode + 0]; //getting the x-coor 
        yi = xyz[3*inode + 1]; //getting the y-coor 
        zi = xyz[3*inode + 2]; //getting the z-coor 
        u_1[inode] = 1.* xi + 2. * yi + 3.* zi;
        u_2[inode] = 1.* xi + 2. * yi + 3.* zi;
        u_3[inode] = 1.* xi + 2. * yi + 3.* zi;
        u_4[inode] = 1.* xi + 2. * yi + 3.* zi;
        u_5[inode] = 1.* xi + 2. * yi + 3.* zi;
    }
    
/*
    //resetting grad arrays
    for(inode = 1; inode <= (*nnodes); inode++)
        for(i = 0; i<3; i++)
        {
                grad_u_1[3*inode + i] = 0.;
                grad_u_2[3*inode + i] = 0.;
                grad_u_3[3*inode + i] = 0.;
                grad_u_4[3*inode + i] = 0.;
                grad_u_5[3*inode + i] = 0.;
        }
*/
    
    
    Calculate_W_e(nnodes,                         //input- number of nodes
                        xyz,                      //input- coordinates
                         nedge,                   //input- number of edges
                         e2n,                     //input- edge to node map
                         &Wx_e, &Wy_e, &Wz_e);    //output- weight function used to compute gradients (edge-based)
    
    Calc_Grad_u(nnodes,                 //input - number of node
                         nedge,             //input- number of edges
                          e2n,              //input- edge to node map
                         &Wx_e, &Wy_e, &Wz_e, //input- weight function used to compute gradients (edge-based)
                         &u_1,&u_2, &u_3, &u_4, &u_5, //input u vectors
                         &grad_u_1, &grad_u_2, &grad_u_3, &grad_u_4, &grad_u_5);   //output grad(u) vectors

    printf("For point 2 => dx = %2.20f, dy = %2.20f, dz = %2.20f\n\r", grad_u_1[3*2 + 0], grad_u_1[3*2 + 1], grad_u_1[3*2 + 2]);

    //cleanup
    free(u_1);
    free(u_2);
    free(u_3);
    free(u_4);
    free(u_5);
    free(grad_u_1);
    free(grad_u_2);
    free(grad_u_3);
    free(grad_u_4);
    free(grad_u_5);

    //completed successfully!
    return 0;
}
///////////////////////////////////////////////////////////////////////////////
//////////END OF TEST SEGMENT OF THE ORIGINAL CODE/////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


int main()
{
  // Variables needed for reading in a grid

  int nnodes;
  int nelem[NUMBER_OF_ELEM_TYPES];
  int *c2n[NUMBER_OF_ELEM_TYPES] = {NULL, NULL, NULL, NULL, NULL, NULL};
  int *factag[NUMBER_OF_ELEM_TYPES] = {NULL, NULL, NULL, NULL, NULL, NULL};
  double *xyz = NULL;

  FILE *fp_in = NULL;
  int i;

  // Grid to read in

  fp_in = fopen("../grids/blasius_fine.ugrid", "r");
  
  // Read in the grid (ugrid format)

  Read_asciiugrid_Grid(fp_in, &nnodes, nelem, &xyz, c2n, factag);

  fclose(fp_in);

///////////////////////////////////////////////////////////////////////////////
//Beginning OF TEST SEGMENT FOR EVALUATING FUNCTIONS ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

  int node_num = 0, indx1 = 0, indx2 = 0, indx = 0;
  int *nesp;
  int *esp;
  
  //Finding elements near the given point
  Gen_Map_Elem__Sur_Point(&nnodes,nelem,c2n,&nesp,&esp);

  /* printf("Please enter the node number: "); */
  /* scanf("%d", &node_num); */
  node_num = 1;

  indx1 = nesp[node_num + 0];
  indx2 = nesp[node_num + 1];
  printf("\n\rHere is the list:\n\r");
  for(indx = indx1; indx < indx2; indx++)
      printf("\t%d\r\n",esp[indx]);
  printf("Completed.\n\r");
  
  //psp is a linked-list-array structure to hold points and their neighbors
  struct int_vec *psp;
  //Getting the neighbors
  Gen_Map_Point__Sur_Point(&nnodes,nelem,c2n,&nesp,&esp,&psp);

  //Echoing the neighbor points
  //for(indx = 1; indx <= psp->length; indx++)
  //{
  indx = node_num;
      printf("\r\nSpecification of point %d:\r\n",indx);
      printf("\tNumber of surrounding points: %d\r\n",psp->nexts[indx].length);
      printf("\tlist of surrounding points:\r\n");
      for(indx1 = 0; indx1 < psp->nexts[indx].length; indx1++)
          printf("\t\t%d\r\n",psp->nexts[indx].int_data[indx1]);
  //}
  printf("\r\n");


////////CRS version ////////////
/*  int ch;
  //scanf("%d",&ch);
  int *npsp1;
  int *psp1;
  int kkk, mmm;

  Gen_Map_Point__Sur_Point_CRS(&nnodes,nelem,c2n,&nesp,&esp,&npsp1,&psp1);

  for(kkk = 1; kkk <= nnodes; kkk++)
  {
      printf("\r\nSpecification of point %d:\r\n",kkk);
      printf("\tNumber of surrounding points: %d\r\n",(npsp1[kkk + 1] - npsp1[kkk + 0]));
      printf("\tlist of surrounding points:\r\n");
      for(mmm = npsp1[kkk + 0]; mmm < npsp1[kkk + 1]; mmm++)
      {
          printf("\t\t%d\r\n",psp1[mmm]);

      }
  }

  printf("\r\n");
 */
  //scanf("%d",&ch);
///////////////////////////////

  //Allocating memory for edges
  int nedge;
  int *e2n;
  //Getting surrounding edges
  Gen_Edg_Point(&psp,&nedge, &e2n);


  //Echoing the edges
  printf("Here is the list of edges:\r\n");
  for(indx = 1; indx <= nedge; indx++)
      printf("\r\n\t%d-%d",e2n[2*indx+0],e2n[2*indx+1]);

  printf("\r\n");


  //Getting surrounding edges surrounding a particular point
  struct int_vec *edsp;
  Create_Edges_Surr_Point(&nnodes, &nedge, &e2n, &edsp);
  for(indx1 = 1; indx1<= nnodes; indx1++)
  {
      printf("edges surrounding points %d:\n\r", indx1);
      for(indx2 = 1; indx2 <= edsp->nexts[indx1].length; indx2++)
          printf("\t%d, %d\n\r", e2n[2*edsp->nexts[indx1].int_data[indx2]+0], e2n[2*edsp->nexts[indx1].int_data[indx2]+1]);
  }

  //////////////////////////////////////////////////
  struct int_vec *ese;
  int ii, jj;

  Gen_Map_Elem__Sur_Elem(nelem,c2n,&nesp,&esp,&ese);
  for(ii = 1; ii <= ese->length; ii++)
  {
 // ii = 100 + node_num;
      printf("elements near element (%d):\n\r",ii);
      printf("\t");
      for(jj = 1; jj<= ese->nexts[ii].length;jj++)
          printf("%d,",ese->nexts[ii].int_data[jj]);
      printf("\n\r");
  }

  //examines the Create_Surr_Face_El_Node(struct int_vec **surr_face_el_node)
  //function by printing all possibilities.
  struct int_vec *surr_face_el_node;
  Create_Surr_Face_El_Node(&surr_face_el_node);
  int indx3;

  for(indx1 = Triangle; indx1 <= Hex; indx1++)
  {
      printf("For element type (%d), we have:\n\r",indx1);
      for(indx2 =1; indx2 < surr_face_el_node->nexts[indx1].length; indx2++)
      {
          printf("\tFor node (%d) we have the following surface tags:\n\r",indx2);
          for(indx3 = 1; indx3 < surr_face_el_node->nexts[indx1].nexts[indx2].length;indx3++)
              printf("\t\t%d\n\r",surr_face_el_node->nexts[indx1].nexts[indx2].int_data[indx3]);
      }
  }

  //examines the Create_Points_In_Face_El(struct int_vec **pnts_face_el)
  //function by printing all possibilities.
  struct int_vec *pnts_face_el;
  Create_Points_In_Face_El(&pnts_face_el);

  for(indx1 = Triangle; indx1 <= Hex; indx1++)
  {
      printf("For element type (%d), we have:\n\r",indx1);
      for(indx2 =1; indx2 < pnts_face_el->nexts[indx1].length; indx2++)
      {
          printf("\tFor face (%d) we have the following local nodes:\n\r",indx2);
          for(indx3 = 1; indx3 < pnts_face_el->nexts[indx1].nexts[indx2].length;indx3++)
              printf("\t\t%d\n\r",pnts_face_el->nexts[indx1].nexts[indx2].int_data[indx3]);
      }
  }

  //Examining function obtain_common_ints which is heavily used for obtaining
  //the number of common points in two integer sets
  int arr1[7] = {0,-1,4,2,89,6,10};
  int arr2[5] = {0,89,6,-1,10};
  int N1= 6, N2 = 4;
  int result = 0;
  obtain_common_ints(arr1, &N1, arr2, &N2, &result);

  printf("The number of common integers is: %d\n\r", result);

  struct double_vec *node_centered;
  struct double_vec *bn_struct;

  Create_Node_Centered_Grid(&nnodes,nelem,
			 &nesp,
                         &esp,
			 c2n,
                         xyz,
                         &nedge,
                         &e2n,
                         &node_centered,
                         &bn_struct);
  int inode, k;
  double voli = 0.0;
  double Areax = 0.0;
  double Areay = 0.0;
  double Areaz = 0.0;

/*

  for(inode = 1; inode <= nnodes; inode++)
      for(k =1; k<= node_centered->nexts[inode].length; k++)
          if(node_centered->nexts[inode].nexts[k].db_data[1] == -1.0)
              continue;
          else
              voli += node_centered->nexts[inode].nexts[k].db_data[1];
*/

  //printf("total vol = %f\n\r", voli);

  
  for(inode = 1; inode <= nnodes; inode++)
  {
      //    Areax = 0.;
      //    Areay = 0.;
      //    Areaz = 0.;
      //    voli = 0.;

      for(k =1; k<= node_centered->nexts[inode].length; k++)
      {
          Areax += node_centered->nexts[inode].nexts[k].db_data[3];
          Areay += node_centered->nexts[inode].nexts[k].db_data[4];
          Areaz += node_centered->nexts[inode].nexts[k].db_data[5];
          if( node_centered->nexts[inode].nexts[k].db_data[1] == -1.0)
              continue;
          else
              voli += node_centered->nexts[inode].nexts[k].db_data[1];
          //voli += node_centered->nexts[inode].nexts[k].db_data[1];
          //printf("k= %d, Ax = %e,Ay = %e,Az = %e, Dv=%e,\n\r",k, Areax, Areay,Areaz,voli);
      }
      //printf("\n\r NODE (%d) , Ax= %e, Ay= %e, Az=%e, Dv= %e\n\r",inode, Areax,Areay,Areaz,voli);
  }
printf("\n\r Ax= %e, Ay= %e, Az=%e, Voli=%e\n\r",Areax,Areay,Areaz, voli);

//doing it for boundary elements
int ibelem = 0;
Areax = 0.0;
Areay = 0.0;
Areaz = 0.0;

for(ibelem = 1; ibelem <= nelem[Triangle]; ibelem++)
      for(k =1; k<= 3; k++)
      {
          Areax +=  (bn_struct)->nexts[ibelem].nexts[k].db_data[1];
          Areay +=  (bn_struct)->nexts[ibelem].nexts[k].db_data[2];
          Areaz +=  (bn_struct)->nexts[ibelem].nexts[k].db_data[3];
      }

for(ibelem = (nelem[Triangle]+1); ibelem <= (nelem[Triangle] + nelem[Quad]); ibelem++)
      for(k =1; k<= 4; k++)
      {
          Areax +=  (bn_struct)->nexts[ibelem].nexts[k].db_data[1];
          Areay +=  (bn_struct)->nexts[ibelem].nexts[k].db_data[2];
          Areaz +=  (bn_struct)->nexts[ibelem].nexts[k].db_data[3];
      }

printf("\n\rat boundaries  Ax= %e, Ay= %e, Az=%e\n\r",Areax,Areay,Areaz);


//array-based conversion for flow solver
double *Vol, *Sx_i, *Sy_i, *Sz_i, *S_i;
double *Sx_b, *Sy_b, *Sz_b, *S_b;
int *bn_edge_to_node_map, *tag_b;
int bn_size;
int iedge;

Node_Based_To_Array_Convertor(&nnodes, &nedge,                             // input variables
                     &e2n, nelem, c2n,                                     // input variables
                     &edsp, factag,                                        // input variables
                     &node_centered,                                       // input variables
                     &bn_struct,                                           // input variables
                     &bn_edge_to_node_map, &Vol,                           // output variables
                     &Sx_i, &Sy_i, &Sz_i, &S_i,                            // output variables
                     &Sx_b, &Sy_b, &Sz_b, &S_b,                            // output variables
                     &bn_size, &tag_b);                                    // output variables

 
 //resetting measurment variables
  voli  = 0.0;
  Areax = 0.0;
  Areay = 0.0;
  Areaz = 0.0;
  iedge = 0;

  //
 for(inode = 1; inode <= nnodes; inode++)
     voli += Vol[inode];

 for(inode = 1; inode <= bn_size; inode++)
 {
     Areax +=  Sx_b[inode] * S_b[inode];
     Areay +=  Sy_b[inode] * S_b[inode];
     Areaz +=  Sz_b[inode] * S_b[inode];
     printf("boundary node %d -> Sx = %e, Sy = %e, Sz = %e, tag = %d, node = %d,\n\r",inode,
     Sx_b[inode] * S_b[inode],
     Sy_b[inode] * S_b[inode],
     Sz_b[inode] * S_b[inode], tag_b[inode], bn_edge_to_node_map[inode]);
 }

  /*
 for(iedge = 1; iedge <= nedge; iedge++)
 {
     printf("interior edge %d -> Sxi = %e, Syi = %e, Szi = %e\n\r",iedge,
     Sx_i[iedge] * S_i[iedge],
     Sy_i[iedge] * S_i[iedge],
     Sz_i[iedge] * S_i[iedge]);
 }
*/

  printf("The sum of area and volumes are (Ax= %e, Ay=%e, Az=%e) and V=%e respectively.\n\r",Areax, Areay, Areaz, voli);

 //this function test the validity of the gradient calculation
 Grad_Benchmarks(&nnodes,                 //input- number of nodes
                         xyz,      //input- coordinates
                         &nedge,             //input- number of edges
                         &e2n);              //input- edge to node map

/*
T = 20  + 273.1;
rho = 1.205;
P = 101e3;
gamma = 1.4;
R = 287.0;
mu = 15.11*1.0e-6 * rho;
Pr = 0.713;
Cp = 1005;
kappa = mu*Cp/Pr;
L = 0.32; //flat plate length
u = mu/(rho * L) * (5./0.03)^2
*/

 //////////////////////

 double R = 287.0, gamma = 1.4; //Nitrogen
 double Cp = gamma * R / (gamma -1.) , Cv = R/ (gamma -1.);
 double CFL = 1.2;

  int N_ITR =60000;
  short init_flag = 1, messages_on =1, time_accur_on = 0; int ITR_Anim = 500;
  double *Q1i, *Q2i, *Q3i,*Q4i,*Q5i;
  double *Q1b, *Q2b, *Q3b,*Q4b,*Q5b;

  double *qinf = (double *)calloc(6, sizeof(double));
  qinf[1] = 1.205;
  qinf[2] = 1.3116*40.;
  qinf[3] = 0.;
  qinf[4] = 0.;
  qinf[5] = 101300.;

  //Init_Tags(&nnodes, &bn_size, &tag_b, &bn_edge_to_node_map, &Q1i, &Q2i, &Q3i, &Q4i, &Q5i); // input/output variables (conservatives vars for interior nodes)


  Solve_RK5_Jameson(&nnodes, &nedge, &e2n, c2n, nelem,                                       // input vars - edge and node spec
                        &bn_edge_to_node_map,&bn_size,                                        // input vars - boundary spec
                         &Sx_i, &Sy_i, &Sz_i, &S_i, &Vol,  // input vars- volume
                         &Sx_b, &Sy_b, &Sz_b, &S_b,                // and boundary and interior normals
                         &Cv, &Cp, &R, &gamma,                         //input vars (fluid characteristics)
                         &CFL, &N_ITR, &init_flag, &messages_on, &time_accur_on, &ITR_Anim,            //input vars - solver parameters
                         &tag_b, &qinf,                                               //input vars - boundary element tags and free stream conds
                         xyz, //coordinates used for gradient method
                         &Q1i, &Q2i, &Q3i, &Q4i, &Q5i,     // output variables (conservatives vars for interior nodes)
                         &Q1b, &Q2b, &Q3b, &Q4b, &Q5b);     // output variables (conservatives vars for boundary nodes)           //input vars ()


  //Init_Tags(&nnodes, &bn_size, &tag_b, &bn_edge_to_node_map, &Q1i, &Q2i, &Q3i, &Q4i, &Q5i);

  /*

  Solve_Euler_Explicit(&nnodes,&nedge, &e2n,                                              // input vars - edge and node spec
                         &bn_node_number,&bn_size,                                        // input vars - boundary spec
                         &Sx_i, &Sy_i, &Sz_i, &S_i, &Vol,                                 // input vars- volume
                         &Sx_b, &Sy_b, &Sz_b, &S_b,                                       // and boundary and interior normals
                         &Cv, &Cp, &R, &gamma,                                            // input vars (fluid characteristics)
                         &CFL, &N_ITR, &init_flag, &messages_on,                       //input vars - solver parameters
                         &Q1i, &Q2i, &Q3i, &Q4i, &Q5i,                                    // output variables (conservatives vars for interior nodes)
                         &Q1b, &Q2b, &Q3b, &Q4b, &Q5b, &is_it_boundary_node);                                   // output variables (conservatives vars for boundary nodes)           //input vars ()
  //////////
*/


  Export_Sol_To_VTK("SteadyState", &Q1i, &Q2i, &Q3i, &Q4i, &Q5i, &nnodes, c2n, xyz, nelem);
/*
  /////////

 
  //Tester for matrix multipication ///
  double M[6][6], P[6][6], T[6][6];
  M[1][1] = 1.;
  M[1][2] = 2.;
  M[1][3] = 3.;
  M[1][4] = 4.;
  M[1][5] = 5.;
  M[2][1] = 6.;
  M[2][2] = 7.;
  M[2][3] = 8.;
  M[2][4] = 9.;
  M[2][5] = 10.;
  M[3][1] = 11.;
  M[3][2] = 12.;
  M[3][3] = 13.;
  M[3][4] = 14.;
  M[3][5] = 15.;
  M[4][1] = 16.;
  M[4][2] = 17.;
  M[4][3] = 18.;
  M[4][4] = 19.;
  M[4][5] = 20.;
  M[5][1] = 21.;
  M[5][2] = 22.;
  M[5][3] = 23.;
  M[5][4] = 24.;
  M[5][5] = 25.;

  P[1][1] = 1.;
  P[1][2] = -2.;
  P[1][3] = 3.;
  P[1][4] = 4.;
  P[1][5] = 5.;
  P[2][1] = 6.;
  P[2][2] = 7.;
  P[2][3] = 8.;
  P[2][4] = 9.;
  P[2][5] = 10.;
  P[3][1] = 11.;
  P[3][2] = 12.;
  P[3][3] = 13.;
  P[3][4] = 14.;
  P[3][5] = 15.;
  P[4][1] = 16.;
  P[4][2] = 17.;
  P[4][3] = 18.;
  P[4][4] = 19.;
  P[4][5] = 20.;
  P[5][1] = 21.;
  P[5][2] = 22.;
  P[5][3] = 23.;
  P[5][4] = .12;
  P[5][5] = 25.;

  Sq_Matrix_Product(M, P, T,5);
  int j =0;
  printf("\n\rMatrix multipication results:");

  for(i = 1; i<= 5; i++)
  { printf("\n\r");
      for(j=1; j<=5; j++)
          printf("%f\t",T[i][j]);
  }
  ///// end of matrix multipication tester
*/

  //clean-up for test section of the code
  free(e2n);
  free(psp);
  free(nesp);
  free(esp);
  printf("Cleanup at test-section of the code completed successfully!\n\r");

///////////////////////////////////////////////////////////////////////////////
//END OF TEST SEGMENT FOR EVALUATING FUNCTIONS ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

  // Write out the grid (Fieldview format)

//  Write_Fieldview_Grid("13-node.b4.fv", NULL, nnodes, nelem, xyz, c2n, factag,
//		       0, NULL, NULL, 0, -1, NULL, NULL, 0, NULL);

  // Free memory

  free(xyz);
  for (i = 0; i < NUMBER_OF_ELEM_TYPES; i++)
    {
      free(c2n[i]);
      free(factag[i]);
    }

  return 0;
}
