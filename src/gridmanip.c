#include "generic.h"

///////////////////////////////////////////////////////////////////////////////
///MY CODE AND WRITING
///THE FOLLOWING DECLERATIONS ARE NOT IN THE ORIGINAL CODE BUT ADDED BY ME
///////////////////////////////////////////////////////////////////////////////

int Local_Elem_To_Global_Elem(int ielem,
                              int nelem[],
                              int etype)
{
    //Declaring index variables
    int e;

    //Declaring output and initialization
    int gelem = ielem;

    for(e = Triangle; e < etype; e++)
        gelem += nelem[e];

    return gelem;

}

int Gen_Map_Elem__Sur_Point(int *nnodes,
			 int nelem[],
			 int*c2n[],
			 int **nesp,
                         int **esp
			 )
{
    //Declaring index variables
    int ii, ielem, inode, mnode, node, indx;
    ElemType etype;

    //Allocating memory for nesp
    *nesp = (int *) calloc((*nnodes) + 2, sizeof(int));

    //1 - Initializing nesp to ZERO
    for(ii = 0; ii < (*nnodes) + 2; ii++)
        *(*nesp + ii) = 0;

    //2- Filling in the entries for nesp
    for (etype = Triangle; etype <= Hex; etype++)
        for (ielem = 1; ielem <= nelem[etype]; ielem++)
        {
            mnode = ElemTypeToMnode(etype);
            for(inode = 0; inode < mnode; inode++)
            {
                node = c2n[etype][mnode*ielem + inode];
                (*nesp)[node + 1]++;
            }

        }
    //3- Generate the offset for the CRS array

    for(ii = 2; ii <= ((*nnodes) + 1); ii++)
        (*nesp)[ii] += (*nesp)[ii - 1];

    //4- Allocating memory for esp
    *esp = (int *) calloc((*nesp)[(*nnodes) + 1], sizeof(int));

    //5- Reapeating (2) , but storing the data
    for (etype = Triangle; etype <= Hex; etype++)
        for (ielem = 1; ielem <= nelem[etype]; ielem++)
        {
            mnode = ElemTypeToMnode(etype);
            for(inode = 0; inode < mnode; inode++)
            {
                node = c2n[etype][mnode*ielem + inode];
                indx = (*nesp)[node];
                (*esp)[indx] = Local_Elem_To_Global_Elem(ielem,nelem,etype);
                (*nesp)[node]++;
            }

        }

    //6- The CRS indices are shifted by 1 after Step (5)
    for( ii = ((*nnodes) + 1); ii >= 1; ii--)
        (*nesp)[ii] = (*nesp)[ii-1];


    //Operations completed successfully!
    return 0;

}

void Global_Elem_To_Local_Elem(int *ielem,
                              int nelem[],
                              int gelem,
                              int *etype)
{

    //Declaring index variables
    int e;

    //Initializing the local element index
    *ielem = gelem;
    for( e = Triangle; e <= Hex; e++)
        if((*ielem) <= nelem[e]) //if the element fits within a type
        {
            *etype = e; //everything is done!
            break;
        }
        else //otherwise we should subtract the previous element indices
            (*ielem) -= nelem[e];

}

int Gen_Map_Point__Sur_Point_LookUp(int *nnodes,
			 int nelem[],
			 int*c2n[],
			 int **nesp,
                         int **esp,
			 int your_point, int **surr_points, struct int_vec **enct)
{
    //Defining local variables
    int indx,indx1,indx2;
    int mnode,inode,node;
    int gelem, ielem , etype;
    int point_index;

    //For the given point, find all surrounding elements
    indx1 = (*nesp)[your_point + 0];
    indx2 = (*nesp)[your_point + 1];
    for(indx = indx1; indx < indx2; indx++)
    {
        gelem = (*esp)[indx]; //get the golobal index of the neighbur
        //find the local parameters of that element
        Global_Elem_To_Local_Elem(&ielem,nelem,gelem,&etype);
        mnode = ElemTypeToMnode(etype); //find the number of nodes in element
        //find the index of given node (one-based index) in that
        //element so that we can recognize neighbouring points
        for(inode = 0; inode < mnode; inode++)
            if( your_point == c2n[etype][mnode*ielem + inode])
                 point_index = inode+1;
        //detect appropriate surrounding nodes
        for(inode = 1; inode < (*enct)->nexts[etype].nexts[point_index].length; inode++)
        {
            //find the local number of surrounding points in that element
            node = (*enct)->nexts[etype].nexts[point_index].int_data[inode];
            //convert the the local node of line above to the corresponding
            //global number of node.
            node = c2n[etype][mnode*ielem + node-1];
            //put a tick in the position where surrounding node appeares
            (*surr_points)[node]++;
        }
    }

    return 0; //completed successfully
}

int Gen_Map_Point__Sur_Point(int *nnodes,
			 int nelem[],
			 int*c2n[],
			 int **nesp,
                         int **esp,
			 struct int_vec **psp)
{
    //declaring local variables
    int innodes = 0, j = 0, ngb_size = 0, k = 0;

    //decleration of element-node map holder pointer
    struct int_vec *enct;

    Create_Element_Node_Conn_Table(&enct);

    //initializing the ROOT of output vector containing neighbours for all
    //points in the domain.
    *psp = (struct int_vec *)calloc(1,sizeof(struct int_vec));
    (*psp)->length = *nnodes;
    (*psp)->nexts = (struct int_vec *)calloc((*nnodes+1),sizeof(struct int_vec));

    //Allocating surr_points array with the same size of *nnodes + 1
    int *surr_points = (int *)calloc(*nnodes + 1,sizeof(int));
    //Initializing surr_points to zero
    for(j = 1; j <= *nnodes; j++)
        *(surr_points+j) = 0;

    //Proceeding to determine surrounding points of all points
    for(innodes = 1; innodes <=  *nnodes; innodes++)
    {
        //Get the neighbour points for innodes
        Gen_Map_Point__Sur_Point_LookUp(nnodes,nelem,c2n,nesp,esp,innodes,&surr_points,&enct);

        //Calculate the number of neighbouring points
        for(j = 1; j <= *nnodes; j++)
            if(surr_points[j] != 0)
                ngb_size++;

        //Creating appropriate array with size "ngb_size" for storing the
        //neighburs of "innodes"
        (*psp)->nexts[innodes].int_data = (int *)calloc(ngb_size,sizeof(int));
        (*psp)->nexts[innodes].length = ngb_size;

        //Now, proceed to store neighbour points
        for(j = 1; j <= *nnodes; j++)
            if(surr_points[j] != 0)
            {
                (*psp)->nexts[innodes].int_data[k] = j;
                k++;
            }
        //resetting "k" for the next array
        k = 0;
        //resetting "ngb_size" for the next array
        ngb_size = 0;
        //resetting surr_points
        for(j = 1; j <= *nnodes; j++)
            surr_points[j] = 0;

    } //done for all points!

    //Performing clean-up for arrays allocated inside the function
    free(surr_points);
    free(enct);

return 0; //return successfully
}

int Gen_Map_Point__Sur_Point_CRS(int *nnodes,
			 int nelem[],
			 int *c2n[],
			 int **nesp,
                         int **esp,
			 int **npsp,
                         int **psp)
{
    //declaring local variables
    int innodes = 0, j = 0, ngb_size = 0, CRS_i = 0;
    //decleration of element-node map holder pointer
    struct int_vec *enct;

    Create_Element_Node_Conn_Table(&enct);

    //allocating memory for CRS array indexer
    *npsp = (int *)calloc(*nnodes+2, sizeof(int));

    //Initializing CRS indexer
    for(CRS_i = 0; CRS_i < (*nnodes+2); CRS_i++)
        (*npsp)[CRS_i] = 0;

    //Allocating surr_points array with the same size of *nnodes + 1
    int *surr_points = (int *)calloc(*nnodes + 1,sizeof(int));
    //Initializing surr_points to zero
    for(j = 1; j <= *nnodes; j++)
        *(surr_points+j) = 0;

    //(*) Proceeding to determine the number surrounding points of all points
    for(innodes = 1; innodes <=  *nnodes; innodes++)
    {

        //Get the neighbour points for innodes
        Gen_Map_Point__Sur_Point_LookUp(nnodes,nelem,c2n,nesp,esp,innodes,&surr_points,&enct);

        //Calculate the number of neighbouring points
        for(j = 1; j <= *nnodes; j++)
            if(surr_points[j] != 0)
                ngb_size++;

        (*npsp)[innodes+1] = ngb_size; //storing the number of surr points
        //for node "innodes" in the CRS indexer

        //resetting "ngb_size" for the next array
        ngb_size = 0;
        //resetting surr_points
        for(j = 1; j <= *nnodes; j++)
            surr_points[j] = 0;


    }
    //Generating the offset for CRS array
    for(CRS_i = 2; CRS_i <= (*nnodes+1); CRS_i++)
        (*npsp)[CRS_i] += (*npsp)[CRS_i-1];
    //Allocating CRS array for psp based on indexer
    *psp = (int *)calloc((*npsp)[*nnodes+1], sizeof(int));
    //Reapiting the operation in (*), but storing the points at this stage
    for(innodes = 1; innodes <=  *nnodes; innodes++)
    {
        //Get the neighbour points for innodes
        Gen_Map_Point__Sur_Point_LookUp(nnodes,nelem,c2n,nesp,esp,innodes,&surr_points,&enct);
        //Now, proceed to store neighbour points
        for(j = 1; j <= *nnodes; j++)
            if(surr_points[j] != 0)
            {
                (*psp)[(*npsp)[innodes]] = j;
                (*npsp)[innodes]++;
            }
        //resetting surr_points
        for(j = 1; j <= *nnodes; j++)
            surr_points[j] = 0;
    }
    //Shifting CRS indices by 1
    for(CRS_i = (*nnodes+1); CRS_i >= 1; CRS_i--)
        (*npsp)[CRS_i] = (*npsp)[CRS_i-1];

/*    //Printing the output
    for(CRS_i = 0; CRS_i < (*npsp)[*nnodes+1]; CRS_i++)
        printf("%d\n\r",(*psp)[CRS_i]);
*/

    //Performing clean-up for arrays allocated inside the function
    free(surr_points);
    free(enct);

return 0; //return successfully
}
int Gen_Edg_Point(struct int_vec **psp,int *nedge, int **e2n)
{
    //initializing local variables
    int inode = 0, jnode =0; //for indexing domain nodes
    int estim = *nedge = 0; //for measuring the number of edges

    //Estimating the number of edges
    for(inode = 1; inode <= (*psp)->length; inode++)
        for(jnode = 0; jnode < (*psp)->nexts[inode].length; jnode++)
            if((*psp)->nexts[inode].int_data[jnode] > inode)
                estim++;
    //Allocating the edge array
    *e2n = (int *) calloc(2*estim+2, sizeof(int));

    //Filling the edge array
    for(inode = 1; inode <= (*psp)->length; inode++)
        for(jnode = 0; jnode < (*psp)->nexts[inode].length; jnode++)
            if((*psp)->nexts[inode].int_data[jnode] > inode)
            {
                (*nedge)++;
                //Putting the edge in ascending format
                (*e2n)[2*(*nedge)+0] = inode;
                (*e2n)[2*(*nedge)+1] = (*psp)->nexts[inode].int_data[jnode];
            }
  //return successfully
  return 0;
}

int Gen_Map_Elem__Sur_Elem(int nelem[],
			 int *c2n[],
			 int **nesp,
                         int **esp,
			 struct int_vec **ese)
{
    //declaring local variables
    int indx, gelem, gelem2, surrN, ielem,etype, mnode, inode, node_num;
    int indx1, indx2;
    int *elem_repeatation;
    int allelem = 0; //the total number of all elements in the block
    //calculating the total number of elements in block
    for(indx = Triangle; indx <= Hex; indx++)
        allelem += nelem[indx];
    //allocating the element repeatation array for indexing surrounding elements
    elem_repeatation = (int *)calloc(allelem + 1, sizeof(int));

    //Allocating element-surrounding-ellement linked-list array structure
    *ese = (struct int_vec *)calloc(1,sizeof(struct int_vec));
    (*ese)->length = allelem;
    (*ese)->nexts = (struct int_vec *)calloc(allelem + 1,sizeof(struct int_vec));

    //finding ese for all element in the block
    for(gelem = 1; gelem <= allelem; gelem++)
    {
        //find the local parameters of that element
        Global_Elem_To_Local_Elem(&ielem,nelem,gelem,&etype);
        mnode = ElemTypeToMnode(etype); //find the number of nodes in element
        //Initializing elem_repeatation to zero
        for(gelem2 = 1; gelem2 <= allelem; gelem2++)
            elem_repeatation[gelem2] = 0;

        for(inode = 0; inode < mnode; inode++)
        {
        //find the index of given node (one-based index) in that
        //element so that we can recognize neighbouring elements
               node_num = c2n[etype][mnode*ielem + inode];
               //find esp
               indx1 = (*nesp)[node_num + 0];
               indx2 = (*nesp)[node_num + 1];
               for(indx = indx1; indx < indx2; indx++)
                   elem_repeatation[(*esp)[indx]]++; //mark that the esp[indx]^th
               //element COULD be a neighour!
        }

    //counting the number of surrounding elements
    surrN = 0;
    for(gelem2 = 1; gelem2 <= allelem; gelem2++)
        if((elem_repeatation[gelem2] >= 3) && (gelem != gelem2))
            surrN ++;
    //Now proceed to create the entry of ese for the given element and store ese
    (*ese)->nexts[gelem].length = surrN;
    (*ese)->nexts[gelem].int_data = (int *)calloc(surrN+1,sizeof(int));
    int j = 1;
    for(gelem2 = 1; gelem2 <= allelem;gelem2++)
        if((elem_repeatation[gelem2] >= 3) && (gelem != gelem2))
        {
            (*ese)->nexts[gelem].int_data[j] = gelem2;
            j++;
        }

    }
//clean up
free(elem_repeatation);


return 0; //completed successfully!

}

//The following creates edges surrounding points structure based on
//Gen_Edg_Point() map that is already ctreated before.

// edsp --------->node1-------------->index of edge 1 in e2n belonging to node1
//           |                  |
//           ---->node2         ----->index of edge 2 in e2n belonging to node1
//           |                  |
//           ---->node3         ----->index of edge m in e2n belonging to node1
//           |
//          ...
//           |
//           ---->nodeN

int Create_Edges_Surr_Point(int *nnodes, int *nedge, int **e2n, struct int_vec **edsp)
{
    //Declaring local variables
    int indx, inode, n_ed_belong_inode = 0, k;
    //Allocating memory for edge surrounding point structure
    (*edsp) = (struct int_vec *)calloc(1,sizeof(struct int_vec));
    (*edsp)->length = *nnodes;
    (*edsp)->nexts = (struct int_vec *)calloc(*nnodes+1,sizeof(struct int_vec));

    //Filling up the edge surr point structure
    for(inode = 1; inode <= (*nnodes); inode++)
    {
        //resetting "n_ed_belong_inode" variable
        n_ed_belong_inode = 0;
        //Calculating the number of edge belonging to the inode point
        for(indx = 1; indx<= *nedge; indx++)
            if(((*e2n)[2*indx+0] == inode) || ((*e2n)[2*indx+1] == inode))
                n_ed_belong_inode++;

        //allocating memory for holding the indices of edges belonging to inode
        (*edsp)->nexts[inode].length = n_ed_belong_inode;
        (*edsp)->nexts[inode].int_data = (int *)calloc(n_ed_belong_inode+1 , sizeof(int));

        //filling the array of edges belonging to the inode point with indices
        //of edges (from 1 to nedge)
        //resetting variable "k",
        k = 1;
        for(indx = 1; indx<= *nedge; indx++)
            if(((*e2n)[2*indx+0] == inode) || ((*e2n)[2*indx+1] == inode))
            {
                (*edsp)->nexts[inode].int_data[k] = indx;
                k++;
            }
    }//end of finding the edges surrounding node for each node

    //completed successfully
    return 0;
}

///////////////////////////////////////////////////////////////////
//The following function accepts two integer arrays {arr1[1], ..., arr1[N1]}
//and {arr2[1], ..., arr2[N2]} and returns the number of common integers that
//exist in both arrays. Note that each integer in each array should be repeated
//only ONCE in that array. The final result is put in result
int obtain_common_ints(int *arr1, int *N1, int *arr2, int *N2, int *result)
{
    //introducing local variables
    int ii, jj;
    //saftey resetting operation
    *result = 0;

    for(ii = 1; ii<= (*N1); ii++)
        for(jj=1; jj<= (*N2); jj++)
            if(arr1[ii] == arr2[jj])
            {
                (*result)++;
                break;
            }

    //completed successfully!
    return 0;
}

int Create_Node_Centered_Grid(int *nnodes,int nelem[],
			 int **nesp,
                         int **esp,
			 int *c2n[],
                         const double *xyz,
                         int *nedge,
                         int **e2n,
                         struct double_vec **node_centered, //output
                         struct double_vec **bn_struct)     //output
{
    //Declaring local variables
    int inode = 0, indx1 = 0, indx2 = 0, indx = 0;
    int gelem = 0, ielem = 0, etype = 0, mnode = 0;
    int k,l,p;
    int this_edge;
    int obtain_common_ints_result = 0, N1, N2;
    //free indices
    int ii;
    double x1,y1,z1; //coordinate of current point in the main loop
    double x2,y2,z2; //centroid of the edge
    double x3,y3,z3; //centroid of the face
    double x4,y4,z4; //centroid of the element
    double xp,yp,zp; //centroid of the volume element surrounding boundary element.
    int p_elem_index = 0; //The index of the volume element surrounding boundary element.
    int p_elem_etype = 0, p_elem_ielem = 0; //The local properties of the volume element surrounding boundary element.
    int p_elem_mnode = 0; //The number of nodes available in it.
    int *p_elem_nodes = (int *)calloc(1,sizeof(int));  //The array containing the nodes of the p-element (just initialized)


    double DeltaV = 0.0;  //Volume contribution of each tet formed by joining centroids
    double DeltaSx = 0.0; //Area_x contribution of each tet formed by joining centroids
    double DeltaSy = 0.0; //Area_y contribution of each tet formed by joining centroids
    double DeltaSz = 0.0; //Area_z contribution of each tet formed by joining centroids
    //the area vector at the bottom of the contributor tet is S = Ax e1+ Ay e2 + Az e3;
    double sign_dot = 0.0; //the result of dot product of edge direction to the contributor area normal
    int ghost_edge_index = 0; //The index of ghost edge.

    int local_node = 0;
    struct int_vec *edsp; //A structure for holding edges surr points
    int *lst_elem_nodes = (int *)calloc(1,sizeof(int));
    int *lst_face_nodes = (int *)calloc(1,sizeof(int)); //local arrays(the size is variable,
    //so we don't REALLY allocate them here.
    int *lst_edge_nodes = (int *)calloc(3, sizeof(int));

    //The following structures are must-have maps
    struct int_vec *surr_face_el_node; //to find faces of an element surrounding a point in that element
    struct int_vec *pnts_face_el; // to find the local index of nodes on each surface tags
    struct int_vec *ese; //for boundary elements to find the corresponding interior
    //VOLUME element so that we can understand what the normal direction of boundary face is.
    short *is_it_boundary_node = (short *)calloc((*nnodes +1),sizeof (short)); //a vector
    //of shorts with length *nnodes+1 which each element is filled with
    //one: showing that the corresponding node is boundary node and zero: meaning that it's
    //not a boundary node.

    //Creating face-tag structure which accept element type and a node in element
    //and returns the surface tags
    Create_Surr_Face_El_Node(&surr_face_el_node);

    //Creating local index of points which are in a face-tag of an element
    Create_Points_In_Face_El(&pnts_face_el);

    //Find element surrounding elements for all elements inside the domain.
    Gen_Map_Elem__Sur_Elem(nelem,c2n,nesp,esp,&ese);

    //Detecting which node is a boundary node.
    //initializing
    for(inode = 1; inode <= (*nnodes); inode++)
    {
        is_it_boundary_node[inode] = 0; //all nodes are initially interrior
    }

    for(inode = 1; inode <= (*nnodes); inode++)
    {
        indx1 = (*nesp)[inode + 0];
        indx2 = (*nesp)[inode + 1];
        //Here is the list of surr elements with indices in (*esp)[indx]
        for(indx = indx1; indx < indx2; indx++)
        {
            //Extracting global points in element
            gelem = (*esp)[indx]; //get the golobal index of the neighbur elem
            //find the local parameters of that element
            Global_Elem_To_Local_Elem(&ielem,nelem,gelem,&etype);
            if(etype <= Quad)
            {
              is_it_boundary_node[inode] = 1;
              break;
            }

        }

    }//finished specifying boundary nodes

    //Since we already know the number of nodes, the number of edges,
    //boundary nodes and their ghost edges and also
    //the structure that volume and faces must be saved in, we can create
    //the "struct int_vec **node_centered" right here.
    //The format of "node_centered" is given below,
    //node_centered->[point]->[edge surr point index]->{vol(i),S(i),nx(i),ny(i),nz(i)}
    //otherwise node_centered->[point]->[ghost edge]->{-1.0,S(i)(boundary),nx(i)(boundary),ny(i)(boundary),nz(i)(boundary)}
    //BUILDING THE ROOT OF TREE (NODES INFO)
    (*node_centered) = (struct double_vec*)calloc(1,sizeof (struct double_vec));
    (*node_centered)->length = (*nnodes);
    (*node_centered)->nexts = (struct double_vec *)calloc((*nnodes +1),sizeof (struct double_vec));
    //BUILDING THE JUNCTIONS (EDGE INFO)
    //first finding edges surr points
    Create_Edges_Surr_Point(nnodes, nedge, e2n, &edsp);
    //Proceeding to make a branch for each edge surr point for all points
    for(inode = 1; inode <= (*nnodes); inode++)
    {
        if(!is_it_boundary_node[inode]) //not boundary node
        {
            (*node_centered)->nexts[inode].length = edsp->nexts[inode].length;
            (*node_centered)->nexts[inode].nexts =(struct double_vec *)calloc((edsp->nexts[inode].length+1),sizeof (struct double_vec));
            //Proceeding to make volume-area-normal structure for each eadges
            for(indx = 1; indx <= (*node_centered)->nexts[inode].length; indx++)
            {
                (*node_centered)->nexts[inode].nexts[indx].length = 5; //vol(i),S(i),nx(i),ny(i),nz(i)
                (*node_centered)->nexts[inode].nexts[indx].db_data = (double *)calloc(6 , sizeof(double));
                //initializing the entry to zero.
                for(indx1 = 1; indx1 <=5; indx1++)
                    (*node_centered)->nexts[inode].nexts[indx].db_data[indx1] = 0.0;
            }
        }else //it's a boundary node up there
        {
            (*node_centered)->nexts[inode].length = edsp->nexts[inode].length + 1; //one additional for ghost node
            (*node_centered)->nexts[inode].nexts =(struct double_vec *)calloc((edsp->nexts[inode].length+2),sizeof (struct double_vec));
            //Proceeding to make volume-area-normal structure for each eadges
            //doing the action for node 1 ... length which are interior-pointed edges
            for(indx = 1; indx < (*node_centered)->nexts[inode].length; indx++)
            {

                (*node_centered)->nexts[inode].nexts[indx].length = 5; //vol(i),S(i),nx(i),ny(i),nz(i)
                (*node_centered)->nexts[inode].nexts[indx].db_data = (double *)calloc(6 , sizeof(double));
                //initializing the entry to zero.
                for(indx1 = 1; indx1 <=5; indx1++)
                    (*node_centered)->nexts[inode].nexts[indx].db_data[indx1] = 0.0;

            }
            //completing boundary node edge specification for ghost edge (the last one)
            (*node_centered)->nexts[inode].nexts[indx].length = 5; //vol(i),S(i),nx(i),ny(i),nz(i)
            (*node_centered)->nexts[inode].nexts[indx].db_data = (double *)calloc(6 , sizeof(double));
            //initializing the entry to zero.
            for(indx1 = 2; indx1 <=5; indx1++)
                (*node_centered)->nexts[inode].nexts[indx].db_data[indx1] = 0.0;
            //Initializing the VOLUME to -1 to show that is a ghost edge.
                (*node_centered)->nexts[inode].nexts[indx].db_data[1] = -1.;
        }

    }
    //Proceeding to create boundary structure containing all areas and tags for each boundary edge per boundary element.
    int bn_elem_size = nelem[Triangle] + nelem[Quad];
    (*bn_struct) = (struct double_vec*)calloc(1,sizeof (struct double_vec));
    (*bn_struct)->length = bn_elem_size;
    (*bn_struct)->nexts = (struct double_vec *)calloc((bn_elem_size +1),sizeof (struct double_vec));
    //Proceeding to make a branch for each boundary element
    int ibelem = 0; //global index of boundary element
    for(ibelem = 1; ibelem <= nelem[Triangle]; ibelem++) //assigning branchs for triangles
    {
            (*bn_struct)->nexts[ibelem].length = 3; //there are three nodes in a triangle and hence three ghost edges
            (*bn_struct)->nexts[ibelem].nexts =(struct double_vec *)calloc(4,sizeof (struct double_vec));
            //Proceeding to assign edges info for each nodes in triangles
            for(indx = 1; indx <= (*bn_struct)->nexts[ibelem].length; indx++)
            {
                (*bn_struct)->nexts[ibelem].nexts[indx].length = 3; //Sx, Sy, Sz
                (*bn_struct)->nexts[ibelem].nexts[indx].db_data = (double *)calloc(4 , sizeof(double));
                //initializing the entry to zero.
                for(indx1 = 1; indx1 <=3; indx1++)
                    (*bn_struct)->nexts[ibelem].nexts[indx].db_data[indx1] = 0.0;
            }
    }//end of edge allocation for triangles

    for(ibelem = (nelem[Triangle] +1); ibelem <= bn_elem_size; ibelem++) //assigning branchs for quads
    {
            (*bn_struct)->nexts[ibelem].length = 4; //there are four in a quad and hence four ghost edges
            (*bn_struct)->nexts[ibelem].nexts =(struct double_vec *)calloc(5,sizeof (struct double_vec));
            //Proceeding to assign edges info for each nodes in quad
            for(indx = 1; indx <= (*bn_struct)->nexts[ibelem].length; indx++)
            {
                (*bn_struct)->nexts[ibelem].nexts[indx].length = 3; //Sx, Sy, Sz
                (*bn_struct)->nexts[ibelem].nexts[indx].db_data = (double *)calloc(4 , sizeof(double));
                //initializing the entry to zero.
                for(indx1 = 1; indx1 <=3; indx1++)
                    (*bn_struct)->nexts[ibelem].nexts[indx].db_data[indx1] = 0.0;
            }
    }//end of edge allocation for quad


    // For all nodes in the grid create the node centered map.
    for(inode = 1; inode <= (*nnodes); inode++)
    {

        //calculating the coordinate of current point
        x1 = xyz[3*inode + 0]; //getting the x-coor of this_node in element
        y1 = xyz[3*inode + 1]; //getting the y-coor of this_node in element
        z1 = xyz[3*inode + 2]; //getting the z-coor of this_node in element

        //we already have the edges surr point inode so we don't need to calc it
        //loop over elements surrounding point "inode",
        indx1 = (*nesp)[inode + 0];
        indx2 = (*nesp)[inode + 1];
        //Here is the list of surr elements with indices in (*esp)[indx]
        for(indx = indx1; indx < indx2; indx++)
        {
            //Element-Specific reset
            //Resetting centroid arrays for the next operations
            x4 = 0.0;
            y4 = 0.0;
            z4 = 0.0;
            free(lst_elem_nodes); //Freeing the allocated memory for
            //the previous elements, if any

            //calculate the centroid of the element.
            ////////////////////////////////////////
            //Extracting global points in element
            gelem = (*esp)[indx]; //get the golobal index of the neighbur elem
            //find the local parameters of that element
            Global_Elem_To_Local_Elem(&ielem,nelem,gelem,&etype);
            mnode = ElemTypeToMnode(etype); //find the number of nodes in element
            //once the number of nodes in element (mnode) is known, we proceed
            //to allocate memory for holding the global index of nodes in the
            //allocated array.
            lst_elem_nodes = (int *)calloc(mnode + 1, sizeof(int));


            //collecting element nodes in the "lst_elem_nodes" array.
            for(ii = 0; ii < mnode; ii++)
                lst_elem_nodes[ii+1] = c2n[etype][mnode*ielem + ii];

            //Proceeding to SUM the X, Y and Z components of element points coordinates
            //which we after use to find the average which is the centroid.
            for(ii = 1; ii <= mnode; ii++)
            {

                x4 += xyz[3*lst_elem_nodes[ii] + 0]; //getting the x-coor of this_node in element
                y4 += xyz[3*lst_elem_nodes[ii] + 1]; //getting the y-coor of this_node in element
                z4 += xyz[3*lst_elem_nodes[ii] + 2]; //getting the z-coor of this_node in element
            }
            //finding the averages wich is centroid of the element
            x4 /= ((double)mnode);
            y4 /= ((double)mnode);
            z4 /= ((double)mnode);

            ////////////////////////////////////////////////////////
            //p-element related resettings /////////////////////////
            ////////////////////////////////////////////////////////
            if(etype <= Quad) //This sufficies to prove that current node,edge and element and face are on boundary
            {
                xp = 0.0;
                yp = 0.0;
                zp = 0.0;
                free(p_elem_nodes); //if any
                //Do boundary operation
                //finding the index of boundary ghost edge
                ghost_edge_index = (*node_centered)->nexts[inode].length; // The ghost edge index
                //Calculating the area vector
                //1- finding the centeroid of volume element surr boundary element (p-element) to construct
                //the vector showing the direction twoard the interior.
                //finding the interior volume element.
                p_elem_index = ese->nexts[gelem].int_data[1];
                //find the local parameters of that element
                Global_Elem_To_Local_Elem(&p_elem_ielem,nelem,p_elem_index,&p_elem_etype);
                p_elem_mnode = ElemTypeToMnode(p_elem_etype); //find the number of nodes in element
                //once the number of nodes in p-element (p_elem_mnode) is known, we proceed
                //to allocate memory for holding the global index of nodes in the
                //allocated array.
                p_elem_nodes = (int *)calloc(p_elem_mnode + 1, sizeof(int));
                //collecting element nodes in the "p_elem_nodes" array.
                for(ii = 0; ii < p_elem_mnode; ii++)
                    p_elem_nodes[ii+1] = c2n[p_elem_etype][p_elem_mnode*p_elem_ielem + ii];
                //Proceeding to SUM the X, Y and Z components of p-element points coordinates
                //which we after use to find the average which is the centroid of the p-element.
                for(ii = 1; ii <= p_elem_mnode; ii++)
                {

                    xp += xyz[3*p_elem_nodes[ii] + 0]; //getting the x-coor of the p-element
                    yp += xyz[3*p_elem_nodes[ii] + 1]; //getting the y-coor of the p-element
                    zp += xyz[3*p_elem_nodes[ii] + 2]; //getting the z-coor of the p-element
                }
                //finding the averages which is centroid of the element
                xp /= ((double)p_elem_mnode);
                yp /= ((double)p_elem_mnode);
                zp /= ((double)p_elem_mnode);

            //Now if we connect the center of the BOUNDARY element (x4,y4,z4) to the center of p-element
            //then we have the vector in the direction toward the interior
            //our final boundary vector should be in opposite direction.
            }
            ////////////////////////////////////////////////////////

            //Now, we need two remaining averages; center of edge and center of face
            //first, filter the edges surr point to find edges(k) which are in that elemen
            for(k=1; k<= edsp->nexts[inode].length; k++)
            {
                //Extracting global index of each edge
                this_edge = edsp->nexts[inode].int_data[k];
                //Finding two points of each edge.
                lst_edge_nodes[1] = (*e2n)[2*this_edge+ 0];
                lst_edge_nodes[2] = (*e2n)[2*this_edge+ 1];
                //Proceeding to check if this edge IS in the element(indx)
                N1 = 2; //number of ints in edge array
                N2 = mnode; //number of ints in element array
                obtain_common_ints(lst_edge_nodes, &N1, lst_elem_nodes, &N2, &obtain_common_ints_result);
                if(obtain_common_ints_result == 2) //that means two-point edge IS in the element
                {
                    //loop over edges(k)
                    //edge- specific resseting commands
                    x2 = 0.0;
                    y2 = 0.0;
                    z2 = 0.0;

                    //Proceeding to SUM the X, Y and Z components of edge points coordinates
                    //which we after use to find the average which is the centroid of the edge
                    for(ii = 1; ii <= 2; ii++)
                    {

                        x2 += xyz[3*lst_edge_nodes[ii] + 0]; //getting the x-coor of this_node in element
                        y2 += xyz[3*lst_edge_nodes[ii] + 1]; //getting the y-coor of this_node in element
                        z2 += xyz[3*lst_edge_nodes[ii] + 2]; //getting the z-coor of this_node in element
                    }
                    //finding the averages
                    x2 /= 2.;
                    y2 /= 2.;
                    z2 /= 2.;

                    //Now, proceeding to find the last thing: face centroids
                    //first we use local info such as element type, local node
                    //to find the number and tag of faces surrounding node inode.
                    //ones the tags are obtained we use another map to find the
                    //local index of nodes inside each face and then find the
                    //corresponding global node by using transformation.
                    //after obtaining all nodes (global index) in each face
                    //we simply check if the current edge is in that face
                    //or not. If yes, we proceed to calculate area and volumes
                    //else we change the face until we find the face that contain
                    //the edge
                    //looping over face(l)
                    //since all operation in face are dependent to LOCAL index of
                    //node, we should convert inode to the corresponding position
                    //in the element using ONE-BASED convention
                    //
                    for(ii = 0; ii < mnode; ii++)
                        if(inode == c2n[etype][mnode*ielem + ii])
                            local_node = ii+1; //one-based


                    for(l=1; l < surr_face_el_node->nexts[etype].nexts[local_node].length; l++)
                    {
                        //Face-specific resetting and cleanup
                        x3 = 0.0;
                        y3 = 0.0;
                        z3 = 0.0;
                        free(lst_face_nodes); //if any.

                        //allocating "lst_face_nodes" based on the number of nodes in te face(l)
                        lst_face_nodes = (int *)calloc(pnts_face_el->nexts[etype].nexts[surr_face_el_node->nexts[etype].nexts[local_node].int_data[l]].length,sizeof(int));

                        //extracting the local nodes of face(l)
                        for(p=1; p< pnts_face_el->nexts[etype].nexts[surr_face_el_node->nexts[etype].nexts[local_node].int_data[l]].length; p++)
                            lst_face_nodes[p] = pnts_face_el->nexts[etype].nexts[surr_face_el_node->nexts[etype].nexts[local_node].int_data[l]].int_data[p];
                        //converting local node index in each face to global node index
                        for(p=1; p< pnts_face_el->nexts[etype].nexts[surr_face_el_node->nexts[etype].nexts[local_node].int_data[l]].length; p++)
                            lst_face_nodes[p] = c2n[etype][mnode*ielem + lst_face_nodes[p]-1];
                        //check if the edge(k) is in the face(l)
                        N2 = pnts_face_el->nexts[etype].nexts[surr_face_el_node->nexts[etype].nexts[local_node].int_data[l]].length-1;
                        obtain_common_ints(lst_edge_nodes, &N1, lst_face_nodes, &N2, &obtain_common_ints_result);
                        if(obtain_common_ints_result == 2) //face(l) includes edge(k)
                        {
                            //Proceeding to SUM the X, Y and Z components of face points coordinates
                            //which we after use to find the average which is the centroid of the face
                            for(ii = 1; ii <= N2; ii++)
                            {

                                x3 += xyz[3*lst_face_nodes[ii] + 0]; //getting the x-coor of this_node in element
                                y3 += xyz[3*lst_face_nodes[ii] + 1]; //getting the y-coor of this_node in element
                                z3 += xyz[3*lst_face_nodes[ii] + 2]; //getting the z-coor of this_node in element
                            }
                            //Finding the averages
                            x3 /= ((double)N2);
                            y3 /= ((double)N2);
                            z3 /= ((double)N2);

                            //Now, (I'm so tired of this complex algorithm)
                            //Now, we have three centroids plus point coordinate.
                            //That means we should start calculating volume and
                            //face contibution and then put them in the right
                            //position in the node_centered structure. Note that if
                            //the etype is tri or quad that means that the element
                            //is a boundary element by definition convention of
                            //ugrid input format. So we expect that it doesn't have
                            //any volume contrubution and also surface contribution
                            //to the INTERIOR ELEMENTS. So we should filter boundary elements
                            //and in the case of boundary element CONTRIBUTOR, we should calculate
                            //the area normal vector that is NORMAL to boundary element, and just
                            //put it in the ghost node index of that boundary node for Ax, Ay, Az.
                            //A post-processing routine WILL be written that extract the data accordingly and makes
                            //the final structure that we want when developing solver.
                            //printf("\t%d\n\r",inode);
                            if(etype <= Quad) //This sufficies to prove that current node,edge and element and face are on boundary
                            {

                                //Proceeding to calculate the area of the boundary contributor.
                                //Calculating the area vector contribution for boundary node
                                DeltaSx = 1./2.*( (y1-y2)*(z3-z2) - (z1-z2)*(y3-y2) );
                                DeltaSy = 1./2.*( (z1-z2)*(x3-x2) - (x1-x2)*(z3-z2) );
                                DeltaSz = 1./2.*( (x1-x2)*(y3-y2) - (y1-y2)*(x3-x2) );

                                //Cheking if the normal vector is in the opposite direction of the p-element vector
                                //then applying modification of the sign.
                                //lemma: if we dot-product vector V1= (x3,y3,z3)->(xp,yp,zp)
                                //and S = DeltaSx i + DeltaSy j + DeltaSz k and the result is -NEGATIVE-
                                //then the normal is in the right direction and sign is OK
                                //otherwise we should multiply S with -1 to correct the direction
                                sign_dot = (xp - x3)*DeltaSx + (yp - y3)*DeltaSy + (zp - z3)*DeltaSz;

                                if( sign_dot > 0.) //a sign correction needed
                                {
                                    DeltaSx *= -1.;
                                    DeltaSy *= -1.;
                                    DeltaSz *= -1.;
                                }
                                else if ( (sign_dot == 0.0) && ( (DeltaSx != 0.0) || (DeltaSy != 0.0) || (DeltaSz != 0.0) ) )    //That's wierd , exit
                                {
                                    printf("fatal error!\n\r bad boundary element or boundary element facing an interior element with zero volume V = 0. Check your grid input.\n\r");
                                    exit(0);
                                }

                                //Adding the contributation to the ghost edge for that boundary node
                                //Note that "ghost_edge_index" is already calculated in the element-loop for etype<= Quad
                                //(*node_centered)->nexts[inode].nexts[ghost_edge_index].db_data[3] += DeltaSx;
                                //(*node_centered)->nexts[inode].nexts[ghost_edge_index].db_data[4] += DeltaSy;
                                //(*node_centered)->nexts[inode].nexts[ghost_edge_index].db_data[5] += DeltaSz;
                                (*bn_struct)->nexts[gelem].nexts[local_node].db_data[1] += DeltaSx;
                                (*bn_struct)->nexts[gelem].nexts[local_node].db_data[2] += DeltaSy;
                                (*bn_struct)->nexts[gelem].nexts[local_node].db_data[3] += DeltaSz;


                            }
                            else //the node is in the interrior
                            {
                                //Calculating the volume contribution
                                DeltaV = 1.0/6.0*fabs( (x1-x4)* ((y2-y4)*(z3-z4)-(z2-z4)*(y3-y4)) + (y1-y4) *((z2-z4)*(x3-x4) - (x2-x4)*(z3-z4)) + (z1-z4) * ((x2-x4)*(y3-y4) - (y2-y4)*(x3-x4)));
                                //Adding the contributation to the right place in node-centered database
                                (*node_centered)->nexts[inode].nexts[k].db_data[1] += DeltaV;

                                //Calculating the area vector contribution
                                DeltaSx = 1./2.*( (y4-y3)*(z2-z3) - (z4-z3)*(y2-y3) );
                                DeltaSy = 1./2.*( (z4-z3)*(x2-x3) - (x4-x3)*(z2-z3) );
                                DeltaSz = 1./2.*( (x4-x3)*(y2-y3) - (y4-y3)*(x2-x3) );

                                //cheking if the normal vector sign is righ-hand or not
                                //then applying modification of the sign.
                                //lemma: if we dot-product vector V1= (x1,y1,z1)->(x2,y2,z2)
                                //and S = DeltaSx i + DeltaSy j + DeltaSz k and the result is positive
                                //then the normal is in the right direction and sign is OK
                                //otherwise we should multiply S with -1 to correct the direction
                                sign_dot = (x2 - x1)*DeltaSx + (y2 - y1)*DeltaSy + (z2 - z1)*DeltaSz;

                                if( sign_dot < 0.) //a sign correction needed
                                {
                                    DeltaSx *= -1.;
                                    DeltaSy *= -1.;
                                    DeltaSz *= -1.;
                                }
                                else if ( (sign_dot == 0.0) && ( (DeltaSx != 0.0) || (DeltaSy != 0.0) || (DeltaSz != 0.0) ) )    //That's wierd , exit
                                {
                                    printf("fatal error.\n\r The direction of edge is prependicular to contributor face normal. Check your grid input.\n\r");
                                    exit(0);
                                }

                                //Adding the contributation to the right place in node-centered database
                                (*node_centered)->nexts[inode].nexts[k].db_data[3] += DeltaSx;
                                (*node_centered)->nexts[inode].nexts[k].db_data[4] += DeltaSy;
                                (*node_centered)->nexts[inode].nexts[k].db_data[5] += DeltaSz;
                            }


                            //printf("for node=%d, gelem = %d, etype=%d, edge=%d, face=%d,  DeltaV = %f\n\r",inode,gelem, etype, edsp->nexts[inode].int_data[k], l, DeltaV);
                            //printf("%f\n\r",DeltaV);



                        }//end of calculations for face+edge+element

                    }//end of face(l)
                }//end of edge in element
            }//end of edge(k)
        }//end of elements surr point (esp[indx]) loop
    }// end of main loop for all points in grid ("inode")

    //Proceeding to cleanup locally generated structures
    free(edsp);
    free(lst_elem_nodes);
    free(lst_face_nodes);
    free(lst_edge_nodes);
    free(surr_face_el_node);
    free(pnts_face_el);

    //completed successfully!
    return 0;
}

int Node_Based_To_Array_Convertor(int *nnodes, int *nedge,                              // input variables
                         int **e2n, int nelem[], int *c2n[],                            // input variables
                         struct int_vec **edsp, int*factag[],                           // input variables
                         struct double_vec **node_centered,                             // input variables
                         struct double_vec **bn_struct,                                 // input variables
                         int **bn_edge_to_node_map, double **Vol,                       // output variables
                         double **Sx_i, double **Sy_i, double **Sz_i, double **S_i,     // output variables
                         double **Sx_b, double **Sy_b, double **Sz_b, double **S_b,     // output variables
                         int *bn_size, int **tag_b)                                     // output variables
{
    //allocating local variables
    int inode = 0, iedge = 0, this_edge = 0, k = 0;
    int *is_this_edge_processed; //a [1,nedge] vector containing global edge indices that are zeros for edges
    //that their area has not been computed sofar and 1 when it is calculated.
    double the_edge_sign = 1.; //the sign of area vector formed by each interior edges

    //Initializing output variables
    (*bn_size) = 3*nelem[Triangle] + 4*nelem[Quad];

    (*Vol) = (double *)calloc((*nnodes +1),sizeof (double)); //The volume of each node-centered control volume.

    //initializing node-based values
    for(inode = 1; inode <= (*nnodes); inode++)
        (*Vol)[inode] = 0.;


    (*Sx_i) = (double *)calloc((*nedge +1),sizeof (double));
    //The x-component of the areavector for the interior edges Sx_i[1..nedges]
    (*Sy_i) = (double *)calloc((*nedge +1),sizeof (double));
    //The y-component of the areavector for the interior edges Sy_i[1..nedges]
    (*Sz_i) = (double *)calloc((*nedge +1),sizeof (double));
    //The z-component of the areavector for the interior edges Sz_i[1..nedges]
    (*S_i) = (double *)calloc((*nedge +1),sizeof (double));
    //The signed magnitude of the areavector for the interior edges. S_i[1..nedges].
    //It is positive when the edge is going from lower node index to higher node index.
    is_this_edge_processed = (int *)calloc((*nedge +1),sizeof (int));
    //initializing edge-based values
    for(iedge = 1; iedge <= (*nedge); iedge++)
    {
        (*Sx_i)[iedge] = 0.;
        (*Sy_i)[iedge] = 0.;
        (*Sz_i)[iedge] = 0.;
        (*S_i)[iedge] = 0.;
        is_this_edge_processed[iedge] = 0;
    }


    //Initializing all fictious-edge-based arrays based on the index of boundary point
    (*bn_edge_to_node_map) = (int *)calloc((*bn_size +1),sizeof (int));
    (*tag_b) = (int *)calloc((*bn_size +1),sizeof (int));
    (*Sx_b) = (double *)calloc((*bn_size +1),sizeof (double)); //The x-component of the fictious areavector for the boundary nodes Sx_b[1..nbnodes]
    (*Sy_b) = (double *)calloc((*bn_size +1),sizeof (double)); //The y-component of the fictious areavector for the boundary nodes Sy_b[1..nbnodes]
    (*Sz_b) = (double *)calloc((*bn_size +1),sizeof (double)); //The z-component of the fictious areavector for the boundary nodes Sz_b[1..nbnodes]
    (*S_b)  = (double *)calloc((*bn_size +1),sizeof (double));  //The magnitude of the fictious areavector for the boundary nodes S_b[1..nbnodes]
    //Note that we don't specify the sign of S_b. instead we let the area vector normal to specify the
    //sign of this vector without changing it.

    //Storing boundary edge info
    int bgelem = 0; //global index of boundary element
    iedge = 0; //resetting iedge
    int ielem = 0; //local index of element
    int etype = 0; //element type
    int mnode = 0; //number of nodes in element.
    //proceeding to store data for triangle boundary elements.
    for(bgelem = 1; bgelem <= nelem[Triangle]; bgelem++)
        for(inode = 1; inode <=3 ; inode++)
        {
            iedge++;
            (*Sx_b)[iedge] = (*bn_struct)->nexts[bgelem].nexts[inode].db_data[1];
            (*Sy_b)[iedge] = (*bn_struct)->nexts[bgelem].nexts[inode].db_data[2];
            (*Sz_b)[iedge] = (*bn_struct)->nexts[bgelem].nexts[inode].db_data[3];
            //storing boundary area (fictious) vector magnitude
            (*S_b)[iedge] = sqrt(pow((*Sx_b)[iedge], 2.) + pow((*Sy_b)[iedge] , 2.0) + pow((*Sz_b)[iedge], 2.));
            if((*S_b)[iedge] == 0.) //putting a warning -divided by zero prevention
            {
                printf("Warning, the magnitude |S| of the area vector at a boundary node on a ghost edge is zero, division by zero danger.\n\r");
                exit(0);
            }
            //normalizing the area vector to area normal without changing the sign.
            (*Sx_b)[iedge] /= (*S_b)[iedge];
            (*Sy_b)[iedge] /= (*S_b)[iedge];
            (*Sz_b)[iedge] /= (*S_b)[iedge];
            //storing boundary fictious edge to real node map
            Global_Elem_To_Local_Elem(&ielem,nelem,bgelem,&etype); //global to local transformation
            mnode = ElemTypeToMnode(etype);
            (*bn_edge_to_node_map)[iedge] = c2n[etype][mnode*ielem + inode -1]; //storing the global number of real node which has that ghost edge
            //storin boundary tag
            (*tag_b)[iedge] = factag[etype][ielem];
        }
    //proceeding to store data for quadlateral boundary elements.
    for(bgelem = (nelem[Triangle]+1); bgelem <= (nelem[Triangle] + nelem[Quad]); bgelem++)
        for(inode = 1; inode <=4 ; inode++)
        {
            iedge++;
            (*Sx_b)[iedge] = (*bn_struct)->nexts[bgelem].nexts[inode].db_data[1];
            (*Sy_b)[iedge] = (*bn_struct)->nexts[bgelem].nexts[inode].db_data[2];
            (*Sz_b)[iedge] = (*bn_struct)->nexts[bgelem].nexts[inode].db_data[3];
            //storing boundary area (fictious) vector magnitude
            (*S_b)[iedge] = sqrt(pow((*Sx_b)[iedge], 2.) + pow((*Sy_b)[iedge] , 2.0) + pow((*Sz_b)[iedge], 2.));
            if((*S_b)[iedge] == 0.) //putting a warning -divided by zero prevention
            {
                printf("Warning, the magnitude |S| of the area vector at a boundary node on a ghost edge is zero, division by zero danger.\n\r");
                exit(0);
            }
            //normalizing the area vector to area normal without changing the sign.
            (*Sx_b)[iedge] /= (*S_b)[iedge];
            (*Sy_b)[iedge] /= (*S_b)[iedge];
            (*Sz_b)[iedge] /= (*S_b)[iedge];
            //storing boundary fictious edge to real node map
            Global_Elem_To_Local_Elem(&ielem,nelem,bgelem,&etype); //global to local transformation
            mnode = ElemTypeToMnode(etype);
            (*bn_edge_to_node_map)[iedge] = c2n[etype][mnode*ielem + inode -1]; //storing the global number of real node which has that ghost edge
            //storin boundary tag
            (*tag_b)[iedge] = factag[etype][ielem];
        }

    //Done boundary conversion to array type
    //Proceeding to convert the interior edges area and volume contribution
    for(inode = 1; inode <= (*nnodes); inode++) //search for all nodes in the domain. and for all edges except the fictious boundary edges (V=-1.0)
        for(k=1; k<= (*edsp)->nexts[inode].length; k++)
        {
            //adding volume data
            (*Vol)[inode] += (*node_centered)->nexts[inode].nexts[k].db_data[1];
            //obtainin the global edge information
            //Extracting global index of each edge
            this_edge = (*edsp)->nexts[inode].int_data[k];
            if(!is_this_edge_processed[this_edge]) //this edge is not processed
            {
                //Storing interior area vector component
                (*Sx_i)[this_edge] = (*node_centered)->nexts[inode].nexts[k].db_data[3];
                (*Sy_i)[this_edge] = (*node_centered)->nexts[inode].nexts[k].db_data[4];
                (*Sz_i)[this_edge] = (*node_centered)->nexts[inode].nexts[k].db_data[5];
                //computing the area magnitude at interior points
                (*S_i)[this_edge] = sqrt( pow((*Sx_i)[this_edge], 2.) + pow((*Sy_i)[this_edge], 2.) + pow((*Sz_i)[this_edge], 2.) );

                //calculating the sign of area vector built on the edge "this-edge"
                if( (*e2n)[2*this_edge+ 0] < (*e2n)[2*this_edge+ 1] )
                    the_edge_sign = 1.;
                else
                    the_edge_sign = -1.;


                //putting a warning message
                if((*S_i)[this_edge] == 0.) //putting a warning -divided by zero prevention
                {
                    printf("Warning, the magnitude |S| of the area vector at a interior edge is zero, division by zero danger.\n\r");
                    exit(0);
                }

                (*Sx_i)[this_edge] /= (the_edge_sign * (*S_i)[this_edge]);
                (*Sy_i)[this_edge] /= (the_edge_sign * (*S_i)[this_edge]);
                (*Sz_i)[this_edge] /= (the_edge_sign * (*S_i)[this_edge]);

                //Mark the edge as a processed one so in the next node it is not repeated.
                is_this_edge_processed[this_edge] = 1;
            }
        }

    //clean-up
    free(is_this_edge_processed);

    //completed successfully!
    return 0;
}

