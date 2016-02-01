#include "generic.h"

///////////////////////////////////////////////////////////////////////////////
///MY CODE AND WRITING
///THE FOLLOWING DECLERATIONS ARE NOT IN THE ORIGINAL CODE BUT ADDED BY ME
///////////////////////////////////////////////////////////////////////////////

int Create_Element_Node_Conn_Table(struct int_vec **elem_node_conn)
{
    ///////////////////////////////////////////////////////////////////////////////
//TEST SEGMENT FOR EVALUATING FUNCTIONS ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//Defining element connected points map
//root-------->[pointer to elements]---
//      -Number of possible elements  -
//                                    -
//                                    -
//                                    --->[given node 1]-----
//                                    -N1_num of connecs to1-
//                                    -                     -
//                                    -                     --->point1(gelem)
//                                    -                     --->point2(gelem)
//                                    -                          ...
//                                    -                     --->pointN1
//                                    -
//                                    -
//                                    --->[given node 2]-----
//                                    -N2_num of connecs to2-
//                                    -                     -
//                                    -                     --->point1(gelem)
//                                    -                     --->point2(gelem)
//                                    -                          ...
//                                    -                     --->pointN2
//                                    -
//                                    -
//                                    so on .............
//Note: we can still use our linked-list array structure for this problem so,


  (*elem_node_conn) = (struct int_vec *)calloc(1,sizeof(struct int_vec));
  //input: element type, the index of "given node" in that element
  (*elem_node_conn)->length = 6+1; //tri,quad,tet,pyramid,prism,hex
  //allocating memory for different element types
  (*elem_node_conn)->nexts = (struct int_vec *)calloc((*elem_node_conn)->length,
                                                sizeof(struct int_vec));
  //////defining maps for triangle//////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  (*elem_node_conn)->nexts[1].length = 3+1; //there are three possibilites for
  //input points in triangles
  (*elem_node_conn)->nexts[1].nexts =
      (struct int_vec *)calloc((*elem_node_conn)->nexts[1].length,
      sizeof(struct int_vec));
  (*elem_node_conn)->nexts[1].nexts[1].length = 2 +1; //there are other 2 points
          //near given points 1
  (*elem_node_conn)->nexts[1].nexts[1].int_data = (int *)calloc((2+1),
      sizeof(int));
  (*elem_node_conn)->nexts[1].nexts[1].int_data[1] = 2;
  (*elem_node_conn)->nexts[1].nexts[1].int_data[2] = 3;

  //for given point 2 in triangle
  (*elem_node_conn)->nexts[1].nexts[2].length = 2 +1; //there are other 2 points
          //near given points 2
  (*elem_node_conn)->nexts[1].nexts[2].int_data = (int *)calloc((2+1),
      sizeof(int));
  (*elem_node_conn)->nexts[1].nexts[2].int_data[1] = 1;
  (*elem_node_conn)->nexts[1].nexts[2].int_data[2] = 3;

  //for given point 3 in triangle
  (*elem_node_conn)->nexts[1].nexts[3].length = 2 +1; //there are other 2 points
          //near given points 3
  (*elem_node_conn)->nexts[1].nexts[3].int_data = (int *)calloc((2+1),
      sizeof(int));
  (*elem_node_conn)->nexts[1].nexts[3].int_data[1] = 1;
  (*elem_node_conn)->nexts[1].nexts[3].int_data[2] = 2;

  //////defining maps for quad//////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  (*elem_node_conn)->nexts[2].length = 4+1; //there are four possibilites for
  //input points in quad
  (*elem_node_conn)->nexts[2].nexts =
      (struct int_vec *)calloc((*elem_node_conn)->nexts[2].length,
      sizeof(struct int_vec));
  (*elem_node_conn)->nexts[2].nexts[1].length = 2 +1; //there are 2 other points
          //near given points 1
  (*elem_node_conn)->nexts[2].nexts[1].int_data = (int *)calloc((2+1),
      sizeof(int));
  (*elem_node_conn)->nexts[2].nexts[1].int_data[1] = 2;
  (*elem_node_conn)->nexts[2].nexts[1].int_data[2] = 4;

  //for given point 2 in quad
  (*elem_node_conn)->nexts[2].nexts[2].length = 2 +1; //there are 2 other points
          //near given points 2
  (*elem_node_conn)->nexts[2].nexts[2].int_data = (int *)calloc((2+1),
      sizeof(int));
  (*elem_node_conn)->nexts[2].nexts[2].int_data[1] = 1;
  (*elem_node_conn)->nexts[2].nexts[2].int_data[2] = 3;

  //for given point 3 in quad
  (*elem_node_conn)->nexts[2].nexts[3].length = 2 +1; //there are 2 other points
          //near given points 3
  (*elem_node_conn)->nexts[2].nexts[3].int_data = (int *)calloc((2+1),
      sizeof(int));
  (*elem_node_conn)->nexts[2].nexts[3].int_data[1] = 2;
  (*elem_node_conn)->nexts[2].nexts[3].int_data[2] = 4;

  //for given point 4 in quad
  (*elem_node_conn)->nexts[2].nexts[4].length = 2 +1; //there are 2 other points
          //near given points 4
  (*elem_node_conn)->nexts[2].nexts[4].int_data = (int *)calloc((2+1),
      sizeof(int));
  (*elem_node_conn)->nexts[2].nexts[4].int_data[1] = 1;
  (*elem_node_conn)->nexts[2].nexts[4].int_data[2] = 3;

  //////defining maps for tet///////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  (*elem_node_conn)->nexts[3].length = 4+1; //there are four possibilites for
  //input points in tet
  (*elem_node_conn)->nexts[3].nexts =
      (struct int_vec *)calloc((*elem_node_conn)->nexts[3].length,
      sizeof(struct int_vec));

  //for given point 1 in tet
  (*elem_node_conn)->nexts[3].nexts[1].length = 3 +1; //there are 3 other points
          //near given points 1
  (*elem_node_conn)->nexts[3].nexts[1].int_data = (int *)calloc((3+1),
      sizeof(int));
  (*elem_node_conn)->nexts[3].nexts[1].int_data[1] = 2;
  (*elem_node_conn)->nexts[3].nexts[1].int_data[2] = 3;
  (*elem_node_conn)->nexts[3].nexts[1].int_data[3] = 4;

  //for given point 2 in tet
  (*elem_node_conn)->nexts[3].nexts[2].length = 3 +1; //there are 3 other points
          //near given points 2
  (*elem_node_conn)->nexts[3].nexts[2].int_data = (int *)calloc((3+1),
      sizeof(int));
  (*elem_node_conn)->nexts[3].nexts[2].int_data[1] = 3;
  (*elem_node_conn)->nexts[3].nexts[2].int_data[2] = 4;
  (*elem_node_conn)->nexts[3].nexts[2].int_data[3] = 1;

  //for given point 3 in tet
  (*elem_node_conn)->nexts[3].nexts[3].length = 3 +1; //there are 3 other points
          //near given points 3
  (*elem_node_conn)->nexts[3].nexts[3].int_data = (int *)calloc((3+1),
      sizeof(int));
  (*elem_node_conn)->nexts[3].nexts[3].int_data[1] = 1;
  (*elem_node_conn)->nexts[3].nexts[3].int_data[2] = 2;
  (*elem_node_conn)->nexts[3].nexts[3].int_data[3] = 4;

  //for given point 4 in tet
  (*elem_node_conn)->nexts[3].nexts[4].length = 3 +1; //there are 3 other points
          //near given points 4
  (*elem_node_conn)->nexts[3].nexts[4].int_data = (int *)calloc((3+1),
      sizeof(int));
  (*elem_node_conn)->nexts[3].nexts[4].int_data[1] = 1;
  (*elem_node_conn)->nexts[3].nexts[4].int_data[2] = 2;
  (*elem_node_conn)->nexts[3].nexts[4].int_data[3] = 3;

  //////defining maps for pyramid///////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  (*elem_node_conn)->nexts[4].length = 5+1; //there are five possibilites for
  //input points in pyramid
  (*elem_node_conn)->nexts[4].nexts =
      (struct int_vec *)calloc((*elem_node_conn)->nexts[4].length,
      sizeof(struct int_vec));

  //for given point 1 in pyramid
  (*elem_node_conn)->nexts[4].nexts[1].length = 3 +1; //there are 3 other points
          //near given points 1
  (*elem_node_conn)->nexts[4].nexts[1].int_data = (int *)calloc((3+1),
      sizeof(int));
  (*elem_node_conn)->nexts[4].nexts[1].int_data[1] = 2;
  (*elem_node_conn)->nexts[4].nexts[1].int_data[2] = 5;
  (*elem_node_conn)->nexts[4].nexts[1].int_data[3] = 4;

  //for given point 2 in pyramid
  (*elem_node_conn)->nexts[4].nexts[2].length = 3 +1; //there are 3 other points
          //near given points 2
  (*elem_node_conn)->nexts[4].nexts[2].int_data = (int *)calloc((3+1),
      sizeof(int));
  (*elem_node_conn)->nexts[4].nexts[2].int_data[1] = 3;
  (*elem_node_conn)->nexts[4].nexts[2].int_data[2] = 5;
  (*elem_node_conn)->nexts[4].nexts[2].int_data[3] = 1;

  //for given point 3 in pyramid
  (*elem_node_conn)->nexts[4].nexts[3].length = 3 +1; //there are 3 other points
          //near given points 3
  (*elem_node_conn)->nexts[4].nexts[3].int_data = (int *)calloc((3+1),
      sizeof(int));
  (*elem_node_conn)->nexts[4].nexts[3].int_data[1] = 4;
  (*elem_node_conn)->nexts[4].nexts[3].int_data[2] = 5;
  (*elem_node_conn)->nexts[4].nexts[3].int_data[3] = 2;

  //for given point 4 in pyramid
  (*elem_node_conn)->nexts[4].nexts[4].length = 3 +1; //there are 3 other points
          //near given points 4
  (*elem_node_conn)->nexts[4].nexts[4].int_data = (int *)calloc((3+1),
      sizeof(int));
  (*elem_node_conn)->nexts[4].nexts[4].int_data[1] = 1;
  (*elem_node_conn)->nexts[4].nexts[4].int_data[2] = 5;
  (*elem_node_conn)->nexts[4].nexts[4].int_data[3] = 3;

  //for given point 5 in pyramid
  (*elem_node_conn)->nexts[4].nexts[5].length = 4 +1; //there are 4 other points
          //near given points 5
  (*elem_node_conn)->nexts[4].nexts[5].int_data = (int *)calloc((4+1),
      sizeof(int));
  (*elem_node_conn)->nexts[4].nexts[5].int_data[1] = 1;
  (*elem_node_conn)->nexts[4].nexts[5].int_data[2] = 2;
  (*elem_node_conn)->nexts[4].nexts[5].int_data[3] = 3;
  (*elem_node_conn)->nexts[4].nexts[5].int_data[4] = 4;

  //////defining maps for prism/////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  (*elem_node_conn)->nexts[5].length = 6+1; //there are six possibilites for
  //input points in prism
  (*elem_node_conn)->nexts[5].nexts =
      (struct int_vec *)calloc((*elem_node_conn)->nexts[5].length,
      sizeof(struct int_vec));

  //for given point 1 in prism
  (*elem_node_conn)->nexts[5].nexts[1].length = 3 +1; //there are 3 other points
          //near given points 1
  (*elem_node_conn)->nexts[5].nexts[1].int_data = (int *)calloc((3+1),
      sizeof(int));
  (*elem_node_conn)->nexts[5].nexts[1].int_data[1] = 2;
  (*elem_node_conn)->nexts[5].nexts[1].int_data[2] = 3;
  (*elem_node_conn)->nexts[5].nexts[1].int_data[3] = 4;

  //for given point 2 in prism
  (*elem_node_conn)->nexts[5].nexts[2].length = 3 +1; //there are 3 other points
          //near given points 2
  (*elem_node_conn)->nexts[5].nexts[2].int_data = (int *)calloc((3+1),
      sizeof(int));
  (*elem_node_conn)->nexts[5].nexts[2].int_data[1] = 3;
  (*elem_node_conn)->nexts[5].nexts[2].int_data[2] = 5;
  (*elem_node_conn)->nexts[5].nexts[2].int_data[3] = 1;

  //for given point 3 in prism
  (*elem_node_conn)->nexts[5].nexts[3].length = 3 +1; //there are 3 other points
          //near given points 3
  (*elem_node_conn)->nexts[5].nexts[3].int_data = (int *)calloc((3+1),
      sizeof(int));
  (*elem_node_conn)->nexts[5].nexts[3].int_data[1] = 6;
  (*elem_node_conn)->nexts[5].nexts[3].int_data[2] = 2;
  (*elem_node_conn)->nexts[5].nexts[3].int_data[3] = 1;

  //for given point 4 in prism
  (*elem_node_conn)->nexts[5].nexts[4].length = 3 +1; //there are 3 other points
          //near given points 4
  (*elem_node_conn)->nexts[5].nexts[4].int_data = (int *)calloc((3+1),
      sizeof(int));
  (*elem_node_conn)->nexts[5].nexts[4].int_data[1] = 1;
  (*elem_node_conn)->nexts[5].nexts[4].int_data[2] = 5;
  (*elem_node_conn)->nexts[5].nexts[4].int_data[3] = 6;

  //for given point 5 in prism
  (*elem_node_conn)->nexts[5].nexts[5].length = 3 +1; //there are 3 other points
          //near given points 5
  (*elem_node_conn)->nexts[5].nexts[5].int_data = (int *)calloc((3+1),
      sizeof(int));
  (*elem_node_conn)->nexts[5].nexts[5].int_data[1] = 2;
  (*elem_node_conn)->nexts[5].nexts[5].int_data[2] = 6;
  (*elem_node_conn)->nexts[5].nexts[5].int_data[3] = 4;

  //for given point 6 in prism
  (*elem_node_conn)->nexts[5].nexts[6].length = 3 +1; //there are 3 other points
          //near given points 6
  (*elem_node_conn)->nexts[5].nexts[6].int_data = (int *)calloc((3+1),
      sizeof(int));
  (*elem_node_conn)->nexts[5].nexts[6].int_data[1] = 4;
  (*elem_node_conn)->nexts[5].nexts[6].int_data[2] = 5;
  (*elem_node_conn)->nexts[5].nexts[6].int_data[3] = 3;

  //////defining maps for hex/////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  (*elem_node_conn)->nexts[6].length = 8+1; //there are eight possibilites for
  //input points in hex
  (*elem_node_conn)->nexts[6].nexts =
      (struct int_vec *)calloc((*elem_node_conn)->nexts[6].length,
      sizeof(struct int_vec));

  //for given point 1 in hex
  (*elem_node_conn)->nexts[6].nexts[1].length = 3 +1; //there are 3 other points
          //near given points 1
  (*elem_node_conn)->nexts[6].nexts[1].int_data = (int *)calloc((3+1),
      sizeof(int));
  (*elem_node_conn)->nexts[6].nexts[1].int_data[1] = 2;
  (*elem_node_conn)->nexts[6].nexts[1].int_data[2] = 4;
  (*elem_node_conn)->nexts[6].nexts[1].int_data[3] = 5;

  //for given point 2 in hex
  (*elem_node_conn)->nexts[6].nexts[2].length = 3 +1; //there are 3 other points
          //near given points 2
  (*elem_node_conn)->nexts[6].nexts[2].int_data = (int *)calloc((3+1),
      sizeof(int));
  (*elem_node_conn)->nexts[6].nexts[2].int_data[1] = 3;
  (*elem_node_conn)->nexts[6].nexts[2].int_data[2] = 6;
  (*elem_node_conn)->nexts[6].nexts[2].int_data[3] = 1;

  //for given point 3 in hex
  (*elem_node_conn)->nexts[6].nexts[3].length = 3 +1; //there are 3 other points
          //near given points 3
  (*elem_node_conn)->nexts[6].nexts[3].int_data = (int *)calloc((3+1),
      sizeof(int));
  (*elem_node_conn)->nexts[6].nexts[3].int_data[1] = 2;
  (*elem_node_conn)->nexts[6].nexts[3].int_data[2] = 4;
  (*elem_node_conn)->nexts[6].nexts[3].int_data[3] = 7;

  //for given point 4 in hex
  (*elem_node_conn)->nexts[6].nexts[4].length = 3 +1; //there are 3 other points
          //near given points 4
  (*elem_node_conn)->nexts[6].nexts[4].int_data = (int *)calloc((3+1),
      sizeof(int));
  (*elem_node_conn)->nexts[6].nexts[4].int_data[1] = 1;
  (*elem_node_conn)->nexts[6].nexts[4].int_data[2] = 3;
  (*elem_node_conn)->nexts[6].nexts[4].int_data[3] = 8;

  //for given point 5 in hex
  (*elem_node_conn)->nexts[6].nexts[5].length = 3 +1; //there are 3 other points
          //near given points 5
  (*elem_node_conn)->nexts[6].nexts[5].int_data = (int *)calloc((3+1),
      sizeof(int));
  (*elem_node_conn)->nexts[6].nexts[5].int_data[1] = 1;
  (*elem_node_conn)->nexts[6].nexts[5].int_data[2] = 6;
  (*elem_node_conn)->nexts[6].nexts[5].int_data[3] = 8;

  //for given point 6 in hex
  (*elem_node_conn)->nexts[6].nexts[6].length = 3 +1; //there are 3 other points
          //near given points 6
  (*elem_node_conn)->nexts[6].nexts[6].int_data = (int *)calloc((3+1),
      sizeof(int));
  (*elem_node_conn)->nexts[6].nexts[6].int_data[1] = 5;
  (*elem_node_conn)->nexts[6].nexts[6].int_data[2] = 2;
  (*elem_node_conn)->nexts[6].nexts[6].int_data[3] = 7;

   //for given point 7 in hex
  (*elem_node_conn)->nexts[6].nexts[7].length = 3 +1; //there are 3 other points
          //near given points 7
  (*elem_node_conn)->nexts[6].nexts[7].int_data = (int *)calloc((3+1),
      sizeof(int));
  (*elem_node_conn)->nexts[6].nexts[7].int_data[1] = 8;
  (*elem_node_conn)->nexts[6].nexts[7].int_data[2] = 6;
  (*elem_node_conn)->nexts[6].nexts[7].int_data[3] = 3;

  //for given point 8 in hex
  (*elem_node_conn)->nexts[6].nexts[8].length = 3 +1; //there are 3 other points
          //near given points 8
  (*elem_node_conn)->nexts[6].nexts[8].int_data = (int *)calloc((3+1),
      sizeof(int));
  (*elem_node_conn)->nexts[6].nexts[8].int_data[1] = 5;
  (*elem_node_conn)->nexts[6].nexts[8].int_data[2] = 7;
  (*elem_node_conn)->nexts[6].nexts[8].int_data[3] = 4;

  ///definition of element-node look-up table completed.

  return 0; //return successfully!
}

///////////////////////////////////////////////////////////////////
//    Surrounding faces to element node (surr_face_el_node)
//The following creates a structure containing the tags of faces surr
//a node in the given element type. Everything is regarded local
// surr_face_el_node -> etype[]->[node].{face tag1,face tag2, ..., face tagn}
int Create_Surr_Face_El_Node(struct int_vec **surr_face_el_node)
{
  (*surr_face_el_node) = (struct int_vec *)calloc(1,sizeof(struct int_vec));
  (*surr_face_el_node)->length = 6+1; //tri,quad,tet,pyramid,prism,hex
  //allocating memory for different element types
  (*surr_face_el_node)->nexts = (struct int_vec *)calloc((*surr_face_el_node)->length,
                                                sizeof(struct int_vec));
  //////defining maps for triangle//////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  //------Triangle
  (*surr_face_el_node)->nexts[1].length = 3+1; //there are three possibilites for
  //input points in triangles
  (*surr_face_el_node)->nexts[1].nexts =
      (struct int_vec *)calloc((*surr_face_el_node)->nexts[1].length,
      sizeof(struct int_vec));
  //there is one face near point 1
  (*surr_face_el_node)->nexts[1].nexts[1].length = 1+1;
  (*surr_face_el_node)->nexts[1].nexts[1].int_data = (int *)calloc((2),
      sizeof(int));
  (*surr_face_el_node)->nexts[1].nexts[1].int_data[1] = 1; //the tag of the face

  //there is one face near point 2 in triangle
  (*surr_face_el_node)->nexts[1].nexts[2].length = 1+1;
  (*surr_face_el_node)->nexts[1].nexts[2].int_data = (int *)calloc((2),
      sizeof(int));
  (*surr_face_el_node)->nexts[1].nexts[2].int_data[1] = 1; //the tag of the face

  //there is one face near point 3 in triangle
  (*surr_face_el_node)->nexts[1].nexts[3].length = 1+1;
  (*surr_face_el_node)->nexts[1].nexts[3].int_data = (int *)calloc((2),
      sizeof(int));
  (*surr_face_el_node)->nexts[1].nexts[3].int_data[1] = 1; //the tag of the face

  //------Quad
  (*surr_face_el_node)->nexts[2].length = 4+1; //there are four possibilites for
  //input points in quad
  (*surr_face_el_node)->nexts[2].nexts =
      (struct int_vec *)calloc((*surr_face_el_node)->nexts[2].length,
      sizeof(struct int_vec));
  //there is one face near point 1
  (*surr_face_el_node)->nexts[2].nexts[1].length = 1+1;
  (*surr_face_el_node)->nexts[2].nexts[1].int_data = (int *)calloc((2),
      sizeof(int));
  (*surr_face_el_node)->nexts[2].nexts[1].int_data[1] = 1; //the tag of the face

  //there is one face near point 2
  (*surr_face_el_node)->nexts[2].nexts[2].length = 1+1;
  (*surr_face_el_node)->nexts[2].nexts[2].int_data = (int *)calloc((2),
      sizeof(int));
  (*surr_face_el_node)->nexts[2].nexts[2].int_data[1] = 1; //the tag of the face

  //there is one face near point 3
  (*surr_face_el_node)->nexts[2].nexts[3].length = 1+1;
  (*surr_face_el_node)->nexts[2].nexts[3].int_data = (int *)calloc((2),
      sizeof(int));
  (*surr_face_el_node)->nexts[2].nexts[3].int_data[1] = 1; //the tag of the face

  //there is one face near point 4
  (*surr_face_el_node)->nexts[2].nexts[4].length = 1+1;
  (*surr_face_el_node)->nexts[2].nexts[4].int_data = (int *)calloc((2),
      sizeof(int));
  (*surr_face_el_node)->nexts[2].nexts[4].int_data[1] = 1; //the tag of the face

//-----Tet
    (*surr_face_el_node)->nexts[3].length = 4+1; //there are four possibilites for
  //input points in tet
  (*surr_face_el_node)->nexts[3].nexts =
      (struct int_vec *)calloc((*surr_face_el_node)->nexts[3].length,
      sizeof(struct int_vec));
  //there are three faces near point 1
  (*surr_face_el_node)->nexts[3].nexts[1].length = 3+1;
  (*surr_face_el_node)->nexts[3].nexts[1].int_data = (int *)calloc((4),
      sizeof(int));
  (*surr_face_el_node)->nexts[3].nexts[1].int_data[1] = 1; //the tag of the face
  (*surr_face_el_node)->nexts[3].nexts[1].int_data[2] = 2; //the tag of the face
  (*surr_face_el_node)->nexts[3].nexts[1].int_data[3] = 3; //the tag of the face

  //there are three faces near point 2
  (*surr_face_el_node)->nexts[3].nexts[2].length = 3+1;
  (*surr_face_el_node)->nexts[3].nexts[2].int_data = (int *)calloc((4),
      sizeof(int));
  (*surr_face_el_node)->nexts[3].nexts[2].int_data[1] = 1; //the tag of the face
  (*surr_face_el_node)->nexts[3].nexts[2].int_data[2] = 2; //the tag of the face
  (*surr_face_el_node)->nexts[3].nexts[2].int_data[3] = 4; //the tag of the face

  //there are three faces near point 3
  (*surr_face_el_node)->nexts[3].nexts[3].length = 3+1;
  (*surr_face_el_node)->nexts[3].nexts[3].int_data = (int *)calloc((4),
      sizeof(int));
  (*surr_face_el_node)->nexts[3].nexts[3].int_data[1] = 1; //the tag of the face
  (*surr_face_el_node)->nexts[3].nexts[3].int_data[2] = 3; //the tag of the face
  (*surr_face_el_node)->nexts[3].nexts[3].int_data[3] = 4; //the tag of the face

  //there are three faces near point 4
  (*surr_face_el_node)->nexts[3].nexts[4].length = 3+1;
  (*surr_face_el_node)->nexts[3].nexts[4].int_data = (int *)calloc((4),
      sizeof(int));
  (*surr_face_el_node)->nexts[3].nexts[4].int_data[1] = 2; //the tag of the face
  (*surr_face_el_node)->nexts[3].nexts[4].int_data[2] = 3; //the tag of the face
  (*surr_face_el_node)->nexts[3].nexts[4].int_data[3] = 4; //the tag of the face

  //-----pyramid
  (*surr_face_el_node)->nexts[4].length = 5+1; //there are five possibilites for
  //input points in pyramid
  (*surr_face_el_node)->nexts[4].nexts =
      (struct int_vec *)calloc((*surr_face_el_node)->nexts[4].length,
      sizeof(struct int_vec));
  //there are three faces near point 1
  (*surr_face_el_node)->nexts[4].nexts[1].length = 3+1;
  (*surr_face_el_node)->nexts[4].nexts[1].int_data = (int *)calloc((4),
      sizeof(int));
  (*surr_face_el_node)->nexts[4].nexts[1].int_data[1] = 1; //the tag of the face
  (*surr_face_el_node)->nexts[4].nexts[1].int_data[2] = 2; //the tag of the face
  (*surr_face_el_node)->nexts[4].nexts[1].int_data[3] = 3; //the tag of the face

  //there are three faces near point 2
  (*surr_face_el_node)->nexts[4].nexts[2].length = 3+1;
  (*surr_face_el_node)->nexts[4].nexts[2].int_data = (int *)calloc((4),
      sizeof(int));
  (*surr_face_el_node)->nexts[4].nexts[2].int_data[1] = 1; //the tag of the face
  (*surr_face_el_node)->nexts[4].nexts[2].int_data[2] = 2; //the tag of the face
  (*surr_face_el_node)->nexts[4].nexts[2].int_data[3] = 4; //the tag of the face

  //there are three faces near point 3
  (*surr_face_el_node)->nexts[4].nexts[3].length = 3+1;
  (*surr_face_el_node)->nexts[4].nexts[3].int_data = (int *)calloc((4),
      sizeof(int));
  (*surr_face_el_node)->nexts[4].nexts[3].int_data[1] = 1; //the tag of the face
  (*surr_face_el_node)->nexts[4].nexts[3].int_data[2] = 4; //the tag of the face
  (*surr_face_el_node)->nexts[4].nexts[3].int_data[3] = 5; //the tag of the face

  //there are three faces near point 4
  (*surr_face_el_node)->nexts[4].nexts[4].length = 3+1;
  (*surr_face_el_node)->nexts[4].nexts[4].int_data = (int *)calloc((4),
      sizeof(int));
  (*surr_face_el_node)->nexts[4].nexts[4].int_data[1] = 1; //the tag of the face
  (*surr_face_el_node)->nexts[4].nexts[4].int_data[2] = 3; //the tag of the face
  (*surr_face_el_node)->nexts[4].nexts[4].int_data[3] = 5; //the tag of the face

  //there are three faces near point 5
  (*surr_face_el_node)->nexts[4].nexts[5].length = 4+1;
  (*surr_face_el_node)->nexts[4].nexts[5].int_data = (int *)calloc((5),
      sizeof(int));
  (*surr_face_el_node)->nexts[4].nexts[5].int_data[1] = 2; //the tag of the face
  (*surr_face_el_node)->nexts[4].nexts[5].int_data[2] = 3; //the tag of the face
  (*surr_face_el_node)->nexts[4].nexts[5].int_data[3] = 4; //the tag of the face
  (*surr_face_el_node)->nexts[4].nexts[5].int_data[4] = 5; //the tag of the face

  //-----prism
  (*surr_face_el_node)->nexts[5].length = 6+1; //there are six possibilites for
  //input points in prism
  (*surr_face_el_node)->nexts[5].nexts =
      (struct int_vec *)calloc((*surr_face_el_node)->nexts[5].length,
      sizeof(struct int_vec));
  //there are three faces near point 1
  (*surr_face_el_node)->nexts[5].nexts[1].length = 3+1;
  (*surr_face_el_node)->nexts[5].nexts[1].int_data = (int *)calloc((4),
      sizeof(int));
  (*surr_face_el_node)->nexts[5].nexts[1].int_data[1] = 1; //the tag of the face
  (*surr_face_el_node)->nexts[5].nexts[1].int_data[2] = 2; //the tag of the face
  (*surr_face_el_node)->nexts[5].nexts[1].int_data[3] = 3; //the tag of the face

  //there are three faces near point 2
  (*surr_face_el_node)->nexts[5].nexts[2].length = 3+1;
  (*surr_face_el_node)->nexts[5].nexts[2].int_data = (int *)calloc((4),
      sizeof(int));
  (*surr_face_el_node)->nexts[5].nexts[2].int_data[1] = 1; //the tag of the face
  (*surr_face_el_node)->nexts[5].nexts[2].int_data[2] = 2; //the tag of the face
  (*surr_face_el_node)->nexts[5].nexts[2].int_data[3] = 4; //the tag of the face

  //there are three faces near point 3
  (*surr_face_el_node)->nexts[5].nexts[3].length = 3+1;
  (*surr_face_el_node)->nexts[5].nexts[3].int_data = (int *)calloc((4),
      sizeof(int));
  (*surr_face_el_node)->nexts[5].nexts[3].int_data[1] = 1; //the tag of the face
  (*surr_face_el_node)->nexts[5].nexts[3].int_data[2] = 3; //the tag of the face
  (*surr_face_el_node)->nexts[5].nexts[3].int_data[3] = 4; //the tag of the face

  //there are three faces near point 4
  (*surr_face_el_node)->nexts[5].nexts[4].length = 3+1;
  (*surr_face_el_node)->nexts[5].nexts[4].int_data = (int *)calloc((4),
      sizeof(int));
  (*surr_face_el_node)->nexts[5].nexts[4].int_data[1] = 2; //the tag of the face
  (*surr_face_el_node)->nexts[5].nexts[4].int_data[2] = 3; //the tag of the face
  (*surr_face_el_node)->nexts[5].nexts[4].int_data[3] = 5; //the tag of the face

  //there are three faces near point 5
  (*surr_face_el_node)->nexts[5].nexts[5].length = 3+1;
  (*surr_face_el_node)->nexts[5].nexts[5].int_data = (int *)calloc((4),
      sizeof(int));
  (*surr_face_el_node)->nexts[5].nexts[5].int_data[1] = 4; //the tag of the face
  (*surr_face_el_node)->nexts[5].nexts[5].int_data[2] = 2; //the tag of the face
  (*surr_face_el_node)->nexts[5].nexts[5].int_data[3] = 5; //the tag of the face

  //there are three faces near point 6
  (*surr_face_el_node)->nexts[5].nexts[6].length = 3+1;
  (*surr_face_el_node)->nexts[5].nexts[6].int_data = (int *)calloc((4),
      sizeof(int));
  (*surr_face_el_node)->nexts[5].nexts[6].int_data[1] = 5; //the tag of the face
  (*surr_face_el_node)->nexts[5].nexts[6].int_data[2] = 4; //the tag of the face
  (*surr_face_el_node)->nexts[5].nexts[6].int_data[3] = 3; //the tag of the face

  //-----hex
  (*surr_face_el_node)->nexts[6].length = 8+1;//there are eight possibilites for
  //input points in hex
  (*surr_face_el_node)->nexts[6].nexts =
      (struct int_vec *)calloc((*surr_face_el_node)->nexts[6].length,
      sizeof(struct int_vec));
  //there are three faces near point 1
  (*surr_face_el_node)->nexts[6].nexts[1].length = 3+1;
  (*surr_face_el_node)->nexts[6].nexts[1].int_data = (int *)calloc((4),
      sizeof(int));
  (*surr_face_el_node)->nexts[6].nexts[1].int_data[1] = 1; //the tag of the face
  (*surr_face_el_node)->nexts[6].nexts[1].int_data[2] = 2; //the tag of the face
  (*surr_face_el_node)->nexts[6].nexts[1].int_data[3] = 3; //the tag of the face

  //there are three faces near point 2
  (*surr_face_el_node)->nexts[6].nexts[2].length = 3+1;
  (*surr_face_el_node)->nexts[6].nexts[2].int_data = (int *)calloc((4),
      sizeof(int));
  (*surr_face_el_node)->nexts[6].nexts[2].int_data[1] = 1; //the tag of the face
  (*surr_face_el_node)->nexts[6].nexts[2].int_data[2] = 2; //the tag of the face
  (*surr_face_el_node)->nexts[6].nexts[2].int_data[3] = 4; //the tag of the face

  //there are three faces near point 3
  (*surr_face_el_node)->nexts[6].nexts[3].length = 3+1;
  (*surr_face_el_node)->nexts[6].nexts[3].int_data = (int *)calloc((4),
      sizeof(int));
  (*surr_face_el_node)->nexts[6].nexts[3].int_data[1] = 1; //the tag of the face
  (*surr_face_el_node)->nexts[6].nexts[3].int_data[2] = 4; //the tag of the face
  (*surr_face_el_node)->nexts[6].nexts[3].int_data[3] = 5; //the tag of the face

  //there are three faces near point 4
  (*surr_face_el_node)->nexts[6].nexts[4].length = 3+1;
  (*surr_face_el_node)->nexts[6].nexts[4].int_data = (int *)calloc((4),
      sizeof(int));
  (*surr_face_el_node)->nexts[6].nexts[4].int_data[1] = 1; //the tag of the face
  (*surr_face_el_node)->nexts[6].nexts[4].int_data[2] = 3; //the tag of the face
  (*surr_face_el_node)->nexts[6].nexts[4].int_data[3] = 5; //the tag of the face

  //there are three faces near point 5
  (*surr_face_el_node)->nexts[6].nexts[5].length = 3+1;
  (*surr_face_el_node)->nexts[6].nexts[5].int_data = (int *)calloc((4),
      sizeof(int));
  (*surr_face_el_node)->nexts[6].nexts[5].int_data[1] = 2; //the tag of the face
  (*surr_face_el_node)->nexts[6].nexts[5].int_data[2] = 3; //the tag of the face
  (*surr_face_el_node)->nexts[6].nexts[5].int_data[3] = 6; //the tag of the face

  //there are three faces near point 6
  (*surr_face_el_node)->nexts[6].nexts[6].length = 3+1;
  (*surr_face_el_node)->nexts[6].nexts[6].int_data = (int *)calloc((4),
      sizeof(int));
  (*surr_face_el_node)->nexts[6].nexts[6].int_data[1] = 2; //the tag of the face
  (*surr_face_el_node)->nexts[6].nexts[6].int_data[2] = 4; //the tag of the face
  (*surr_face_el_node)->nexts[6].nexts[6].int_data[3] = 6; //the tag of the face

  //there are three faces near point 7
  (*surr_face_el_node)->nexts[6].nexts[7].length = 3+1;
  (*surr_face_el_node)->nexts[6].nexts[7].int_data = (int *)calloc((4),
      sizeof(int));
  (*surr_face_el_node)->nexts[6].nexts[7].int_data[1] = 5; //the tag of the face
  (*surr_face_el_node)->nexts[6].nexts[7].int_data[2] = 4; //the tag of the face
  (*surr_face_el_node)->nexts[6].nexts[7].int_data[3] = 6; //the tag of the face

  //there are three faces near point 8
  (*surr_face_el_node)->nexts[6].nexts[8].length = 3+1;
  (*surr_face_el_node)->nexts[6].nexts[8].int_data = (int *)calloc((4),
      sizeof(int));
  (*surr_face_el_node)->nexts[6].nexts[8].int_data[1] = 5; //the tag of the face
  (*surr_face_el_node)->nexts[6].nexts[8].int_data[2] = 3; //the tag of the face
  (*surr_face_el_node)->nexts[6].nexts[8].int_data[3] = 6; //the tag of the face

  //completed successfully!
    return 0;
}

///////////////////////////////////////////////////////////////////
//    points in a face in an element map (pnts_face_el)
//The following creates a structure containing the points of faces in elem type
//Everything is regarded local
// pnts_face_el -> etype[]->[face].{local ptn1,local ptn2, ..., local ptn,m}
int Create_Points_In_Face_El(struct int_vec **pnts_face_el)
{
  (*pnts_face_el) = (struct int_vec *)calloc(1,sizeof(struct int_vec));
  (*pnts_face_el)->length = 6+1; //tri,quad,tet,pyramid,prism,hex
  //allocating memory for different element types
  (*pnts_face_el)->nexts = (struct int_vec *)calloc((*pnts_face_el)->length,
                                                sizeof(struct int_vec));
  //////defining maps for triangle//////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  //------Triangle
  (*pnts_face_el)->nexts[1].length = 1+1; //there are one possibility for
  //input points in triangles because there is only one surface
  (*pnts_face_el)->nexts[1].nexts =
      (struct int_vec *)calloc((*pnts_face_el)->nexts[1].length,
      sizeof(struct int_vec));
  //there are three points on the face 1
  (*pnts_face_el)->nexts[1].nexts[1].length = 3+1;
  (*pnts_face_el)->nexts[1].nexts[1].int_data = (int *)calloc((4),
      sizeof(int));
  (*pnts_face_el)->nexts[1].nexts[1].int_data[1] = 1; //the local index of point
  (*pnts_face_el)->nexts[1].nexts[1].int_data[2] = 2; //the local index of point
  (*pnts_face_el)->nexts[1].nexts[1].int_data[3] = 3; //the local index of point

  //------Quad
  (*pnts_face_el)->nexts[2].length = 1+1; //there are one possibility for
  //input points in quadlateral because there is only one face
  (*pnts_face_el)->nexts[2].nexts =
      (struct int_vec *)calloc((*pnts_face_el)->nexts[2].length,
      sizeof(struct int_vec));
  //there are four points on the face 1
  (*pnts_face_el)->nexts[2].nexts[1].length = 4+1;
  (*pnts_face_el)->nexts[2].nexts[1].int_data = (int *)calloc((5),
      sizeof(int));
  (*pnts_face_el)->nexts[2].nexts[1].int_data[1] = 1; //the local index of point
  (*pnts_face_el)->nexts[2].nexts[1].int_data[2] = 2; //the local index of point
  (*pnts_face_el)->nexts[2].nexts[1].int_data[3] = 3; //the local index of point
  (*pnts_face_el)->nexts[2].nexts[1].int_data[4] = 4; //the local index of point

  //------tet
  (*pnts_face_el)->nexts[3].length = 4+1; //there are four possibility for
  //input points in tetrahedral because there are four faces
  (*pnts_face_el)->nexts[3].nexts =
      (struct int_vec *)calloc((*pnts_face_el)->nexts[3].length,
      sizeof(struct int_vec));
  //there are three points on the face 1
  (*pnts_face_el)->nexts[3].nexts[1].length = 3+1;
  (*pnts_face_el)->nexts[3].nexts[1].int_data = (int *)calloc((4),
      sizeof(int));
  (*pnts_face_el)->nexts[3].nexts[1].int_data[1] = 1; //the local index of point
  (*pnts_face_el)->nexts[3].nexts[1].int_data[2] = 2; //the local index of point
  (*pnts_face_el)->nexts[3].nexts[1].int_data[3] = 3; //the local index of point

  //there are three points on the face 2
  (*pnts_face_el)->nexts[3].nexts[2].length = 3+1;
  (*pnts_face_el)->nexts[3].nexts[2].int_data = (int *)calloc((4),
      sizeof(int));
  (*pnts_face_el)->nexts[3].nexts[2].int_data[1] = 1; //the local index of point
  (*pnts_face_el)->nexts[3].nexts[2].int_data[2] = 2; //the local index of point
  (*pnts_face_el)->nexts[3].nexts[2].int_data[3] = 4; //the local index of point

  //there are three points on the face 3
  (*pnts_face_el)->nexts[3].nexts[3].length = 3+1;
  (*pnts_face_el)->nexts[3].nexts[3].int_data = (int *)calloc((4),
      sizeof(int));
  (*pnts_face_el)->nexts[3].nexts[3].int_data[1] = 1; //the local index of point
  (*pnts_face_el)->nexts[3].nexts[3].int_data[2] = 3; //the local index of point
  (*pnts_face_el)->nexts[3].nexts[3].int_data[3] = 4; //the local index of point

  //there are three points on the face 4
  (*pnts_face_el)->nexts[3].nexts[4].length = 3+1;
  (*pnts_face_el)->nexts[3].nexts[4].int_data = (int *)calloc((4),
      sizeof(int));
  (*pnts_face_el)->nexts[3].nexts[4].int_data[1] = 2; //the local index of point
  (*pnts_face_el)->nexts[3].nexts[4].int_data[2] = 3; //the local index of point
  (*pnts_face_el)->nexts[3].nexts[4].int_data[3] = 4; //the local index of point

  //------Pyramid
  (*pnts_face_el)->nexts[4].length = 5+1; //there are five possibilities for
  //input faces in tetrahedral because there are fives faces
  (*pnts_face_el)->nexts[4].nexts =
      (struct int_vec *)calloc((*pnts_face_el)->nexts[4].length,
      sizeof(struct int_vec));
  //there are four points on the face 1
  (*pnts_face_el)->nexts[4].nexts[1].length = 4+1;
  (*pnts_face_el)->nexts[4].nexts[1].int_data = (int *)calloc((5),
      sizeof(int));
  (*pnts_face_el)->nexts[4].nexts[1].int_data[1] = 1; //the local index of point
  (*pnts_face_el)->nexts[4].nexts[1].int_data[2] = 2; //the local index of point
  (*pnts_face_el)->nexts[4].nexts[1].int_data[3] = 3; //the local index of point
  (*pnts_face_el)->nexts[4].nexts[1].int_data[4] = 4; //the local index of point

  //there are three points on the face 2
  (*pnts_face_el)->nexts[4].nexts[2].length = 3+1;
  (*pnts_face_el)->nexts[4].nexts[2].int_data = (int *)calloc((4),
      sizeof(int));
  (*pnts_face_el)->nexts[4].nexts[2].int_data[1] = 1; //the local index of point
  (*pnts_face_el)->nexts[4].nexts[2].int_data[2] = 2; //the local index of point
  (*pnts_face_el)->nexts[4].nexts[2].int_data[3] = 5; //the local index of point

  //there are three points on the face 3
  (*pnts_face_el)->nexts[4].nexts[3].length = 3+1;
  (*pnts_face_el)->nexts[4].nexts[3].int_data = (int *)calloc((4),
      sizeof(int));
  (*pnts_face_el)->nexts[4].nexts[3].int_data[1] = 1; //the local index of point
  (*pnts_face_el)->nexts[4].nexts[3].int_data[2] = 4; //the local index of point
  (*pnts_face_el)->nexts[4].nexts[3].int_data[3] = 5; //the local index of point

  //there are three points on the face 4
  (*pnts_face_el)->nexts[4].nexts[4].length = 3+1;
  (*pnts_face_el)->nexts[4].nexts[4].int_data = (int *)calloc((4),
      sizeof(int));
  (*pnts_face_el)->nexts[4].nexts[4].int_data[1] = 2; //the local index of point
  (*pnts_face_el)->nexts[4].nexts[4].int_data[2] = 3; //the local index of point
  (*pnts_face_el)->nexts[4].nexts[4].int_data[3] = 5; //the local index of point

  //there are three points on the face 5
  (*pnts_face_el)->nexts[4].nexts[5].length = 3+1;
  (*pnts_face_el)->nexts[4].nexts[5].int_data = (int *)calloc((4),
      sizeof(int));
  (*pnts_face_el)->nexts[4].nexts[5].int_data[1] = 3; //the local index of point
  (*pnts_face_el)->nexts[4].nexts[5].int_data[2] = 4; //the local index of point
  (*pnts_face_el)->nexts[4].nexts[5].int_data[3] = 5; //the local index of point

  //------Prism
  (*pnts_face_el)->nexts[5].length = 5+1; //there are five possibilities for
  //input faces in tetrahedral because there are fives faces
  (*pnts_face_el)->nexts[5].nexts =
      (struct int_vec *)calloc((*pnts_face_el)->nexts[5].length,
      sizeof(struct int_vec));
  //there are four points on the face 1
  (*pnts_face_el)->nexts[5].nexts[1].length = 3+1;
  (*pnts_face_el)->nexts[5].nexts[1].int_data = (int *)calloc((4),
      sizeof(int));
  (*pnts_face_el)->nexts[5].nexts[1].int_data[1] = 1; //the local index of point
  (*pnts_face_el)->nexts[5].nexts[1].int_data[2] = 2; //the local index of point
  (*pnts_face_el)->nexts[5].nexts[1].int_data[3] = 3; //the local index of point

  //there are four points on the face 2
  (*pnts_face_el)->nexts[5].nexts[2].length = 4+1;
  (*pnts_face_el)->nexts[5].nexts[2].int_data = (int *)calloc((5),
      sizeof(int));
  (*pnts_face_el)->nexts[5].nexts[2].int_data[1] = 1; //the local index of point
  (*pnts_face_el)->nexts[5].nexts[2].int_data[2] = 2; //the local index of point
  (*pnts_face_el)->nexts[5].nexts[2].int_data[3] = 4; //the local index of point
  (*pnts_face_el)->nexts[5].nexts[2].int_data[4] = 5; //the local index of point

  //there are four points on the face 3
  (*pnts_face_el)->nexts[5].nexts[3].length = 4+1;
  (*pnts_face_el)->nexts[5].nexts[3].int_data = (int *)calloc((5),
      sizeof(int));
  (*pnts_face_el)->nexts[5].nexts[3].int_data[1] = 1; //the local index of point
  (*pnts_face_el)->nexts[5].nexts[3].int_data[2] = 4; //the local index of point
  (*pnts_face_el)->nexts[5].nexts[3].int_data[3] = 6; //the local index of point
  (*pnts_face_el)->nexts[5].nexts[3].int_data[4] = 3; //the local index of point

  //there are four points on the face 4
  (*pnts_face_el)->nexts[5].nexts[4].length = 4+1;
  (*pnts_face_el)->nexts[5].nexts[4].int_data = (int *)calloc((5),
      sizeof(int));
  (*pnts_face_el)->nexts[5].nexts[4].int_data[1] = 2; //the local index of point
  (*pnts_face_el)->nexts[5].nexts[4].int_data[2] = 3; //the local index of point
  (*pnts_face_el)->nexts[5].nexts[4].int_data[3] = 5; //the local index of point
  (*pnts_face_el)->nexts[5].nexts[4].int_data[4] = 6; //the local index of point

  //there are four points on the face 5
  (*pnts_face_el)->nexts[5].nexts[5].length = 3+1;
  (*pnts_face_el)->nexts[5].nexts[5].int_data = (int *)calloc((4),
      sizeof(int));
  (*pnts_face_el)->nexts[5].nexts[5].int_data[1] = 4; //the local index of point
  (*pnts_face_el)->nexts[5].nexts[5].int_data[2] = 5; //the local index of point
  (*pnts_face_el)->nexts[5].nexts[5].int_data[3] = 6; //the local index of point

  //------Hex
  (*pnts_face_el)->nexts[6].length = 6+1; //there are six possibilities for
  //input faces in hex because there are six faces
  (*pnts_face_el)->nexts[6].nexts =
      (struct int_vec *)calloc((*pnts_face_el)->nexts[6].length,
      sizeof(struct int_vec));
  //there are four points on the face 1
  (*pnts_face_el)->nexts[6].nexts[1].length = 4+1;
  (*pnts_face_el)->nexts[6].nexts[1].int_data = (int *)calloc((5),
      sizeof(int));
  (*pnts_face_el)->nexts[6].nexts[1].int_data[1] = 1; //the local index of point
  (*pnts_face_el)->nexts[6].nexts[1].int_data[2] = 2; //the local index of point
  (*pnts_face_el)->nexts[6].nexts[1].int_data[3] = 3; //the local index of point
  (*pnts_face_el)->nexts[6].nexts[1].int_data[4] = 4; //the local index of point

  //there are four points on the face 2
  (*pnts_face_el)->nexts[6].nexts[2].length = 4+1;
  (*pnts_face_el)->nexts[6].nexts[2].int_data = (int *)calloc((5),
      sizeof(int));
  (*pnts_face_el)->nexts[6].nexts[2].int_data[1] = 1; //the local index of point
  (*pnts_face_el)->nexts[6].nexts[2].int_data[2] = 2; //the local index of point
  (*pnts_face_el)->nexts[6].nexts[2].int_data[3] = 5; //the local index of point
  (*pnts_face_el)->nexts[6].nexts[2].int_data[4] = 6; //the local index of point

  //there are four points on the face 3
  (*pnts_face_el)->nexts[6].nexts[3].length = 4+1;
  (*pnts_face_el)->nexts[6].nexts[3].int_data = (int *)calloc((5),
      sizeof(int));
  (*pnts_face_el)->nexts[6].nexts[3].int_data[1] = 1; //the local index of point
  (*pnts_face_el)->nexts[6].nexts[3].int_data[2] = 4; //the local index of point
  (*pnts_face_el)->nexts[6].nexts[3].int_data[3] = 5; //the local index of point
  (*pnts_face_el)->nexts[6].nexts[3].int_data[4] = 8; //the local index of point

  //there are four points on the face 4
  (*pnts_face_el)->nexts[6].nexts[4].length = 4+1;
  (*pnts_face_el)->nexts[6].nexts[4].int_data = (int *)calloc((5),
      sizeof(int));
  (*pnts_face_el)->nexts[6].nexts[4].int_data[1] = 2; //the local index of point
  (*pnts_face_el)->nexts[6].nexts[4].int_data[2] = 3; //the local index of point
  (*pnts_face_el)->nexts[6].nexts[4].int_data[3] = 6; //the local index of point
  (*pnts_face_el)->nexts[6].nexts[4].int_data[4] = 7; //the local index of point

  //there are four points on the face 5
  (*pnts_face_el)->nexts[6].nexts[5].length = 4+1;
  (*pnts_face_el)->nexts[6].nexts[5].int_data = (int *)calloc((5),
      sizeof(int));
  (*pnts_face_el)->nexts[6].nexts[5].int_data[1] = 3; //the local index of point
  (*pnts_face_el)->nexts[6].nexts[5].int_data[2] = 4; //the local index of point
  (*pnts_face_el)->nexts[6].nexts[5].int_data[3] = 7; //the local index of point
  (*pnts_face_el)->nexts[6].nexts[5].int_data[4] = 8; //the local index of point

  //there are four points on the face 6
  (*pnts_face_el)->nexts[6].nexts[6].length = 4+1;
  (*pnts_face_el)->nexts[6].nexts[6].int_data = (int *)calloc((5),
      sizeof(int));
  (*pnts_face_el)->nexts[6].nexts[6].int_data[1] = 5; //the local index of point
  (*pnts_face_el)->nexts[6].nexts[6].int_data[2] = 6; //the local index of point
  (*pnts_face_el)->nexts[6].nexts[6].int_data[3] = 7; //the local index of point
  (*pnts_face_el)->nexts[6].nexts[6].int_data[4] = 8; //the local index of point

//completed successfully!
return 0;
}
