#include "generic.h"

//
// table converts from ugrid --> our format
//
static int translation[NUMBER_OF_ELEM_TYPES][MAX_MNODE] = {
  {0,1},   // Edge
  {2,1,0}, // Triangle
  {3,2,1,0}, // Quad
  {0,1,2,3},  // Tet
  {0,3,4,1,2},  // Pyramid
  {0,1,2,3,4,5},  // Prism
  {0,1,2,3,4,5,6,7}  // Hex
};


int Read_asciiugrid_Grid(FILE *fp,
			 int *nnodes,
			 int nelem[],
			 double**xyz,
			 int*c2n[],
			 int*factag[]
			 )
{
  int i;
  int k,mnode;
  ElemType etype;

  memset(nelem, 0, NUMBER_OF_ELEM_TYPES*sizeof(int));
  //tblank(nelem,NUMBER_OF_ELEM_TYPES);

  //
  // read the header
  //
  fscanf(fp,"%d %d %d %d %d %d %d\n",
	 nnodes,
	 &nelem[Triangle],
	 &nelem[Quad],
	 &nelem[Tet],
	 &nelem[Pyramid],
	 &nelem[Prism],
	 &nelem[Hex]);

  printf("PLOT3D IO: nodes     = %d\n",*nnodes);
  printf("PLOT3D IO: triangles = %d\n",nelem[Triangle]);
  printf("PLOT3D IO: quads     = %d\n",nelem[Quad]);
  printf("PLOT3D IO: tets      = %d\n",nelem[Tet]);
  printf("PLOT3D IO: pyramids  = %d\n",nelem[Pyramid]);
  printf("PLOT3D IO: prisms    = %d\n",nelem[Prism]);
  printf("PLOT3D IO: hexes     = %d\n",nelem[Hex]);

  //
  // read the coordinates
  //
  printf("PLOT3D IO: reading coordinates...\n");

  *xyz = (double *) calloc(3*(*nnodes + 1), sizeof(double));

  for (i = 1; i <= *nnodes; i++)
    {
      fscanf(fp,"%lf %lf %lf",
	     &(*xyz)[3*i + 0],
	     &(*xyz)[3*i + 1],
	     &(*xyz)[3*i + 2]);
    }

  //
  // read the boundary faces
  //
  printf("PLOT3D IO: reading surface elements...\n");
  for (etype = Triangle; etype <= Quad; etype++) 
    {
      mnode = ElemTypeToMnode(etype);

      c2n[etype] = (int *) calloc((nelem[etype] + 1)*mnode, sizeof(int));

      for (i = 1; i <= nelem[etype]; i++)
	{
	  for (k = 0; k < mnode; k++)
	    {
	      fscanf(fp,"%d",&c2n[etype][mnode*i + k]);
	    }
	}

      for (i = 1; i <= nelem[etype]; i++)
	{
	  TranslateElementWinding(c2n[etype] + mnode*i,translation,etype,0);
	}
    }



  //
  // read the surface tags
  //
   for (etype = Triangle; etype <= Quad; etype++)
    {
      factag[etype] = (int *) calloc(nelem[etype] + 1, sizeof(int));

      for (i = 1; i <= nelem[etype]; i++)
	fscanf(fp,"%d",&factag[etype][i]);
    }


  //
  // read the elements
  //
  printf("PLOT3D IO: reading volume elements...\n");
  for (etype = Tet; etype <= Hex; etype++)
    {
      mnode = ElemTypeToMnode(etype);

      c2n[etype] = (int *) calloc((nelem[etype] + 1)*mnode, sizeof(int));

      for (i = 1; i <= nelem[etype]; i++)
	for (k = 0; k < mnode; k++)
	  fscanf(fp,"%d",&c2n[etype][mnode*i + k]);

      if (etype == Pyramid)
	{
	  for (i = 1; i <= nelem[etype]; i++)
	    {
	      TranslateElementWinding(c2n[etype] + mnode*i, translation, etype,0);
	    }
	}
    }

  printf("PLOT3D IO: Done.\n");


  //
  // return success
  //
  return(0);
}

