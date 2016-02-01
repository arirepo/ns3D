#include "generic.h"

static int translation[NUMBER_OF_ELEM_TYPES][MAX_MNODE] = {
  {0,1},   // Edge (not used in fieldview)
  {0,1,2}, // Triangle
  {0,1,2,3}, // Quad
  {0,1,3,2},  // Tet
  //{0,1,2,3},  // Tet (NASA Langley)
  {0,1,2,3,4},  // Pyramid
  {0,3,4,1,5,2},  // Prism
  {0,1,3,2,4,5,7,6}  // Hex
};

// currently writes a FieldView 2.5 file

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
			 )
{
  int ibuf[100];
  char strbuf[80];
  float fbuf[12];
  int itag,ntag,maxtag;
  int i,j,k,kk;
  float *tempxyz = NULL;
  float *tmpvar  = NULL;
  int ntmpfaces;
  int *tmpfaces = NULL;
  int majorversion = 2;
  int minorversion = 5;
  int *taglist = NULL;
  int mnode;
  int code,inode;
  int *tmpelem = NULL;
  int flowtype = flowtype_;
  int did_open = 0;
  int gface, nfaces;
  ElemType etype;
  char const ** nm = NULL;

  switch(flowtype)
    {
    case UNSIO_FLOWTYPE_INCOMPRESSIBLE:
      if (nvariables != 4) flowtype = UNSIO_FLOWTYPE_UNKNOWN;
      break;

    case UNSIO_FLOWTYPE_COMPRESSIBLE:
       if (nvariables != 5) flowtype = UNSIO_FLOWTYPE_UNKNOWN;
       break;

    case UNSIO_FLOWTYPE_VARMACHRUP:
      if (nvariables != 6 || nvariables != 5)
	flowtype = UNSIO_FLOWTYPE_UNKNOWN;
      break;
    }

  if (fp == NULL)
    {
      fp = fopen(filename,"w");
      did_open = 1;
    }

  printf("FIELDVIEW IO: WRITE\n");
  printf("FIELDVIEW IO: file      = %s\n",filename);
  printf("FIELDVIEW IO: version   = %d.%d\n",majorversion,minorversion);

  if (fp == NULL)
    {
      printf("FIELDVIEW IO: Error: unable to open file.\n");
      return(1);
    }

  // write the magic number

  ibuf[0] = 66051;
  fwrite(ibuf,sizeof(int),1,fp);

  // write the FIELDVIEW moniker

  memset(strbuf,0,80);
  strcpy(strbuf,"FIELDVIEW");
  fwrite(strbuf,sizeof(char),80,fp);
 
  // write the FIELDVIEW file version

  fwrite(&majorversion,sizeof(int),1,fp);
  fwrite(&minorversion,sizeof(int),1,fp);

  // write the flow header information

  if (flp != NULL)
    {
      fbuf[0] = (float)flp->elapsed_time;
      fbuf[1] = (float)flp->mach_number;
      fbuf[2] = (float)flp->angle_of_attack;
      fbuf[3] = (float)flp->reynolds_number;
    }
  else
    {
      fbuf[0] = fbuf[1] = fbuf[2] = fbuf[3] = 0.0;
    }

  fwrite(fbuf,sizeof(float),4,fp);
 
  // write the number of grids

  ibuf[0] = 1;
  fwrite(ibuf,sizeof(int),1,fp);

  // detect the number of surface tags

  ntag = DetectNumberOfSurfaceTags(nelem,factag,&taglist,&maxtag);

  printf("FIELDVIEW IO: nsurfs_actual = %d\n",ntag);
  printf("FIELDVIEW IO: nsurfs        = %d\n",maxtag);

  // write the number of boundary types ("surface tags")

  fwrite(&maxtag,sizeof(int),1,fp);

  if (nID != 0 && nID != ntag)
    {
      printf("Warning: number of surface ID's does not match actual number in grid!\n");
      printf("  ntag_grid = %d\n",ntag);
      printf("  ntag      = %d\n",nID);
    }
  
  nm = (const char **) calloc(maxtag, sizeof(char *));

  for (i = 0; i < ntag; i++)
    {
      itag = taglist[i];
      for (k = 0; k < nID; k++) 
	{
	  if (itag == IDs[k] && itag <= maxtag) nm[itag-1] = names[k];
	}
    }

  // write the names of each boundary type ("surface")

  // we make the assumption that if a boundary is being written out, then
  // we are going to be writing boundary data associated with it, i.e., if
  // nbvariables is not equal to 0, the result flag is 1 for all surfaces

  for (i = 0; i < maxtag; i++)
    {
      memset(strbuf,0,80);
      if (nm[i])
	{
	  strcpy(strbuf,nm[i]);
	}
      else
	{
	  sprintf(strbuf,"Surface%d",i+1);
	}

      if (nbvariables > 0)
	{
	  ibuf[0] = 1;
	}
      else
	{
	  ibuf[0] = 0;
	}
      ibuf[1] = 0;

      if (minorversion >= 5)
	{
	  fwrite(&ibuf[0],sizeof(int),1,fp); // result flag
	  fwrite(&ibuf[1],sizeof(int),1,fp); // normal flag
	}

      fwrite(strbuf,sizeof(char),80,fp);

      printf("FIELDVIEW IO: surface   = ID=%d|NAME=\"%s\"\n",i+1,strbuf);
    }

   free(nm);

   // write the number of variables in the file; construct variable
   // names and write those also.

   ibuf[0] = nvariables;

   fwrite(ibuf,sizeof(int),1,fp);  
   printf("FIELDVIEW IO: variables = %d\n",ibuf[0]);
   switch(flowtype)
     {
     case UNSIO_FLOWTYPE_INCOMPRESSIBLE:
       _ASSERT(nvariables == 4);
       strcpy(strbuf,"Pressure");
       fwrite(strbuf,sizeof(char),80,fp);
       printf("FIELDVIEW IO:\t%s\n",strbuf);
       strcpy(strbuf,"XVelocity; Velocity");
       fwrite(strbuf,sizeof(char),80,fp);
       printf("FIELDVIEW IO:\t%s\n",strbuf);
       strcpy(strbuf,"YVelocity");
       fwrite(strbuf,sizeof(char),80,fp);
       printf("FIELDVIEW IO:\t%s\n",strbuf);
       strcpy(strbuf,"ZVelocity");
       fwrite(strbuf,sizeof(char),80,fp);
       printf("FIELDVIEW IO:\t%s\n",strbuf);

       break;

     case UNSIO_FLOWTYPE_COMPRESSIBLE:
       _ASSERT(nvariables == 5);
       strcpy(strbuf,"Density");
       fwrite(strbuf,sizeof(char),80,fp);
       printf("FIELDVIEW IO:\t%s\n",strbuf);
       strcpy(strbuf,"XMomentum; Momentum");
       fwrite(strbuf,sizeof(char),80,fp);
       printf("FIELDVIEW IO:\t%s\n",strbuf);
       strcpy(strbuf,"YMomentum");
       fwrite(strbuf,sizeof(char),80,fp);
       printf("FIELDVIEW IO:\t%s\n",strbuf);
       strcpy(strbuf,"ZMomentum");
       fwrite(strbuf,sizeof(char),80,fp);
       printf("FIELDVIEW IO:\t%s\n",strbuf);
       strcpy(strbuf,"TotalEnergy");
       fwrite(strbuf,sizeof(char),80,fp);
       printf("FIELDVIEW IO:\t%s\n",strbuf);

       break;

     case UNSIO_FLOWTYPE_VARMACHRUP:
       _ASSERT(nvariables == 6 || nvariables == 5);
       strcpy(strbuf,"Density");
       fwrite(strbuf,sizeof(char),80,fp);
       printf("FIELDVIEW IO:\t%s\n",strbuf);
       strcpy(strbuf,"XVelocity; Velocity");
       fwrite(strbuf,sizeof(char),80,fp);
       printf("FIELDVIEW IO:\t%s\n",strbuf);
       strcpy(strbuf,"YVelocity");
       fwrite(strbuf,sizeof(char),80,fp);
       printf("FIELDVIEW IO:\t%s\n",strbuf);
       strcpy(strbuf,"ZVelocity");
       fwrite(strbuf,sizeof(char),80,fp);
       printf("FIELDVIEW IO:\t%s\n",strbuf);
       strcpy(strbuf,"Pressure");
       fwrite(strbuf,sizeof(char),80,fp);
       printf("FIELDVIEW IO:\t%s\n",strbuf);
       if (nvariables == 6)
	 {
	   strcpy(strbuf,"Temperature");
	   fwrite(strbuf,sizeof(char),80,fp);
	   printf("FIELDVIEW IO:\t%s\n",strbuf);
	 }

       break;

     default:
       for (i = 0; i < nvariables; i++)
	 {
	   sprintf(strbuf,"QVariable%d",i+1);
	   fwrite(strbuf,sizeof(char),80,fp);
	 }
     }

   // write the number of boundary variables in the file

   if (minorversion >= 5)
     {
       ibuf[0] = nbvariables;
       fwrite(ibuf,sizeof(int),1,fp);

       printf("FIELDVIEW IO: bvariables = %d\n",ibuf[0]);

       if (nbvariables > 0)
	 {
	   // Create boundary variable names
	   strcpy(strbuf,"Cp");
	   fwrite(strbuf,sizeof(char),80,fp);
	   printf("FIELDVIEW IO:\t%s\n",strbuf);
	   if (nbvariables > 1)
	     {
	       strcpy(strbuf,"Cf");
	       fwrite(strbuf,sizeof(char),80,fp);
	       printf("FIELDVIEW IO:\t%s\n",strbuf);
	     }
	   if (nbvariables > 2)
	     {
	       strcpy(strbuf,"YPlus");
	       fwrite(strbuf,sizeof(char),80,fp);
	       printf("FIELDVIEW IO:\t%s\n",strbuf);
	     }
	   if (nbvariables > 3)
	     {
	       strcpy(strbuf,"wallY");
	       fwrite(strbuf,sizeof(char),80,fp);
	       printf("FIELDVIEW IO:\t%s\n",strbuf);
	     }
	 }
     }
   
  printf("FIELDVIEW IO: nnodes    = %d\n",nnodes);

  // write the FV_NODES tag

  ibuf[0] = 1001;
  fwrite(ibuf,sizeof(int),1,fp);

  // write the number of nodes

  printf("FIELDVIEW IO: writing coordinates...\n");
  fwrite(&nnodes,sizeof(int),1,fp);

  // write the coordinates

  tempxyz = (float *) calloc(nnodes, sizeof(float));

  for (k = 0; k < 3; k++)
    {
      for (i = 1; i <= nnodes; i++) tempxyz[i-1] = (float)xyz[3*i + k];
      fwrite(tempxyz,sizeof(float),nnodes,fp);
    } 
  
  free(tempxyz);

  printf("FIELDVIEW IO: triangles = %d\n",nelem[Triangle]);
  printf("FIELDVIEW IO: quads     = %d\n",nelem[Quad]);

  // write the boundary faces

  for (etype = Triangle; etype <= Quad; etype++)
    {
      mnode = ElemTypeToMnode(etype);

      tmpfaces = (int *) calloc(4*nelem[etype], sizeof(int));

      // pick out the element types individually for the surface tag we 
      // are on.

      for (k = 0; k < maxtag; k++)
	{
	  itag = k+1;

          memset(tmpfaces, 0, 4*nelem[etype]*sizeof(int));
	  //tblank(tmpfaces,4*nelem[etype]);
	  
	  // write FV_FACES tag

	  ibuf[0] = 1002;
	  fwrite(ibuf,sizeof(int),1,fp);

	  // build a list of boundary faces that belong to a particular surface
	  // as well as being of the same type.

	  for (i = 1, ntmpfaces = 0; i <= nelem[etype]; i++)
	    {
	      if (factag[etype][i] == itag)
		{
		  for (kk = 0; kk < mnode; kk++) 
		    tmpfaces[4*ntmpfaces + kk] = c2n[etype][mnode*i + kk];

		  TranslateElementWinding(&tmpfaces[4*ntmpfaces],
					  translation,etype,1);

		  ntmpfaces++;
		}
	    }

	  if (ntmpfaces > 0)
	    {
	      printf("FIELDVIEW IO: writing %d faces for surface ID %d\n",
		     ntmpfaces,itag);
	    }

	  ibuf[0] = itag;
	  ibuf[1] = ntmpfaces;
	  fwrite(ibuf,sizeof(int),2,fp);
	  fwrite(tmpfaces,sizeof(int),4*ntmpfaces,fp);
	}

      free(tmpfaces);
    }

  printf("FIELDVIEW IO: tets      = %d\n",nelem[Tet]);
  printf("FIELDVIEW IO: pyramids  = %d\n",nelem[Pyramid]);
  printf("FIELDVIEW IO: prisms    = %d\n",nelem[Prism]);
  printf("FIELDVIEW IO: hexes     = %d\n",nelem[Hex]);
  printf("FIELDVIEW IO: writing volume elements...\n");

  // write the elements

  ibuf[0] = 1003;
  fwrite(ibuf,sizeof(int),1,fp);
  fwrite(&nelem[Tet],sizeof(int),1,fp);
  fwrite(&nelem[Hex],sizeof(int),1,fp);
  fwrite(&nelem[Prism],sizeof(int),1,fp);
  fwrite(&nelem[Pyramid],sizeof(int),1,fp);

  for (etype = Tet; etype <= Hex; etype++)
    {
      mnode = ElemTypeToMnode(etype);
      switch(etype)
	{
	case Tet:     code = 1; break;
	case Pyramid: code = 2; break;
	case Prism:   code = 3; break;
	case Hex:     code = 4; break;
	default:
	  printf("ERROR: unsupported element type: etype = %d\n",etype);
	  return(1);
	}

      code = code << 18;
      
      tmpelem = (int *) calloc((mnode + 1)*nelem[etype], sizeof(int));

      for (i = 1; i <= nelem[etype]; i++)
	{
	  tmpelem[(mnode+1)*(i-1) + 0] = code;
	  for (kk = 0; kk < mnode; kk++)
	    {
	      inode = c2n[etype][mnode*i + kk];
	      tmpelem[(mnode+1)*(i-1) + 1 + kk] = inode;

	      _ASSERT(inode > 0 && inode <= nnodes);
	    }

	  TranslateElementWinding(&tmpelem[(mnode+1)*(i-1) + 1],
				  translation,etype,1);

	  // debug check
	  for (kk = 0; kk < mnode; kk++)
	    {
	      inode = tmpelem[(mnode+1)*(i-1) + 1 + kk];
	      _ASSERT(inode > 0 && inode <= nnodes);
	    }

	}

      fwrite(tmpelem,(mnode+1)*nelem[etype],sizeof(int),fp);

      free(tmpelem);
    }
      
  // write FV_VARIABLES tag

  ibuf[0] = 1004;
  fwrite(ibuf,sizeof(int),1,fp);

  // write the variables

  tmpvar = (float *) calloc(nnodes, sizeof(float));

  for (i = 0; i < nvariables; i++)
    {
      printf("FIELDVIEW IO: Writing solution variable %d\n",i+1);
      for (k = 1; k <= nnodes; k++)
	{
	  tmpvar[k-1] = (float)qvars[k*nvariables + i];
	}
      fwrite(tmpvar,sizeof(float),nnodes,fp);
    }

  free(tmpvar);

  // Write out boundary variables if they exist.  Here we need to jump through
  // some hoops to get the data into a format that FieldView wants.  The data
  // needs to be sorted by the boundary tag number, and we need to do this
  // for every boundary variable we are dealing with.  Right now, the method
  // is very inefficient, it might be better to create a lookup table and reuse
  // it

  // The header is always written

  if (minorversion >= 5)
    {
      // Header

      ibuf[0] = 1006;
      fwrite(ibuf,sizeof(int),1,fp);

      free(taglist);

    }

  if (did_open) fclose(fp);
  
  printf("FIELDVIEW IO: Done.\n");

  // return success

  return(0);
}
