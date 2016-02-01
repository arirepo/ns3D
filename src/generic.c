#include "generic.h"

int ElemTypeToMnode(ElemType e)
{
  static int mnode_[NUMBER_OF_ELEM_TYPES] = {2,3,4,4,5,6,8};
  return mnode_[e];
}

int IsInList(const int *list, const int nlist, const int item, int *loc)
{
  int i;

  if (list == NULL) return(0);

  for (i = 0; i < nlist; i++)
    {
      if (list[i] == item)
      {
        if (loc != NULL) *loc = i;
        return (1);
      }
    }

  return(0);
}

int DetectNumberOfSurfaceTags(const int nelem[],
                              const int* const ifactag[],
                              int **taglist, int *maxtag_)
{
  int i;
  int ntag = 0;
  int itag;
  int etype,e;
  int maxtag = -1;
  int loc = 0;

  if (taglist != NULL) free(*taglist);

  for (etype = Triangle; etype <= Quad; etype++)
    {
      for(i = 1; i <= nelem[etype]; i++)
        {
          itag = ifactag[etype][i];

          if (!IsInList(*taglist,ntag,itag,&loc))
            {
              if (taglist != NULL)
                {
                  *taglist = (int *) realloc(*taglist, (ntag + 1)*sizeof(int));
                  //trealloc(taglist,ntag+1);
                  (*taglist)[ntag] = itag;
                  ntag++;
                }

	      maxtag = maxtag > itag ? maxtag:itag;

            }
        }
    }

  if (maxtag_) *maxtag_ = maxtag;

  return(ntag);

}

void TranslateElementWinding(int*nodes,
				    int translation[NUMBER_OF_ELEM_TYPES][MAX_MNODE],
				    ElemType etype,int to_other_format)
{
  int tempnode[MAX_MNODE];
  int mnode = ElemTypeToMnode(etype);

  int k;

  memcpy(tempnode,nodes,mnode*sizeof(int));

  if (to_other_format) // convert our format to other format 
    {
      for (k = 0; k < mnode; k++)
	{
	  nodes[k] = tempnode[translation[etype][k]];
	}
    }
  else // convert other format to our format
    {
      for (k = 0; k < mnode; k++)
	{
	  // translation[k] gives our local node number
	  // k = plot3d local node number

	  nodes[translation[etype][k]] = tempnode[k];
	}
    }

  return;
}

