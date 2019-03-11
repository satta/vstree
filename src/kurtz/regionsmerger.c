#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "types.h"
#include "debugdef.h"
#include "redblackdef.h"
#include "intbits.h"
#include "arraydef.h"
#include "errordef.h"

#include "rbtree.h"

static int compareregions(const void* key,
                           const void* treeelem, /*@unused@*/ void *info)
{
  if(((const PairUint *) key)->uint0 < ((const PairUint *) treeelem)->uint0)
  {
    return (Sint) -1;
  }
  if(((const PairUint *) key)->uint0 > ((const PairUint *) treeelem)->uint0)
  {
    return (Sint) 1;
  }
  return 0;
}

static Sint checkregioncontainment(RBTree *regiontree,PairUint *region)
{
  PairUint *previouselement;

  NOTSUPPOSEDTOBENULL(regiontree);
  previouselement
    = (PairUint *) rbtree_previous_equal_key(regiontree, (const Keytype) region,
                                             compareregions, NULL);
  NOTSUPPOSEDTOBENULL(previouselement);
  if(region->uint0 < previouselement->uint0 ||
     region->uint1 > previouselement->uint1)
  {
    ERROR2("region (%lu,%lu) does not occur in table",
           (Showuint) region->uint0, (Showuint) region->uint1);
    return (Sint) -1;
  }
  return 0;
}

static int collectregions (const Keytype key,
                            RBTreeContext which,
                            /*@unused@ */ Uint depth,
                            void *info)
{
  if (which == RBTREE_POSTORDER || which == RBTREE_LEAF)
  {
    ArrayPairUint *mergedregions = (ArrayPairUint *) info;
    PairUint *region = (PairUint *) key;

    CHECKARRAYSPACE (mergedregions, PairUint, 128);
    mergedregions->spacePairUint[mergedregions->nextfreePairUint++]
      = *region;
  }
  return 0;
}

static Sint checkforoverlappingregions (ArrayPairUint *regions)
{
  PairUint *pptr;

  pptr = regions->spacePairUint;
  DEBUG2(3,"(%lu,%lu)\n",(Showuint) pptr->uint0, (Showuint) pptr->uint1);
  for (pptr=regions->spacePairUint+1;
       pptr < regions->spacePairUint + 
              regions->nextfreePairUint; 
       pptr++)
  {
    DEBUG2(2,"(%lu,%lu)\n",
            (Showuint) pptr->uint0,
            (Showuint) pptr->uint1);
    if ((pptr-1)->uint1 + 1 >= pptr->uint0)
    {
      ERROR4("(%lu,%lu) overlaps with (%lu,%lu)",
             (Showuint) (pptr-1)->uint0,
             (Showuint) (pptr-1)->uint1,
             (Showuint) pptr->uint0,
             (Showuint) pptr->uint1);
      return (Sint) -1;
    }
  }
  printf("no overlaps for %lu merged regions detected\n",
           (Showuint) regions->nextfreePairUint);
  return 0;
}

static Sint checkallregionsforcontainment(RBTree *regiontree,
                                          ArrayPairUint *regions)
{
  Uint i;

  for (i = 0; i < regions->nextfreePairUint; i++)
  {
    DEBUG3 (2, "%lu: check region (%lu,%lu)\n", 
            (Showuint) i,
            (Showuint) regions->spacePairUint[i].uint0,
            (Showuint) regions->spacePairUint[i].uint1);
    if (checkregioncontainment (regiontree, 
                                &regions->spacePairUint[i]) != 0)
    {
      return (Sint) -1;
    }
  }
  return 0;
}

static void markbitsintable(Uint *marktable,Uint left,Uint right)
{
  Uint i;
  for (i=left; i<=right; i++)
  {
    SETIBIT(marktable,i);
  }
}

static void markallregions(Uint *marktable,ArrayPairUint *regions)
{
  PairUint *pptr;

  for (pptr=regions->spacePairUint;
       pptr < regions->spacePairUint + regions->nextfreePairUint; 
       pptr++)
  {
    markbitsintable(marktable,pptr->uint0,pptr->uint1);
  }
}

static void comparemarkedbits(Uint *marktable1,Uint *marktable2,Uint width)
{
  Uint i;
  BOOL ismarked1, ismarked2;

  for(i=0; i<width; i++)
  {
    ismarked1 = ISIBITSET(marktable1,i) ? True : False;
    ismarked2 = ISIBITSET(marktable2,i) ? True : False;
    if(ismarked1 != ismarked2)
    {
      fprintf(stderr,"index %lu: ismarked1=%s != %s=ismarked2\n",
           (Showuint) i,SHOWBOOL(ismarked1),SHOWBOOL(ismarked2));
      exit(EXIT_FAILURE);
    }
  }
  printf("tables of width %lu are identical\n",(Showuint) width);
}

Sint verifyregiontree(RBTree *regiontree,Uint width,
                      ArrayPairUint *originalregions)
{
  ArrayPairUint mergedregions;
  Uint *markmerged, *markorig;

  INITARRAY (&mergedregions, PairUint);
  if (rbtree_walk (regiontree,
                   collectregions, &mergedregions) != 0)
  {
    return (Sint) -1;
  }
  if(checkforoverlappingregions (&mergedregions) != 0)
  {
    return (Sint) -2;
  }
  if(checkallregionsforcontainment(regiontree,originalregions) != 0)
  {
    return (Sint) -3;
  }
  INITBITTAB(markmerged,width);
  INITBITTAB(markorig,width);
  markallregions(markmerged,&mergedregions);
  markallregions(markorig,originalregions);
  comparemarkedbits(markmerged,markorig,width);
  FREESPACE(markmerged);
  FREESPACE(markorig);
  FREEARRAY (&mergedregions, PairUint);
  return 0;
}

static int mergewithonepreviouselement(RBTree *regiontree, 
                                        PairUint **node2storeregion,
                                        PairUint *region)
{
  PairUint *previouselement;

  if(regiontree == NULL)
  {
    return 0;
  }
  previouselement
    = (PairUint *) rbtree_previous_equal_key(regiontree, (const Keytype) region,
                                             compareregions,
                                             NULL);
  if(previouselement != NULL && previouselement->uint1 + 1 >= region->uint0)
  {
    DEBUG4(2,"left merge of previous (%lu,%lu) with (%lu,%lu)\n",
           (Showuint) previouselement->uint0,
           (Showuint) previouselement->uint1,
           (Showuint) region->uint0,
           (Showuint) region->uint1);
    region->uint0 = previouselement->uint0;
    if(previouselement->uint1 > region->uint1)
    {
      region->uint1 = previouselement->uint1;
    }
    *node2storeregion = previouselement;
  } 
  return 0;
}

static int mergewithallrightelements(RBTree *regiontree, 
                                      PairUint **node2storeregion,
                                      PairUint *region)
{
  PairUint *nextelement;

  while(True)
  {
    nextelement = (PairUint *) rbtree_next_key(regiontree, (const Keytype) region,
                                               compareregions, NULL);
   
    if(nextelement != NULL && region->uint1 + 1 >= nextelement->uint0)
    {
      DEBUG2(2,"right merge with (%lu,%lu)\n",
             (Showuint) nextelement->uint0,(Showuint) nextelement->uint1);
      if(region->uint1 <= nextelement->uint1)
      {
        region->uint1 = nextelement->uint1;
        if(*node2storeregion == NULL)
        {
          *node2storeregion = nextelement;
        } else
        {
          if(rbtree_erase(regiontree, (const Keytype) nextelement) != 0)
          {
            ERROR2("cannot delete (%lu,%lu)",
                   (Showuint) nextelement->uint0,
                   (Showuint) nextelement->uint1);
            return (Sint) -1;
          }
          free(nextelement);
        }
        break;
      }
      
      if(rbtree_erase(regiontree, (const Keytype) nextelement) != 0)
      {
        ERROR2("cannot delete (%lu,%lu)",
               (Showuint) nextelement->uint0,
               (Showuint) nextelement->uint1);
        return (Sint) -1;
      }
      free(nextelement);
    } else
    {
      break;
    }
  }
  return 0;
}

static int dotheregioninsertion(RBTree *regiontree, 
                                 BOOL makedatacopy,
                                 bool *nodecreated, 
                                 PairUint *region)
{
  if(makedatacopy)
  {
    PairUint *copyofregion;
    copyofregion = malloc(sizeof(PairUint));
    if(copyofregion == NULL)
    {
      ERROR0("cannot allocate space for copyofregion");
      return (Sint) -1;
    }
    *copyofregion = *region;
    (void) rbtree_search(regiontree, (const Keytype) copyofregion,
                              nodecreated);
  } else
  {
    (void) rbtree_search(regiontree, (const Keytype) region,
                              nodecreated);
  }
  return 0;
}


Sint insertnewregion(RBTree **regiontree, 
                     BOOL makedatacopy,
                     bool *nodecreated, 
                     PairUint *region)
{
  if(*regiontree == NULL)
  {
    *regiontree = rbtree_new(compareregions, NULL, NULL);
    if(dotheregioninsertion(*regiontree, 
                            makedatacopy,
                            nodecreated, 
                            region) != 0)
    {
      return (Sint) -1;
    }
  } else
  {
    PairUint *node2storeregion = NULL;

    if(mergewithonepreviouselement(*regiontree, 
                                   &node2storeregion, region) != 0)
    {
      return (Sint) -2;
    }
    if(mergewithallrightelements(*regiontree, 
                                 &node2storeregion, region) != 0)
    {
      return (Sint) -3;
    }
    if(node2storeregion == NULL)
    {
      DEBUG2(2, "insert (%lu,%lu)\n",
                 (Showuint) region->uint0, (Showuint) region->uint1);
      if(dotheregioninsertion(*regiontree, 
                              makedatacopy,
                              nodecreated, 
                              region) != 0)
      {
        return (Sint) -4;
      }
      if(!(*nodecreated))
      {
        ERROR2("no node for insert (%lu,%lu) created\n",
                 (Showuint) region->uint0, (Showuint) region->uint1);
        return (Sint) -5;
      }
    } else
    {
      *nodecreated = False;
      *node2storeregion = *region;
    }
  }
  return 0;
}
