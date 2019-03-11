#include "types.h"
#include "debugdef.h"
#include "arraydef.h"
#include "redblackdef.h"

#include "rbtree.h"

typedef struct hnode_t
{
  Uint symbol, count;
  struct hnode_t *leftchild, 
                 *rightchild;
} Huffmannode;

typedef struct
{
  Huffmannode *roothuffmantree;
  Uint insertednodes, largestsymbol,
       totalbits, totalnumofvalues;
  RBTree *rbhuffroot;
} Huffmaninfo;

static int cmpHuffmannode(const void *keya,
                          const void *keyb,
                           /*@unused@*/ void *info)
{
  const Huffmannode *h1 = (const Huffmannode *) keya,
                    *h2 = (const Huffmannode *) keyb;
  
  if(h1->count < h2->count)
  {
     return (Sint) -1;
  }
  if(h1->count > h2->count)
  {
    return (Sint) 1;
  }
  if(h1->symbol < h2->symbol)
  {
    return (Sint) -1;
  }
  if(h1->symbol > h2->symbol)
  {
    return (Sint) 1;
  }
  return 0;
}

static void initialhuffmaninsert(Huffmaninfo *huffmaninfo,
                                 const ArrayUint *distribution)
{
  Huffmannode *huffptr;
  bool nodecreated;
  Uint i; 

  huffmaninfo->rbhuffroot = rbtree_new(cmpHuffmannode, NULL, NULL);
  huffmaninfo->insertednodes = 0;
  for(i = 0; i < distribution->nextfreeUint; i++)
  {
    if(distribution->spaceUint[i] > 0)
    {
      huffptr = malloc(sizeof(Huffmannode));
      if(huffptr == NULL)
      {
        fprintf(stderr,"cannot allocate space for initial Huffmannode\n");
        exit(EXIT_FAILURE);
      }
      huffptr->count = distribution->spaceUint[i];
      huffptr->symbol = i;
      huffptr->leftchild = NULL;
      huffptr->rightchild = NULL;
      (void) rbtree_search(huffmaninfo->rbhuffroot,
		           huffptr,
                           &nodecreated);
      huffmaninfo->insertednodes++;
      huffmaninfo->largestsymbol = i;
    }
  }
}

static Sint makehuffmantree(Huffmaninfo *huffmaninfo)
{
  Huffmannode *n1, *n2, *newnode = NULL;
  Uint i, currentsymbol = huffmaninfo->largestsymbol + 1;
  bool nodecreated;

  if(huffmaninfo->insertednodes == 0)
  {
    huffmaninfo->roothuffmantree = NULL;
    return (Sint) 0;
  }
  if(huffmaninfo->insertednodes == UintConst(1))
  {
    huffmaninfo->roothuffmantree
      = (Huffmannode *) rbtree_root_key(huffmaninfo->rbhuffroot);
  } else
  {
    for(i=0; i<huffmaninfo->insertednodes-1; i++)
    {
      n1 = rbtree_minimum_key(huffmaninfo->rbhuffroot);
      NOTSUPPOSEDTOBENULL(n1);
      if(rbtree_erase(huffmaninfo->rbhuffroot, n1) != 0)
      {
        return (Sint) -1;
      }
      n2 = rbtree_minimum_key(huffmaninfo->rbhuffroot);
      NOTSUPPOSEDTOBENULL(n2);
      if(rbtree_erase(huffmaninfo->rbhuffroot, n2) != 0)
      {
        return (Sint) -1;
      }
      newnode = malloc(sizeof(Huffmannode));
      if(newnode == NULL)
      {
        fprintf(stderr,"cannot allocate space for Huffmannode in trees\n");
        exit(EXIT_FAILURE);
      }
      newnode->count = n1->count + n2->count;
      if(n1->count < n2->count)
      {
        newnode->leftchild = n2;
        newnode->rightchild = n1;
      } else
      {
        newnode->leftchild = n1;
        newnode->rightchild = n2;
      }
      newnode->symbol = currentsymbol++;
      (void) rbtree_search(huffmaninfo->rbhuffroot, newnode,
                           &nodecreated);
    }
    huffmaninfo->roothuffmantree = (Huffmannode *) newnode;
  }
  return 0;
}

#ifdef DEBUG
static void showsinglehuffmancode(Uint length,Uint code)
{
  if(length > 0)
  {
    Uint leftbit;

    for(leftbit = UintConst(1) << (length-1);
        leftbit != 0;
        leftbit >>= 1)
    {
      (void) putchar((code & leftbit) ? '1' : '0');
    }
  }
}
#endif

static void recurseextractallhuffmancodes(Huffmaninfo *huffmaninfo,
                                          Uint length,Uint code,
                                          Huffmannode *hnode)
{
  if(hnode->leftchild == NULL)
  {
#ifdef DEBUG
    printf("symbol %lu, count %lu, codelength %lu: ",
             (Showuint) hnode->symbol,
             (Showuint) hnode->count,
             (Showuint) length);
    showsinglehuffmancode(length,code);
    printf("\n");
#endif
    huffmaninfo->totalbits += length * hnode->count;
    huffmaninfo->totalnumofvalues += hnode->count;
  } else
  {
    recurseextractallhuffmancodes(huffmaninfo,length+1,
                                  code << 1,hnode->leftchild);
    recurseextractallhuffmancodes(huffmaninfo,length+1,
                                  (code << 1) | UintConst(1),
                                  hnode->rightchild);
  }
}

static void extractallhuffmancodes(Huffmaninfo *huffmaninfo)
{
  if(huffmaninfo->roothuffmantree != NULL)
  {
    recurseextractallhuffmancodes(huffmaninfo,0,0,huffmaninfo->roothuffmantree);
  }
}

Sint huffmanencoding(Uint *totalbits,
                     Uint *totalnumofvalues,
                     const ArrayUint *distribution)
{
  Huffmaninfo huffmaninfo; 

  initialhuffmaninsert(&huffmaninfo,distribution);
  if(makehuffmantree(&huffmaninfo) != 0)
  {
    return (Sint) -1;
  }
  huffmaninfo.totalbits = 0;
  huffmaninfo.totalnumofvalues = 0;
  extractallhuffmancodes(&huffmaninfo);
  if(huffmaninfo.rbhuffroot != NULL)
  {
    rbtree_delete(huffmaninfo.rbhuffroot);
    free(huffmaninfo.rbhuffroot);
  }
  huffmaninfo.rbhuffroot = NULL;
  *totalbits = huffmaninfo.totalbits;
  *totalnumofvalues = huffmaninfo.totalnumofvalues;
  return 0;
}
