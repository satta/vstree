/*
  Copyright (c)      2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008-2012 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/
/* RBTree.c - red-black tree
 * Copyright (c) 2009 Charles Wilson
 * This file is part of the C Data Structures Library.
 *
 * Derived from code released into the public domain by Julienne Walker:
 * http://eternallyconfuzzled.com/tuts/datastructures/jsw_tut_rbtree.aspx
 *   > Created (Julienne Walker): August 23, 2003
 *   > Modified (Julienne Walker): March 14, 2008
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
*/

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "rbtree.h"

#ifndef HEIGHT_LIMIT
#define HEIGHT_LIMIT 64
#endif

typedef struct RBTreeNode
{
  int red;
  void *key;
  struct RBTreeNode *link[2];  /* Left (0) and right (1) links */
} RBTreeNode;

struct RBTree
{
  RBTreeNode *root;
  CompareWithData cmp;
  RBTreeFreeFunc free;
  size_t size;
  void *info;
};

struct RBTreeIter
{
  const RBTree *tree;
  RBTreeNode *it;
  RBTreeNode *path[HEIGHT_LIMIT];
  size_t top;
};

#define RBTREE_NODE_IS_RED(node)\
        (node != NULL && node->red == 1)

static inline RBTreeNode *rbtree_single(RBTreeNode *root, int dir)
{
  RBTreeNode *save;
  save = root->link[(int) !dir];
  root->link[(int) !dir] = save->link[dir];
  save->link[dir] = root;

  root->red = 1;
  save->red = 0;

  return save;
}

static inline RBTreeNode *rbtree_double(RBTreeNode *root, int dir)
{
  root->link[(int) !dir] = rbtree_single(root->link[(int) !dir], (int) !dir);
  return rbtree_single (root, dir);
}

static RBTreeNode *rbtree_new_node(void *key)
{
  RBTreeNode *rn = (RBTreeNode *) malloc(sizeof *rn);

  rn->red = 1;
  rn->key = key;
  rn->link[0] = rn->link[1] = NULL;

  return rn;
}

RBTree *rbtree_new(CompareWithData cmp, RBTreeFreeFunc free,
                        void *info)
{
  RBTree *rt = (RBTree *) malloc(sizeof *rt);
  assert(cmp);

  rt->root = NULL;
  rt->cmp = cmp;
  rt->free = free;
  rt->info = info;
  rt->size = 0;

  return rt;
}

void rbtree_delete(RBTree *tree)
{
  rbtree_clear(tree);
  free(tree);
}

void rbtree_clear(RBTree *tree)
{
  RBTreeNode *it;
  RBTreeNode *save;
  if (tree != NULL) {
    it = tree->root;
    if (it != NULL) {

      /* Rotate away the left links so that we can treat this like the
         destruction of a linked list */
      while (it != NULL) {
        if (it->link[0] == NULL) {
          /* No left links, just kill the node and move on */
          save = it->link[1];
          if (tree->free != NULL)
            tree->free(it->key);
          free(it);
        } else {
          /* Rotate away the left link and check again */
          save = it->link[0];
          it->link[0] = save->link[1];
          save->link[1] = it;
        }
        it = save;
      }
      tree->size = 0;
      tree->root = NULL;
    }
  }
}

static inline void *rbtree_find_with_cmp_g(const RBTree *tree, void *key,
                                              CompareWithData cmpfunc,
                                              void *info)
{
  RBTreeNode *it = tree->root;

  while (it != NULL) {
    int cmp = cmpfunc(it->key, key, info);

    if (cmp == 0)
      break;

    it = it->link[(int) (cmp < 0)];
  }

  return it == NULL ? NULL : it->key;
}

void *rbtree_find_with_cmp(const RBTree *tree, void *key,
                              CompareWithData cmpfunc, void *info)
{
  assert(tree);
  assert(cmpfunc);
  assert(key);
  return rbtree_find_with_cmp_g(tree, key, cmpfunc, info);
}

void *rbtree_find(const RBTree *tree, void *key)
{
  assert(tree);
  assert(key);
  return rbtree_find_with_cmp_g(tree, key, tree->cmp, tree->info);
}

static inline void rbtree_insert_g(RBTree *tree, void *key,
                                      bool *nodecreated,
                                      CompareWithData cmpfunc,
                                      void *info)
{
  *nodecreated = false;
  if (tree->root == NULL) {
    /* We have an empty tree; attach the new node directly to the root */
    tree->root = rbtree_new_node(key);
    *nodecreated = true;

  } else {
    RBTreeNode head = {0,0,{NULL, NULL}}; /* False tree root */
    RBTreeNode *g, *t;     /* Grandparent & parent */
    RBTreeNode *p, *q;     /* Iterator & parent */
    int dir = 0, last = 0;

    /* Set up our helpers */
    t = &head;
    g = p = NULL;
    q = t->link[1] = tree->root;

    /* Search down the tree for a place to insert */
    while (true) {
      int cmp;
      if (q == NULL) {
        /* Insert a new node at the first null link */
        assert(p != NULL);
        p->link[dir] = q = rbtree_new_node(key);
        *nodecreated = true;

      } else if (RBTREE_NODE_IS_RED(q->link[0] )
                 && RBTREE_NODE_IS_RED(q->link[1])) {
        /* Simple red violation: color flip */
        q->red = 1;
        q->link[0]->red = 0;
        q->link[1]->red = 0;
      }

      if ( RBTREE_NODE_IS_RED(q) && RBTREE_NODE_IS_RED(p)) {
        /* Hard red violation: rotations necessary */
        int dir2 = (int) (t->link[1] == g);

        if (q == p->link[last])
          t->link[dir2] = rbtree_single(g, (int) !last);
        else
          t->link[dir2] = rbtree_double(g, (int) !last);
      }

      cmp = cmpfunc(q->key, key, info);

      /* Stop working if we inserted a node. This check also disallows
         duplicates in the tree */
      if (cmp == 0)
        break;

      last = dir;
      dir = (int) (cmp < 0);

      /* Move the helpers down */
      if (g != NULL)
        t = g;

      g = p;
      p = q;
      q = q->link[dir];
    }

    /* Update the root (it may be different) */
    tree->root = head.link[1];
  }

  /* Make the root black for simplified logic */
  tree->root->red = 0;
  if (*nodecreated)
    ++tree->size;
}

void rbtree_insert(RBTree *tree, void *key)
{
   bool nodecreated;
  assert(tree);
  assert(key);
  rbtree_insert_g(tree, key, &nodecreated, tree->cmp, tree->info);
}

void rbtree_insert_with_cmp(RBTree *tree, void *key,
                               CompareWithData cmpfunc, void *info)
{
   bool nodecreated;
  assert(tree);
  assert(key);
  assert(cmpfunc);
  rbtree_insert_g(tree, key, &nodecreated, cmpfunc, info);
}

void* rbtree_search(RBTree *tree, void *key, bool *nodecreated)
{
  assert(tree);
  assert(key);
  assert(nodecreated);
  rbtree_insert_g(tree, key, nodecreated, tree->cmp, tree->info);
  return key;
}

void* rbtree_search_with_cmp(RBTree *tree, void *key,
                                CompareWithData cmpfunc, void *info,
                                bool *nodecreated)
{
  assert(tree);
  assert(nodecreated);
  assert(cmpfunc);
  assert(key);
  rbtree_insert_g(tree, key, nodecreated, cmpfunc, info);
  return key;
}

int rbtree_erase(RBTree *tree, void *key)
{
  int rv = -1;
  assert(tree);
  assert(key);
  if (tree->root != NULL) {
    RBTreeNode head = {0,0,{NULL, NULL}}; /* False tree root */
    RBTreeNode *q, *p, *g; /* Helpers */
    RBTreeNode *f = NULL;  /* Found item */
    int dir = 1;

    /* Set up our helpers */
    q = &head;
    p = NULL;
    q->link[1] = tree->root;

    /* Search and push a red node down to fix red violations as we go */
    while (q->link[dir] != NULL) {
      int last = dir;

      /* Move the helpers down */
      g = p, p = q;
      q = q->link[dir];
      dir = (int) (tree->cmp(q->key, key, tree->info) < 0);

      /* Save the node with matching key and keep going; we'll do removal
         tasks at the end */
      if (tree->cmp(q->key, key, tree->info) == 0)
        f = q;

      /* Push the red node down with rotations and color flips */
      if (!RBTREE_NODE_IS_RED(q)
            && !RBTREE_NODE_IS_RED(q->link[dir])) {
        if (RBTREE_NODE_IS_RED(q->link[(int) !dir])) {
          p = p->link[last] = rbtree_single(q, dir);
        } else if (!RBTREE_NODE_IS_RED(q->link[(int) !dir])) {
          RBTreeNode *s = p->link[(int) !last];

          if (s != NULL) {
            if (!RBTREE_NODE_IS_RED(s->link[(int) !last])
                  && !RBTREE_NODE_IS_RED(s->link[last])) {
              /* Color flip */
              p->red = 0;
              s->red = 1;
              q->red = 1;
            } else {
              int dir2;
              assert(g != NULL);
              dir2 = (int) (g->link[1] == p);

              if (RBTREE_NODE_IS_RED(s->link[last]))
                g->link[dir2] = rbtree_double(p, last);
              else if (RBTREE_NODE_IS_RED(s->link[(int) !last]))
                g->link[dir2] = rbtree_single(p, last);

              /* Ensure correct coloring */
              q->red = g->link[dir2]->red = 1;
              g->link[dir2]->link[0]->red = 0;
              g->link[dir2]->link[1]->red = 0;
            }
          }
        }
      }
    }

    /* Replace and remove the saved node */
    if (f != NULL) {
      if (tree->free != NULL)
        tree->free(f->key);
      f->key = q->key;
      assert(p);
      p->link[(int) (p->link[1] == q)] = q->link[(int) (q->link[0] == NULL)];
      free(q);
      rv = 0;
      --tree->size;
    }

    /* Update the root (it may be different) */
    tree->root = head.link[1];

    /* Make the root black for simplified logic */
    if (tree->root != NULL)
      tree->root->red = 0;
  }

  return rv;
}

size_t rbtree_size(const RBTree *tree)
{
  assert(tree);
  return tree->size;
}

static int rbtree_recurse(const RBTreeNode *root, RBTreeAction action,
                             Uint level, void *actinfo)
{
  if (root->link[0] == NULL && root->link[1] == NULL) {
    if (action(root->key, RBTREE_LEAF, level, actinfo) != 0) {
      return -1;
    }
  } else {
    if (action(root->key, RBTREE_PREORDER, level, actinfo) != 0) {
      return -2;
    }
    if (root->link[0] != NULL) {
      if (rbtree_recurse(root->link[0], action, level + 1, actinfo) != 0) {
        return -3;
      }
    }
    if (action(root->key, RBTREE_POSTORDER, level, actinfo) != 0) {
      return -4;
    }
    if (root->link[1] != NULL) {
      if (rbtree_recurse(root->link[1], action, level + 1, actinfo) != 0) {
        return -5;
      }
    }
    if (action(root->key, RBTREE_ENDORDER, level, actinfo) != 0) {
      return -6;
    }
  }
  return 0;
}

#define RBT_CHECK_RETURN_CODE\
        if (retcode < 0 || retcode == 1) {\
          return retcode;\
        }

static int rbtree_recursewithstop(const RBTreeNode *root,
                                     RBTreeAction action,
                                     Uint level, void *actinfo)
{
  int retcode;

  if (root->link[0] == NULL && root->link[1] == NULL) {
    retcode = action(root->key, RBTREE_LEAF, level, actinfo);
    RBT_CHECK_RETURN_CODE;
  } else {
    retcode = action(root->key, RBTREE_PREORDER, level, actinfo);
    RBT_CHECK_RETURN_CODE;
    if (root->link[0] != NULL) {
      retcode = rbtree_recursewithstop(root->link[0], action,
                                      level + 1, actinfo);
      RBT_CHECK_RETURN_CODE;
    }
    retcode = action(root->key, RBTREE_POSTORDER, level, actinfo);
    RBT_CHECK_RETURN_CODE;
    if (root->link[1] != NULL) {
      retcode = rbtree_recursewithstop(root->link[1], action,
                                      level + 1, actinfo);
      RBT_CHECK_RETURN_CODE;
    }
    retcode = action(root->key, RBTREE_ENDORDER, level, actinfo);
    RBT_CHECK_RETURN_CODE;
  }
  return 0;
}

static int rbtree_recursereverseorder(const RBTreeNode *root,
                                         RBTreeAction action,
                                         Uint level, void *actinfo)
{
  if (root->link[0] == NULL && root->link[1] == NULL) {
    if (action(root->key, RBTREE_LEAF, level, actinfo) != 0)
    {
      return -1;
    }
  } else {
    if (action(root->key, RBTREE_PREORDER, level, actinfo) != 0) {
      return -2;
    }
    if (root->link[1] != NULL) {
      if (rbtree_recursereverseorder(root->link[1], action, level + 1,
                                        actinfo) != 0) {
        return -3;
      }
    }
    if (action(root->key, RBTREE_POSTORDER, level, actinfo) != 0) {
      return -4;
    }
    if (root->link[0] != NULL) {
      if (rbtree_recursereverseorder(root->link[0], action, level + 1,
                                        actinfo) != 0) {
        return -5;
      }
    }
    if (action (root->key, RBTREE_ENDORDER, level, actinfo) != 0) {
      return -6;
    }
  }
  return 0;
}

int rbtree_walk(const RBTree *tree, RBTreeAction action, void *actinfo)
{
  assert(tree);
  assert(action);
  if (tree->root != NULL) {
    if (rbtree_recurse(tree->root, action, 0, actinfo) != 0) {
      return -1;
    }
  }
  return 0;
}

int rbtree_walk_stop(const RBTree *tree, RBTreeAction action,
                        void *actinfo)
{
  assert(tree);
  assert(action);
  if (tree->root != NULL) {
    int retcode = rbtree_recursewithstop(tree->root, action, 0, actinfo);
    RBT_CHECK_RETURN_CODE;
  }
  return 0;
}

int rbtree_walk_reverse(const RBTree *tree, RBTreeAction action,
                           void *actinfo)
{
  assert(tree);
  assert(action);
  if (tree->root != NULL) {
    if (rbtree_recursereverseorder(tree->root, action, 0, actinfo) != 0) {
      return -1;
    }
  }
  return 0;
}

static inline void* rbtree_minimum_key_for_node(const RBTreeNode *root)
{
  if (root == NULL) {
    return NULL;
  }
  while (root->link[0] != NULL) {
    root = root->link[0];
  }
  return root->key;
}

static inline void* rbtree_maximum_key_for_node(const RBTreeNode *root)
{
  if (root == NULL) {
    return NULL;
  }
  while (root->link[1] != NULL) {
    root = root->link[1];
  }
  return root->key;
}

void* rbtree_minimum_key(const RBTree *tree)
{
  assert(tree);
  return rbtree_minimum_key_for_node(tree->root);
}

void* rbtree_maximum_key(const RBTree *tree)
{
  assert(tree);
  return rbtree_maximum_key_for_node(tree->root);
}

void* rbtree_root_key(const RBTree *tree)
{
  assert(tree);
  if (tree->size == 0 || tree->root == NULL)
    return NULL;
  else
    return tree->root->key;
}

void* rbtree_previous_key(const RBTree *tree, void *key,
                             CompareWithData cmpfun, void *cmpinfo)
{
  int cmp;
  const RBTreeNode *current,
                     *found = NULL;
  assert(tree);
  assert(key);
  assert(cmpfun);
  current = tree->root;

  while (current != NULL) {
    cmp = cmpfun(key, current->key, cmpinfo);
    if (cmp == 0) {
      if (current->link[0] == NULL) {
        if (found == NULL) {
          return NULL;
        }
        return found->key;
      }
      return rbtree_maximum_key_for_node(current->link[0]);
    } else {
      if (cmp < 0) {
        current = current->link[0];
      } else {
        found = current;
        current = current->link[1];
      }
    }
  }
  if (found == NULL) {
    return NULL;
  }
  return found->key;
}

void* rbtree_previous_equal_key(const RBTree *tree, void *key,
                                   CompareWithData cmpfun, void *cmpinfo)
{
  int cmp;
  const RBTreeNode *current,
                     *found = NULL;
  assert(tree);
  assert(key);
  assert(cmpfun);
  current = tree->root;

  while (current != NULL) {
    cmp = cmpfun(key, current->key, cmpinfo);
    if (cmp == 0) {
      return current->key;
    } else {
      if (cmp < 0) {
        current = current->link[0];
      } else {
        found = current;
        current = current->link[1];
      }
    }
  }
  if (found == NULL) {
    return NULL;
  }
  return found->key;
}

void* rbtree_next_key(const RBTree *tree, void *key,
                         CompareWithData cmpfun, void *cmpinfo)
{
  int cmp;
  const RBTreeNode *current,
                     *found = NULL;
  assert(tree);
  assert(key);
  assert(cmpfun);
  current = tree->root;

  while (current != NULL) {
    cmp = cmpfun(key, current->key, cmpinfo);
    if (cmp == 0) {
      if (current->link[1] == NULL) {
        if (found == NULL) {
          return NULL;
        }
        return found->key;
      }
      return rbtree_minimum_key_for_node(current->link[1]);
    } else {
      if (cmp < 0) {
        found = current;
        current = current->link[0];
      } else {
        current = current->link[1];
      }
    }
  }
  if (found == NULL) {
    return NULL;
  }
  return found->key;
}

void* rbtree_next_equal_key(const RBTree *tree, void *key,
                               CompareWithData cmpfun, void *cmpinfo)
{
  int cmp;
  const RBTreeNode *current,
                     *found = NULL;
  assert(tree);
  assert(key);
  assert(cmpfun);
  current = tree->root;

  while (current != NULL) {
    cmp = cmpfun(key, current->key, cmpinfo);
    if (cmp == 0) {
      return current->key;
    } else {
      if (cmp < 0) {
        found = current;
        current = current->link[0];
      } else {
        current = current->link[1];
      }
    }
  }
  if (found == NULL) {
    return NULL;
  }
  return found->key;
}

static inline void *start(RBTreeIter *trav, const RBTree *tree, int dir)
{
  trav->tree = tree;
  trav->it = tree->root;
  trav->top = 0;

  /* Save the path for later traversal */
  if (trav->it != NULL)
    {
      while (trav->it->link[dir] != NULL)
        {
          trav->path[trav->top++] = trav->it;
          trav->it = trav->it->link[dir];
        }
    }

  return trav->it == NULL ? NULL : trav->it->key;
}

static inline void *move(RBTreeIter *trav, int dir,  bool strict)
{
  if (trav->it->link[dir] != NULL) {
    /* Continue down this branch */
    trav->path[trav->top++] = trav->it;
    trav->it = trav->it->link[dir];

    while (trav->it->link[(int) !dir] != NULL) {
      trav->path[trav->top++] = trav->it;
      trav->it = trav->it->link[(int) !dir];
    }
  } else {
    /* Move to the next branch */
    RBTreeNode *last;

    do {
      if (trav->top == 0) {
        trav->it = NULL;
        break;
      }

      last = trav->it;
      trav->it = trav->path[--trav->top];
    } while (last == trav->it->link[dir]);
  }

  return trav->it == NULL ? NULL : trav->it->key;
}

RBTreeIter* rbtree_iter_new_from_first(const RBTree *tree)
{
  RBTreeIter *trav = malloc(sizeof (RBTreeIter));
  assert(tree);
  (void) start(trav, tree, 0);
  return trav;
}

RBTreeIter* rbtree_iter_new_from_last(const RBTree *tree)
{
  RBTreeIter *trav = malloc(sizeof (RBTreeIter));
  assert(tree);
  (void) start(trav, tree, 1);
  return trav;
}

void rbtree_iter_reset_from_first(RBTreeIter *trav)
{
  assert(trav);
  (void) start(trav, trav->tree, 0);
}

void rbtree_iter_reset_from_last(RBTreeIter *trav)
{
  assert(trav);
  (void) start(trav, trav->tree, 1);
}

void *rbtree_iter_next(RBTreeIter *trav)
{
  assert(trav);
  return move(trav, 1, false); /* Toward larger items */
}

void *rbtree_iter_prev(RBTreeIter *trav)
{
  assert(trav);
  return move(trav, 0, false); /* Toward smaller items */
}

void *rbtree_iter_data (RBTreeIter *trav)
{
  RBTreeNode *node;
  assert(trav);
  node = trav->it;
  return (node != NULL ? node->key : NULL);
}

void rbtree_iter_delete (RBTreeIter *trav)
{
  free(trav);
}
