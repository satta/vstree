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

#ifndef RBTREE_H
#define RBTREE_H

#include "types.h"

typedef enum
{
  RBTREE_PREORDER,
  RBTREE_POSTORDER,
  RBTREE_ENDORDER,
  RBTREE_LEAF
} RBTreeContext;

/* The <RBTree> class. Fast logarithmic data structure. This implementation
   does not allow storage of duplicates. */
typedef struct RBTree RBTree;
typedef struct RBTreeIter RBTreeIter;

typedef int   (*CompareWithData)(const void*, const void*, void *data);
typedef int   (*RBTreeAction)(void *key, RBTreeContext, Uint, void*);
typedef void  (*RBTreeFreeFunc)(void *p);

/* Returns a new <RBTree> object. <free> might be NULL and will be used
   to free key-object otherwise. <info> is the data for the <cmp>-function. */
RBTree*      rbtree_new(CompareWithData cmp, RBTreeFreeFunc free,
                             void *info);
void           rbtree_delete(RBTree *tree);
/* Deletes all tree elements */
void           rbtree_clear(RBTree *tree);

/* Returns <key> if element was found in <tree> and NULL if not */
void*          rbtree_find(const RBTree *tree, void *key);
void*          rbtree_find_with_cmp(const RBTree *tree, void *key,
                                       CompareWithData cmpfunc, void *info);
/* inserts <key> into <tree>. If <key> is already present in <tree>, it will not
   be changed. */
void           rbtree_insert(RBTree *tree, void *key);
void           rbtree_insert_with_cmp(RBTree *tree, void *key,
                                         CompareWithData cmpfunc,
                                         void *info);
/* Returns <key>, if <key> is not present in <tree>
   it will be inserted and <nodecreated> set accordingly */
void*          rbtree_search(RBTree *tree, void *key, bool *nodecreated);
void*          rbtree_search_with_cmp(RBTree *tree, void *key,
                                         CompareWithData cmpfunc,
                                         void *info, bool *nodecreated);
/* Remove <key> from <tree>, returns -1 if no such key exists and 0 on success
*/
int            rbtree_erase(RBTree *tree, void *key);
size_t         rbtree_size(const RBTree *tree);
int            rbtree_walk(const RBTree *tree, RBTreeAction action,
                              void *actinfo);
int            rbtree_walk_stop(const RBTree *tree, RBTreeAction action,
                                   void *actinfo);
int            rbtree_walk_reverse(const RBTree *tree,
                                      RBTreeAction action,
                                      void *actinfo);
void*          rbtree_minimum_key(const RBTree *tree);
void*          rbtree_maximum_key(const RBTree *tree);
void*          rbtree_root_key(const RBTree *tree);
void*          rbtree_next_key(const RBTree *tree, void *key,
                                  CompareWithData cmpfun,
                                  void *cmpinfo);
void*          rbtree_next_equal_key(const RBTree *tree, void *key,
                                        CompareWithData cmpfun,
                                        void *cmpinfo);
void*          rbtree_previous_key(const RBTree *tree, void *key,
                                      CompareWithData cmpfun,
                                      void *cmpinfo);
void*          rbtree_previous_equal_key(const RBTree *tree, void *key,
                                            CompareWithData cmpfun,
                                            void *cmpinfo);

RBTreeIter*  rbtree_iter_new_from_first(const RBTree *tree);
RBTreeIter*  rbtree_iter_new_from_last(const RBTree *tree);
/* Resets the iterator to the first (smallest) element */
void           rbtree_iter_reset_from_first(RBTreeIter *trav);
/* Resets the iterator to the last (largest) element */
void           rbtree_iter_reset_from_last(RBTreeIter *trav);
/* Return next (larger) key, NULL if iterator reached end */
void*          rbtree_iter_next(RBTreeIter *trav);
/* Return previous (smaller) key, NULL if iterator reached end */
void*          rbtree_iter_prev(RBTreeIter *trav);
/* Returns data of the current node the iterator <trav> is positioned on. */
void*          rbtree_iter_data (RBTreeIter *trav);
/* free all memory of <trav> */
void           rbtree_iter_delete(RBTreeIter *trav);

#endif
