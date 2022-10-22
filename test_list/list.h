#ifndef __list_h__
#define __list_h__

typedef struct list_node{
  int data;
  int key;
  struct list_node *next;
  struct list_node *prev;
} list_t;

#endif
