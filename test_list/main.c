//doubly elemented list example adapted from
//  https://www.tutorialspoint.com/data_structures_algorithms/doubly_elemented_list_program_in_c.htm
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "list.h"

//insert element at the first location
void insertFirst(list_t **head, list_t **tail, int key, int data)
{
  list_t *element = (list_t*) malloc(sizeof(list_t)); //create a element
  element->key = key;
  element->data = data;
  if(*head==NULL) {
    element->prev = NULL;
    element->next = NULL;
    *head = element;
    *tail = element;//make it the head element
  }else{
    element->prev = NULL;
    element->next = *head;  //point it to old first element
    (*head)->prev = element;//update first prev element 
    *head = element; //point first to new first element
  }
}

//insert element at the tail location
void insertTail(list_t **head, list_t **tail, int key, int data)
{
  list_t *element = (list_t*) malloc(sizeof(list_t));//create a element
  element->key = key;
  element->data = data;	
  if(*head==NULL){
    element->prev = NULL;
    element->next = NULL;
    *head = element;
    *tail = element;//make it the *tail element
  }else {
    element->prev = *tail;//mark old *tail listnode as prev of new element
    element->next = NULL;
    (*tail)->next = element;//make element a new *tail element    
    *tail = element;//point *tail to new *tail listnode
  }
}

//delete first item
void deleteFirst(list_t **head, list_t **tail)
{
  list_t *element = *head;//save reference to first element
  if(*head==NULL) return;//empty list
  
  if((*head)->next == NULL) {//if only one element, list becomes empty after deletion
    *head = NULL;
    *tail = NULL;
  }else{//more than one element in the list
    (*head)->next->prev = NULL;
    *head = (*head)->next;
  }
  
  if(element != NULL) free(element);
}

//delete element at the tail location
void deleteTail(list_t **head, list_t **tail)
{
  list_t *element = *tail;//save reference to *tail element

  if((*head)->next == NULL) { //if only one element, list becomes empty after deletion
    *head = NULL;
    *tail = NULL;
  }else {//more than one element in the list
    (*tail)->prev->next = NULL;
    *tail = (*tail)->prev;
  }

  if(element != NULL) free(element);
}

//delete a element with given key
void delete(list_t **head, list_t **tail, int key)
{
  list_t* current = *head;
  while(current!=NULL && current->key != key){
    if(current->next == NULL) current = NULL;//if it is tail node
    else current = current->next;//move to next element
  }  

  if(current != NULL){//found a node matching
    if(current == *head) *head = (*head)->next;//change first to point to next element
    else current->prev->next = current->next;//bypass the current element

    if(current == *tail) *tail = current->prev;//change *tail to point to prev element
    else current->next->prev = current->prev;

    free(current);//free the listnode matched
  }
}

void insertAfter(list_t **head, list_t **tail, int key, int newKey, int data)
{
  list_t *current = *head;	

  while(current!=NULL && current->key != key) {//navigate through list
    if(current->next == NULL) current = false;//if it is tail listnode
    else current = current->next;//move to next element
  }
	
  list_t *newElement = (list_t*)malloc(sizeof(list_t));//create a element
  
  newElement->key = newKey;
  newElement->data = data;
  if(current == NULL){
    newElement->prev = *tail;
    newElement->next = NULL; 
    *tail = newElement; 
  }else{
    newElement->next = current->next;
    newElement->prev = current; 
    current->next->prev = newElement;
    current->next = newElement; 
  }	

}

//delete a element with given key, Don't forget to free it outside this routine
void modify_element(list_t **head, list_t **tail, int key, int din)
{
  //start from the first link
  list_t *current, *cursor;

  if(*head == NULL || din==0) current = NULL;

  current = *head;
  while(current->key != key) {		
    if(current->next == NULL) current = NULL;//if it is *tail node
    else current = current->next;//move to next link
  }

  if(current!=NULL){//found a match: current->key==key
    //now remove the current node from the list
    if(current == *head) *head = (*head)->next;//change first to point to next link
    else current->prev->next = current->next;//bypass the current link

    if(current == *tail) *tail = current->prev;//change last to point to prev link
    else current->next->prev = current->prev;

    //--------------------------------------
    current->data += din;//modify current node
    if(din>0){//search backward to the head
      cursor = current->prev;
      while(cursor!=NULL && cursor->data<current->data) cursor = cursor->prev;

      if(cursor==NULL) {//reach the head
	current->prev = NULL;
	current->next = *head;
	(*head)->prev = current;
	*head = current;//make it the head element
      }else{//cursor->data>=current->data, insert now
	current->prev = cursor;
	current->next = cursor->next;  //point it to old first element
	cursor->next->prev = current;
	cursor->next = current;//update first prev element
      }

    }else if(din<0){//search forward to the tail
      cursor = current->next;
      while(cursor!=NULL && cursor->data>current->data) cursor = cursor->next;

      if(cursor==NULL){
	current->prev = *tail;
	current->next = NULL;
	(*tail)->next = current;
	*tail = current;//make it the tail element
      }else{//cursor->data<=current->data, insert now
	current->prev = cursor->prev;
	current->next = cursor;  //point it to old first element
	cursor->prev->next = current;
	cursor->prev = current;//update first prev element
      }
    }
  }  
}


// function to insert a new listnode in sorted way in a sorted doubly elemented list
void sorted_insert(list_t** head, list_t** tail, int key, int data)
{
  list_t **current = head;
  while(*current!=NULL && (*current)->data>data) current = &(*current)->next;

  list_t *elem = (list_t*) malloc(sizeof(list_t)); //create an element
  elem->key = key;
  elem->data = data;
  
  if(*current==NULL){//reach the end of the list
    elem->prev = *tail;
    elem->next = NULL;
    (*tail)->next = elem;
    *tail = elem;//make it the tail element
  }else{//(*current)->data<=data, place elem in front of current node
    elem->prev = (*current)->prev;
    elem->next = *current;//point it to old first element
    (*current)->prev = elem;//update first prev element    
    *current = elem; //point first to new first element
  }

}
 

//display the list in from first to tail
void displayForward(list_t *head)
{
  list_t *ptr = head;
  printf("\n forward: [ ");
  while(ptr != NULL) {        
    printf("(%d,%d) ",ptr->key,ptr->data);
    ptr = ptr->next;
  }
  printf(" ]\n");
}

//display the list from tail to first
void displayBackward(list_t *tail)
{
  list_t *ptr = tail;
	
  printf("\n backward: [ ");
  while(ptr != NULL) {
    printf("(%d,%d) ", ptr->key, ptr->data);
    ptr = ptr ->prev;
  }
  printf(" ]\n");
}

int main()
{
  list_t *head = NULL;//this element always point to first Element
  list_t *tail = NULL;//this element always point to tail Element 

  insertFirst(&head, &tail, 0, 10);
  sorted_insert(&head, &tail, 1, 10);
  sorted_insert(&head, &tail, 2, 72);
  insertFirst(&head, &tail, 7, 200);
  sorted_insert(&head, &tail, 3, 15);
  sorted_insert(&head, &tail, 4, 25);
  sorted_insert(&head, &tail, 5, 17);
  sorted_insert(&head, &tail, 6, 98);

  printf("\n initial list:\n");
  displayForward(head);
  displayBackward(tail);
  
  printf("\n modify key 5:\n");
  modify_element(&head, &tail, 5, +99);
  displayForward(head);
  displayBackward(tail);

  printf("\n modify key 2:\n");
  modify_element(&head, &tail, 2, -100);
  displayForward(head);
  displayBackward(tail);

  printf("\n after deleting first record: \n");
  deleteFirst(&head, &tail);        
  displayForward(head);
  displayBackward(tail);

  printf("\n after deleting tail record: \n");  
  deleteTail(&head, &tail);
  displayForward(head);
  displayBackward(tail);

  printf("\n insert after key(4): \n");  
  insertAfter(&head, &tail, 4, 7, 13);
  displayForward(head);
  displayBackward(tail);

  printf("\n after deleting key(7): \n");  
  delete(&head, &tail, 7);
  displayForward(head);
  displayBackward(tail);
  printf("\n");
  
  printf("\n after deleting key(33, nonexist): \n");  
  delete(&head, &tail, 33);
  displayForward(head);
  displayBackward(tail);

  return 0;
}
