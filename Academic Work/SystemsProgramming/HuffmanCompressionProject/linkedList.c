#include <stdio.h>
#include <stdlib.h>
#include "linkedList.h"

/*int main()
{
  LinkedList* l = llCreate();
  huffTree node;
  /*node.character = 'a';
  node.weight = 1;
  llAdd(&l,node);
  node.character = 'c';
  node.weight = 3;
  llAdd(&l,node);
  node.character = 'd';
  node.weight = 4;
  llAdd(&l,node);
  llDisplay(l);
  node.character = 'b';
  node.weight = 2;
  listAddInOrder(&l,&node);
  llDisplay(l);
  removeFirst(&l);
  llDisplay(l);
  removeFirst(&l);
  llDisplay(l);
  llFree(l);
}
*/


LinkedList* llCreate()
{
  return NULL;
}

int llIsEmpty(LinkedList* ll)
{
  return(ll == NULL);
}

void llDisplay(LinkedList* ll)
{
  LinkedList* p = ll;

  printf("[");

  while(p != NULL)
    {
      printf("%c, ", p->value->character);
      p = p -> next;
    }

  printf("]\n");
}


void llAdd(LinkedList** ll, huffTree newValue)
{
  /* //Create a new node
  LinkedList* newNode = (LinkedList*)malloc(1 * sizeof(LinkedList));
  newNode -> value = newValue;
  newNode -> next = NULL;

  //Find the end of LinkedList
  LinkedList* p = *ll;

  if(p == NULL)
    {
      *ll = newNode;
    }
  else
    {
      while(p -> next != NULL)
	{
	  p = p -> next;
	}
      //Attach to end
      p -> next = newNode;
    } 
}

void llFree(LinkedList* ll)
{
  LinkedList* p = ll;

  while(p != NULL)
    {
      LinkedList* oldP = p;
      p = p -> next;
      free(oldP);
      }*/
}

void listAddInOrder(LinkedList** ll, huffTree* newValue)
{
  //Create a new node
  LinkedList* newNode = (LinkedList*)malloc(1 * sizeof(LinkedList));
  newNode -> value = newValue;
  newNode -> next = NULL;
  int found = 0;

  //Find the end of LinkedList
  LinkedList* p = *ll;

  if(p == NULL)
    {
      *ll = newNode;
    }
  else
    {
      while(p -> next != NULL)
	{
	  //if new value's weight is greater than p's
	  if(newNode->value->weight > p->value->weight)
	    {
	      newNode -> next = p -> next;
	      p -> next = newNode;
	      found = 1;
	    }
	  p = p -> next;
	}
      //Attach to end if not found
      if(found == 0)
	{
	  p -> next = newNode;
	}
    }
}

huffTree* removeFirst(LinkedList** ll)
{
  LinkedList* p = *ll;
  LinkedList* n = *ll;
  n = n -> next;
  *ll = n;
  return p->value;
}



