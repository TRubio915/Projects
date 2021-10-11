#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "linkedList.h"

huffTree* createFreqTable(char* s)
{
  FILE* filePoint = fopen(s, "r"); //Open file
  if(filePoint == NULL)
    {
      printf("CFT Error: File did not open properly!\n");
      exit(0);
    }
  char c; //Holds each character read-in
  int max = 128; //Maximum based on # of ascii characters
  int size = 0; //Track size of table
  int found = 0; //Bool to determine if found
  huffTree* arr = (huffTree*)malloc(max*sizeof(huffTree));//Array of huffTree nodes

  //Beging reading file until EOF
  while(c = fgetc(filePoint) != EOF)
    {
      for(int i = 0; i <= size; i++)
	{
	  //Character exists in table
	  if(arr[i].character == c)
	    {
	      arr[i].weight++;
	      found = 1;
	      break;
	    }
	}
      //Character needs to be added
      if(found == 0)
	{
	  arr[size].character = c;
	  size++;
	}
      found = 0; //Reset found
    }
  for(int i = size-1; i<max; i++)
    {
      arr[i].weight = 0;
    }
  fclose(filePoint); //Close file
  //Return Table
  return arr;
}

huffTree* createHuffmanTree(huffTree* leafNodes)
{
  LinkedList* l = llCreate(); //Create linked list
  //Used for adding weights
  huffTree* n1;
  huffTree* n2;
  //Add leafNodes in-order to LL
  for(int i = 0; i < 128; i++)
    {
      listAddInOrder(&l,&(leafNodes[i]));
    }
  //Begin creating tree
  for(int i = 0; i<127; i++)
    {
      huffTree* new = (huffTree*)malloc(1*sizeof(huffTree)); //Creates new node
      //each pass
      n1 = removeFirst(&l);
      n2 = removeFirst(&l);
      //Add weights into new node
      new->weight = n1->weight + n2->weight;
      //Reset new node as parent of added nodes
      new->left = n1;
      new->right = n2;
      n1->parent = new;
      n2->parent = new;
      //Add new node in-order to LL
      listAddInOrder(&l,new);
    }
  //return root
  return removeFirst(&l);
}

void encodeFile(char* s, huffTree* leaves)
{
  FILE* filePoint = fopen(s, "r"); //Open file
  if(filePoint == NULL)
    {
      printf("EF Error: File did not open properly!\n");
      exit(0);
    }
  FILE* newPoint = fopen("decind.huf","w");//Create encoded file
  char c; //Holds each character read-in
  char byte[] = "00000000"; //Initial empty byte in string form
  int j = 0; //Index for byte char array
  
  //Beging reading file until EOF
  while(c = fgetc(filePoint) != EOF)
    {
      
      //If 8 bits have been processed, string write to file
      if(j > 7)
	{
	  fputs(byte, newPoint);
	  fputc(' ', newPoint);
	  for(int i = 0; i < 8; i++)
	    {
	      byte[i] = '0';
	    }
	}
      //Find leaf containing char and assign its address to new node pointer
      huffTree* node;
      for(int i = 0; i < 127; i++)
	{
	  //Leaf found
	  if(leaves[i].character == c)
	    {
	      node = &leaves[i];
	      break;
	    }
	}
      //Move to leafs parent
      node = node->parent;
      //If node's right child contains the char, bit is a '1'.
      //otherwise it remains '0'
      if(node->right->character == c)
	{
	  byte[j] = '1';
	}
      //Increment byte index
      j++;
    }
  fclose(filePoint);
}

int isLeaf(huffTree* node)
{
  if (node -> left == NULL && node -> right == NULL)
    {
      return 1;
    }
  else
    {
      return 0;
    }
}

void decodeFile(char* s, huffTree* node)
{
  FILE* filePoint = fopen(s, "r"); //Open file
  if(filePoint == NULL)
    {
      printf("DF Error: File did not open properly!\n");
      exit(0);
    }
  FILE* newPoint = fopen("decind.huf.dec","w");//Create decoded file
  char c; //Holds each character read-in
    
  //Beging reading file until EOF
  while(scanf("%c",&c) != EOF)
  {
    //break to beginning of loop if 'white space'
    if(c == ' ') 
      {
	break; 
      }
    //If '1', progress to rightChild
    else if(c == '1')
      { 
	node = node -> right;
	//If leaf, write character to file
	if(isLeaf(node)==1)
	  {
	    fputc(node->character, newPoint); 
	  } 
      }
    //if '0' progress to rightChild
    else if(c == '0') 
      { 
	node = node -> left;
	//If leaf, write character to file
	if(isLeaf(node)==1)
	  {
	    fputc(node->character, newPoint); 
	  } 
      } 
  }
  fclose(filePoint);
}

int main(int argc, char* argv[])
{
  //Check to ensure input parameters are correct
  if(argc != 3)
    {
      printf("Error: The correct format is \'hcompress -e filename\' or \'hcompress-d filename.huf\'\n"); fflush(stdout);

      exit(1);
    }

  //Create the frequency table by reading the generic file
  huffTree* leafNodes = createFreqTable("decind.txt");
  

  //Create the huffman tree from the freq table
  huffTree* treeRoot = createHuffmanTree(leafNodes);

  //encode
  if(strcmp(argv[1],"-e") == 0)
    {
      //Pass the leafNodes since it will process bottom-up
      encodeFile(argv[2],leafNodes);
    }
  //decode
  else
    {
      //Pass tree root since it will process top-down
      decodeFile(argv[2],treeRoot);
    }
  
  return 0;
}
