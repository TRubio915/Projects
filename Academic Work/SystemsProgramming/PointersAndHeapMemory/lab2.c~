#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct freq
{
  char* value;
  int result;
}Histogram;

char** readScores(int* c);
void displayScores(char** a, int* c);
int calcHistogram(char** in, Histogram** out, int c);
void displayHistogram(Histogram* out, int c);
void sortHistogram(Histogram* arr, int* c);
void swap(Histogram* a, Histogram* b);


int main()
{
  int count = 0; // Counter variable for each digit read in

  char** inArr = readScores(&count);
  displayScores(inArr, &count);

  Histogram* outArr;

  count = calcHistogram(inArr,&outArr,count);
  displayHistogram(outArr,count); 
  sortHistogram(outArr,&count);
  displayHistogram(outArr,count);
  for(int i = 0; i < 100; i++)
    {
      free(*(inArr + i));
    }
  free(inArr);
  free(outArr);
}

//Reads scores from input .txt file
char** readScores(int* c)
 {
   char** a = (char**)malloc(100*sizeof(char*)); //Create array on heap for read-in
   FILE* filePoint = fopen("test2.txt", "r"); //Open .txt file
   if(filePoint == NULL)
     {
       printf("Error: File did not open properly!\n");
       exit(0);
     }
   char* word = (char*)malloc(10*sizeof(char)); //hold each value read in
   
   //Begin reading file until EOF
   while(scanf("%s",word) != EOF)
     {
       *(a + *c) = word;
       (*c)++;
       word = (char*)malloc(10*sizeof(char)); //hold each value read in
     }
   fclose(filePoint);
   free(word);
   return a;
 }

//Displays Scores read in from .txt file
void displayScores(char** a, int* c)
{
  printf("\nScores Read:\n");
  for(int i = 0; i < *c; i++)
    {
      printf("Word %d: %s\n", i, a[i]);
    }
  printf("\n");
}

//Calculates histogram and places into array of structs
int calcHistogram(char** in, Histogram** o, int c)
{
  Histogram* out = (Histogram*)malloc(c*sizeof(Histogram)); //creating array on heap
  int j = 0; //index for while loop
  int outSize = 0; //Tracks size of Output array
  Histogram temp; //Temp variable to hold each value in array
  //begin calculating array
  for(int i = 0; i < c; i++)
    {
      j=0; //Reset index
      temp.value = *(in + i); //set temp value as next in input array
      
      //loop while index is less than size, and temp value not found
      while(j<outSize && strcmp(temp.value, out[j].value) != 0)
	{
	  j++;
	}
      //Array end has been reached, add value and set result to 1
      if(j == outSize)
	{
	  temp.result = 1;
	  *(out + j) = temp;
	  outSize++;
	}
      //Value found in array, increment result value
      else
	{
	  out[j].result++;
	}
    }
  *o = out;
  return outSize;
}

//Prints calculated histogram array
void displayHistogram(Histogram* out, int c)
{
  printf("\nCalculated Histogram:\n");
  for(int i = 0; i < c; i++)
    {
      printf("Word - %s: freq %d\n", out[i].value, out[i].result);
    }
  printf("\n");
}

//Swap function used for histogramSort
void swap(Histogram* a, Histogram* b)
{
  Histogram temp = *a;
  *a = *b;
  *b = temp;
}

//Sorts histogram using Selection Sort algorithm modified
//to work with our struct
void sortHistogram(Histogram* arr, int* c)
{
  int minIndex;
  int i;
  int j;

  for(i = 0; i < *c-1; i++)
    {
      minIndex = i;
      for(j = i + 1; j < *c; j++)
	{
	  if(arr[j].result < arr[minIndex].result)
	    {
	      minIndex = j;
	    }
	}
      swap(&arr[minIndex], &arr[i]);
    }
}




      
