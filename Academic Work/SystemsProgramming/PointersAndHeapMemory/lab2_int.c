#include <stdio.h>
#include <stdlib.h>

typedef struct freq
{
  int value;
  int result;
}Histogram;

int* readScores(int* c);
void displayScores(int* a, int* c);
int calcHistogram(int* in, Histogram** out, int c);
void displayHistogram(Histogram* out, int c);
void sortHistogram(Histogram* arr, int* c);
void swap(Histogram* a, Histogram* b);


int main()
{
  int count = 0; // Counter variable for each digit read in

  int* inArr = readScores(&count);
  displayScores(inArr, &count);

  Histogram* outArr;

  count = calcHistogram(inArr,&outArr,count);
  displayHistogram(outArr,count); 
  sortHistogram(outArr,&count);
  displayHistogram(outArr,count);
  free(inArr);
  free(outArr);
}

//Reads scores from input .txt file
int* readScores(int* c)
 {
   int* a = (int*)malloc(100*sizeof(int)); //Create array on heap for read-in vals
   FILE* filePoint = fopen("test.txt", "r"); //Open .txt file
   if(filePoint == NULL)
     {
       printf("Error: File did not open properly!\n");
       exit(0);
     }
   int num = 0; //hold each value read in
   
   //Begin reading file until EOF
   while(scanf("%d",&num) != EOF)
     {
       a[*c] = num;
       (*c)++;
     }
   fclose(filePoint);
   return a;
 }

//Displays Scores read in from .txt file
void displayScores(int* a, int* c)
{
  printf("\nScores Read:\n");
  for(int i = 0; i < *c; i++)
    {
      printf("Score %d: %d\n", i, a[i]);
    }
  printf("\n");
}

//Calculates histogram and places into array of structs
int calcHistogram(int* in, Histogram** o, int c)
{
  Histogram* out = (Histogram*)malloc(c*sizeof(Histogram)); //creating array on heap
  int j = 0; //index for while loop
  int outSize = 0; //Tracks size of Output array
  Histogram temp; //Temp variable to hold each value in array
  //begin calculating array
  for(int i = 0; i < c; i++)
    {
      j=0; //Reset index
      temp.value = in[i]; //set temp value as next in input array
      
      //loop while index is less than size, and temp value not found
      while(j<outSize && temp.value != out[j].value)
	{
	  j++;
	}
      //Array end has been reached, add value and set result to 1
      if(j == outSize)
	{
	  temp.result = 1;
	  out[j] = temp;
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
      printf("Value %d: freq %d\n", out[i].value, out[i].result);
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




      
