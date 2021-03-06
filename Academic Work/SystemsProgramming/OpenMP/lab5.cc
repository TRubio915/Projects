#include <iostream>
#include <omp.h>
#include <mutex>

bool prime(int num);
int primeFind(int beg, int end);

void primeBlock(int beg, int end)
{
  int result = 0;
  double sTime;
  double eTime;
  int tNum;

  printf("BLOCKING\n");
  omp_set_num_threads(5);
#pragma omp parallel reduction(+:result) 
  {
    sTime = omp_get_wtime();
#pragma omp for nowait
    for(int i = beg; i<=end; i++)
      {
        if(prime(i))
         {
           result++;
         }
      }
    tNum = omp_get_thread_num();
    eTime = omp_get_wtime();
    eTime = eTime - sTime;
    printf("Time for %d: %f with %d found\n", tNum, eTime, result);
  }
  printf("From: %d To: %d There is: %d Prime Numbers!\n\n", beg, end, result);
}

void primeStripe(int beg, int end)
{
  int result = 0;
  double sTime;
  double eTime;
  int tNum;
  
  printf("\nSTRIPING\n");
  
  omp_set_num_threads(5);
  
#pragma omp parallel reduction(+:result)
  {
    sTime = omp_get_wtime();
#pragma omp for schedule(static, 1) nowait
    for(int i = beg; i<=end; i++)
      {
	if(prime(i))
	  {
	    result++;
	  }
      }
    tNum = omp_get_thread_num();
    eTime = omp_get_wtime();
    eTime = eTime - sTime;
    printf("Time for %d: %f with %d found\n", tNum, eTime, result);
  }
  printf("From: %d To: %d There is: %d Prime Numbers!\n\n", beg, end, result);
}

int main()
{
  int beg = 1000;
  int end = 1000000;
  primeStripe(beg, end);
  primeBlock(beg, end);
}

bool prime(int num)
{
  int temp = 2;
  bool isPrime = true;
  while(temp != num)
    {
      //is prime
      if(num % temp == 0)
	{
	  isPrime = false;
	  temp = num;
	}
      else
	{
	  temp++;
	}
    }
  return isPrime;
}

int primeFind(int beg, int end)
{
  int result = 0;
  for(int i = beg; i<=end; i++)
    {
      if(prime(i))
	{
	  result++;
	}
    }
  return result;
}
