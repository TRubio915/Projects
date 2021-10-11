#include <iostream>
#include <thread>
#include <mutex>

std::mutex epiLock;

bool prime(int num);
int primeFind(int beg, int end);

int final;

void primeThread(int start, int stop)
{
  int result = 0;
  for(int i = start; i<=stop; i++)
    {
      if(prime(i))
	{
	  result++;
	}
    }
  epiLock.lock();
  final += result;
  epiLock.unlock();
}

int main()
{
  int beg = 10;
  int end = 100;
  printf("Without Threading: %d\n", primeFind(beg, end));
  int n = end-beg;
  int nth = 4;
  final = 0;

  std::thread* ths[nth];
  for(int i = 0; i<nth; i++)
    {
      int start = (n/nth)*i;
      int stop = (n/nth)*(i+1);
      start += beg;
      stop += beg;
      ths[i] = new std::thread(primeThread, start, stop);
    }
  for(int i=0; i<nth; i++)
    {
      ths[i]->join(); 
    }
  printf("With Threading: %d\n",final);
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
