#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int binToDec(char* bin);
char* decToBin(int dec);
int baseToDec(int base, char* value);
char* decToBase(int base, int dec);

int main()
{
  int num = binToDec("11001");
  printf("Bin to Dec: %d\n", num);

  char* bin = decToBin(num);
  printf("Dec to Bin: %s\n", bin);

  int dec = baseToDec(2, "11001");
  printf("Base 2 to Dec: %d\n", dec);

  dec = baseToDec(8, "157");
  printf("Base 8 to Dec: %d\n", dec);

  dec = baseToDec(16, "f8");
  printf("Base 16 to Dec: %d\n", dec);

  char* base = decToBase(8, 111);
  printf("Dec to Base 8: %s", base);
  for(int i = 0; i < 8; i++)
    {
      printf("%c", base[i]);
    }
  base = decToBase(16, 248);
  printf("\nDec to Base 16: %s", base);
  for(int i = 0; i < 8; i++)
    {
      printf("%c", base[i]);
    }
  printf("\n");
  
  free(bin);
  free(base);
}

int binToDec(char* bin)
{
  int i = (strlen(bin)-1); //Set i to length of binary string
  double exp = 0.0; //Exponent to be multiplied by 2
  int result = 0; //Resulting integer to be returned
  
  //Calculate integer from binary string
  while(i >= 0)
    {
      //If digit is a 1, add resulting power to result
      if(bin[i] == '1')
	{
	  result += (int)(pow(2, exp));
	}
      i--;
      exp++;
    }
  
  return result;
}

char* decToBin(int dec)
{
  int size = 8; //Size of our array to be determined by size of dec
  char* bin = (char*)malloc(size*sizeof(char*)); //Binary String to be returned
  double bit; //Tracks current bit's value
  int place = size - 1; //Tracks our current place in the byte

  //Compose binary string from given integer
  for(int i = 0; i<size; i++)
    {
      bit = pow(2,place);
      if(bit <= dec)
	{
	  dec -= bit;
	  bin[i] = '1';
	}
      else
	{
	  bin[i] = '0';
	}
      place--;
    }
  
  return bin;
}

int baseToDec(int base, char* value)
{
  int i = (strlen(value)-1); //Initialize index to the size of string - 1
  double exp = 0.0; //Exponent for conversion equation
  int result = 0; //Result to be returned
  int iVal; //Integer value of each char (bit)
  while(i >= 0)
    {
     iVal = value[i] - '0';
     //If iVal is greater than 16, our value is a-f
     if(iVal > 16)
       {
	 //Switch to determine iVal based on a-f in hex
	 switch(value[i])
	   {
	   case 'a':
	     iVal = 10;
	   case 'b':
	     iVal = 11;
	   case 'c':
	     iVal = 12;
	   case 'd':
	     iVal = 13;
	   case 'e':
	     iVal = 14;
	   case 'f':
	     iVal = 15;
	   }
       }
     //Calculation
     result += (iVal * (int)(pow(base, exp)));
     i--;
     exp++;
    }
  return result;
}

char* decToBase(int base, int dec)
{
  int size = 8; //Size of our array to be determined by size of dec
  char* bin = (char*)malloc(size*sizeof(char*)); //Binary String to be returned
  int bit; //Tracks current bit's value
  int place = size - 1; //Tracks our current place in the byte
  char cVal; //Holds char value to be place in array
  for(int i = 0; i<size; i++)
    {
      bit = pow(base,place);
      if(bit <= dec)
	{
	  if(base == 16 && (dec/bit) > 9)
	    {
	      switch (dec/bit)
		{
		case 10:
		  cVal = 'a';
		case 11:
		  cVal = 'b';
		case 12:
		  cVal = 'c';
		case 13:
		  cVal = 'd';
		case 14:
		  cVal = 'e';
		case 15:
		  cVal = 'f';
		}
	    }
	  else
	    {
	      cVal = ((dec / bit) + '0');
	    }
	  bin[i] = cVal;
	  dec = dec % bit;
	}
      place--;
    }
  
  return bin;
}
