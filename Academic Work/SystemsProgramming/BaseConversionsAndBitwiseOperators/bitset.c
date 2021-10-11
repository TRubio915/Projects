#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef unsigned short bitSet;

bitSet makeBitSet();
void displayBitSet(bitSet bs);
void setBit(bitSet* bs, int index);
void clearBit(bitSet* bs, int index);
int bitValue(bitSet bs, int index);

int main()
{
  bitSet bs = makeBitSet();
  
  displayBitSet(bs);

  setBit(&bs,2);

  displayBitSet(bs);

  //clearBit(&bs,2);

  //displayBitSet(bs);
  
}

bitSet makeBitSet()
{
  return 0;
}

void displayBitSet(bitSet bs)
{
  for(int i = 15; i >= 0; i--)
    {
      printf("%d", bitValue(bs,i));
    }
  printf("\n");
}

void setBit(bitSet* bs, int index)
{
  bitSet bit = 1 << index;
  *bs = *bs | bit;
}

void clearBit(bitSet* bs, int index)
{
  bitSet bit = ~(1 << index);
  *bs = *bs & bit;
}

int bitValue(bitSet bs, int index)
{
  bs = bs >> index;
  bs = bs & 1;
  return bs;
}
