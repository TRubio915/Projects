#ifndef DUCK_H
#define DUCK_H
#include <iostream>
#include <string>
#include "Animal.h"

class Duck : public Animal
{
 private:
 public:
  Duck();
  Duck(std::string color);
  Duck(const Duck& d);
  ~Duck();
};
#endif
