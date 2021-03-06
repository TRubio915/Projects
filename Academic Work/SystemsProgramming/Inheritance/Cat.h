#ifndef CAT_H
#define CAT_H
#include <iostream>
#include <string>
#include "Animal.h"

class Cat : public Animal
{
 private:
  std::string name;                                                                   
 public:
  Cat();
  Cat(std::string name, std::string color);
  Cat(const Cat& c);
  ~Cat();
  std::string getName();
};
#endif
