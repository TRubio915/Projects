#ifndef DOG_H
#define DOG_H
#include <iostream>
#include <string>
#include "Animal.h"

class Dog : public Animal
{
 private:
  std::string name;                                                                    public:
  Dog();
  Dog(std::string name, std::string color);
  Dog(const Dog& d);
  ~Dog();
  std::string getName();
};
#endif
