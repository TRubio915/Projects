#include <iostream>
#include "Dog.h"

//Constructors
Dog::Dog() : Animal()
{
  name = "";
}

Dog::Dog(std::string name, std::string color) : Animal(color, "Bark!")
{
  this->name = name;
}

//Cpy Const
Dog::Dog(const Dog& d) : Animal(d)
{
  this->name = d.name;
}

//Destructor
Dog::~Dog()
{
  //Empty
}

std::string Dog::getName()
{
  return name;
}
