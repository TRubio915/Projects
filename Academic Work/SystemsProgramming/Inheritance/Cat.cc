#include <iostream>
#include "Cat.h"

//Constructors
Cat::Cat() : Animal()
{
  name = "";
}

Cat::Cat(std::string name, std::string color) : Animal(color, "Bark!")
{
  this->name = name;
}

//Cpy Const
Cat::Cat(const Cat& c) : Animal(c)
{
  this->name = c.name;
}

//Destructor
Cat::~Cat()
{
  //Empty
}

std::string Cat::getName()
{
  return name;
}
