#include <iostream>
#include "Duck.h"

//Constructors
Duck::Duck() : Animal()
{
  //Empty
}

Duck::Duck(std::string color) : Animal(color, "Quack!")
{
  //Empty
}

//Cpy
Duck::Duck(const Duck& d) : Animal(d)
{
  //Empty
}

//Destructor
Duck::~Duck()
{
  //Empty
}
