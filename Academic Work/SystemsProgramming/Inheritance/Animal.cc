#include <iostream>
#include "Animal.h"

//Constructors
Animal::Animal()
{
  color = "";
  sound = "";
}

Animal::Animal(std::string color, std::string sound)
{
  this->color = color;
  this->sound = sound;                                                                 }

//Copy constructor
Animal::Animal(const Animal& a)
{
  this->color = a.color;
  this->sound = a.sound;
}

//Destructor
Animal::~Animal()
{
  //Empty
}

std::string Animal::getColor()
{                                                                                        std::cout << "My Color is: " << color << "\n";
  return color;
}

void Animal::makeSound()
{
  std::cout << sound << "\n";
}
