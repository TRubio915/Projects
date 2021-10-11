#ifndef ANIMAL_H
#define ANIMAL_H
#include <iostream>

class Animal
{
 private:
  std::string color;
  std::string sound;
 public:
  Animal();
  Animal(std::string color, std::string sound);
  Animal(const Animal& a);
  virtual ~Animal();
  virtual std::string getColor();
  virtual void makeSound();
};
#endif
