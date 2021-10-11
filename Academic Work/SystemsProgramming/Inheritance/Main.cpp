#include <iostream>
#include "Dog.h"
#include "Cat.h"
#include "Duck.h"

int main()
{
  Dog* d1 = new Dog("Fido", "brown");
  std::cout << d1->getName( );
  Animal* a2 = new Cat("Whiskers", "tabby");
  a2->makeSound( );
  Duck* du1 = new Duck("green");
}
