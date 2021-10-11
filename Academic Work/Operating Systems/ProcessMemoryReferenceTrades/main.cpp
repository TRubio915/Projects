/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: trubio915
 *
 * Created on October 8, 2020, 12:06 PM
 */

#include <cstdlib>

#include "Trace.h"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) 
{
    Trace *trace = new Trace(argv[1]);
    trace->RunTrace();
    return 0;
}

