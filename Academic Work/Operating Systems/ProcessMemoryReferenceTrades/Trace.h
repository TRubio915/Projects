/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Trace.h
 * Author: trubio915
 *
 * Created on October 8, 2020, 12:07 PM
 */

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include <string.h>

using namespace std;

#ifndef TRACE_H
#define TRACE_H

class Trace {
public:
    Trace(string fileName); //Constructor
    Trace(const Trace& other) = delete; //Copy Constructor
    Trace(Trace&& other) = delete; //Move Constructor
    Trace &operator = (const Trace &) = delete; //Copy Assignment
    Trace &operator = (Trace&&) = delete; //Move Assignment
    virtual ~Trace(); //Destructor
    void RunTrace(); //Runs trace of file
private:
    void ParseHex(string code, istringstream &cmnd); //Decides how to handle each bit of hex read
    ifstream in; //Input stream
    vector<uint8_t> memRef; //Vector to store memory references
};

#endif /* TRACE_H */

