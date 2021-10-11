/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Trace.cpp
 * Author: trubio915
 * 
 * Created on October 8, 2020, 12:07 PM
 */

#include "Trace.h"

Trace::Trace(string fileName) 
{
    //Open input file
    in.open(fileName);
}

Trace::~Trace() 
{
    //Close Input file
    in.close();
}

/**
 * RunTrace - Trace text file
 * 
 * @throws      n/a
 * @return      void
 */
void Trace::RunTrace()
{
    string s; //used for parsing lines
    string s2; //Used to hold parsed strings
    int count = 1; //Counter for lines read
    bool isHex = true; //Flags line as hex or regular text to print
    //Check if file opened properly before running
    if(in.is_open())
    {
        //While not at EOF
        while(getline(in, s))
        {
            cout << count << ":";
            istringstream line(s); //stream to hold each line read
            
            //While not at end of line
            while(line >> s2)
            {
                if(s2 == "*")
                {
                    isHex = false;
                }
                if(isHex)
                {
                   ParseHex(s2, line);
                   break;
                }
                else
                {
                    cout << s2 << " ";
                }
            }
            cout << endl;
            isHex = true;
            count++;
        }
    }
    else //If file opened incorrectly
    {
        cout << "File Failed to Open!" << endl;
    }
}

/**
 * ParseHex - parse line identified as hexidecimal values
 * 
 * @param code      code identifying action to take
 * @param cmnd      the rest of line, containing commands to be used
 * @throws          n/a
 * @return          void
 */
void Trace::ParseHex(string code, istringstream &cmnd)
{
    string s; //Holds commands from istream
    stringstream ss; //For conversion
    int x; //Holds integer conversion
    cout << code << " "; 
    
    if(code == "F00" || code == "00f00")
    {
        //Convert to integers and assign to appropriate variables
        cmnd >> s;
        cout << s << " ";
        ss << hex << s;
        ss >> x;
        //Add 0's up to specified count
        for(int i = memRef.size(); i < x; i++)
        {
            memRef.push_back(0);
        }
    }
    
    else if(code == "CB1" || code == "00cb1")
    {
        //Convert to integers and assign to appropriate variables
        int addr;
        uint8_t value;
        for(int i = 0; i < 2; i++)
        {
            cmnd >> s;
            cout << s << " ";
            ss << hex << s;
            ss >> x;
            if(i == 0)
            {
                addr = x;
            }
            else
            {
                value = x;
            }
        }
        //Check if values from addr to end of vector all match value
        //Print error if not
        for(int i = addr; i < memRef.size(); i++)
        {
            if(memRef[i] != value)
            {
                cout << "compare error at address " << addr << ", expected "
                    << value << ", actual is " << memRef[i] << endl;
            }
        }
    }
    
    else if(code == "cba" || code == "CBA" || code == "00000cba")
    {
        //Convert to integers and assign to appropriate variables
        int count;
        int addr;
        uint8_t value;
        for(int i = 0; i < 3; i++)
        {
            cmnd >> s;
            cout << s << " ";
            ss << hex << s;
            ss >> x;
            if(i == 0)
            {
                count = x;
            }
            else if(i == 1)
            {
                addr = x;
            }
            else
            {
                value = x;
            }
        }
        //Starting from addr, check if count # of entries all match value
        //If any do not, print error message
        for(int i = addr; i < count; i++)
        {
            if(memRef[i] != value)
            {
                cout << "compare error at address " << addr << ", expected "
                    << value << ", actual is " << memRef[i] << endl;
            }
        }
    }
    
    else if(code == "301")
    {
        //Convert to integers and assign to appropriate variables
        int addr;
        uint8_t value;
        cmnd >> s;
        cout << s << " ";
        ss << hex << s;
        ss >> x;
        addr = x;
        auto it = memRef.begin() + addr; //Iterator for adding variables
        //While not end of line, add each value provided to vector starting at 
        while(cmnd >> s)
        {
            cout << s << " ";
            ss << hex << s;
            ss >> x;
            if((it == memRef.end()) || (addr > memRef.size()))
            {
                memRef.push_back(value);
            }
            else
            {
                memRef.insert(it,value);
                it++;
            }
        }
    }
    
    else if(code == "30A")
    {
        //Convert to integers and assign to appropriate variables
        int count;
        int addr;
        uint8_t value;
        for(int i = 0; i < 3; i++)
        {
            cmnd >> s;
            cout << s << " ";
            ss << hex << s;
            ss >> x;
            if(i == 0)
            {
                count = x;
            }
            else if(i == 1)
            {
                addr = x;
            }
            else
            {
                value = x;
            }
        }
        for(int i = addr; i < count; i++)
        {
            memRef[i] = value;
        }
    }
    
    else if(code == "31D" || code == "031d")
    {
        //Convert to integers and assign to appropriate variables
        int count;
        int dest;
        int src;
        int j = 0; //Variable to iterate through dest addresses
        for(int i = 0; i < 3; i++)
        {
            cmnd >> s;
            cout << s << " ";
            ss << hex << s;
            ss >> x;
            if(i == 0)
            {
                count = x;
            }
            else if(i == 1)
            {
                dest = x;
            }
            else
            {
                src = x;
            }
        }
        for(int i = src; i < count; i++)
        {
            memRef[dest+j] = memRef[i];
            j++;
        }
    }
    
    else if(code == "4F0" || code == "4f0")
    {
        //Convert to integers and assign to appropriate variables
        int addr;
        int count;
        for(int i = 0; i < 2; i++)
        {
            cmnd >> s;
            cout << s << " ";
            ss << hex << s;
            ss >> x;
            if(i == 0)
            {
                addr = x;
            }
            else
            {
                count = x;
            }
        }
        cout << endl << "Address: " << addr << ": ";
        for(int i = addr; i < memRef.size(); i++)
        {
            cout << memRef[i] << " ";
        }
    }
}

