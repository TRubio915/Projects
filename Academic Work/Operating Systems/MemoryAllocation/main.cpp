#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include <string.h>

#include "BitMapAllocator.h"

using namespace std;


int main(int argc, char** argv) 
{
    vector<vector<uint32_t>> processes(4); //Vector to simulate pages in each process
    ifstream in;
    in.open(argv[0]);
    string s; //used for parsing lines
    string s2; //Used to hold parsed strings
    stringstream stream;
    int x; //Holds int values when reading from file
    
    //create BitMapAllocator object of size x
    BitMapAllocator *bma;
    
    //Check if file opened properly before running
    if(in.is_open())
    {
        int pNum; //Holds process number
        int count; //Holds count for functions
        //First line is number of page frames
        getline(in, s);
        //convert to decimal
        stream << dec << s;
        stream >> x;
        
        //Create new BMA object of size x
        bma = new BitMapAllocator(x);
        
        //While not at EOF
        while(getline(in, s))
        {
            istringstream line(s); //stream to hold each line read
            while(line >> s2)
            {
                if(s2 == "G")
                {
                    line >> s2;
                    stream << dec << s;
                    stream >> pNum;
                    line >> s2;
                    stream << dec << s;
                    stream >> count;
                    bma->GetFrames(count, processes[pNum]);
                }
                else if(s2 == "F")
                {
                    line >> s2;
                    stream << dec << s;
                    stream >> pNum;
                    line >> s2;
                    stream << dec << s;
                    stream >> count;
                    bma->FreeFrames(count, processes[pNum]);
                }
                else if(s2 == "B")
                {
                    string bmS = bma->GetBitMapString();
                    cout << s2 << endl;
                }
            }
        }
    }
    else //If file opened incorrectly
    {
        cout << "File Failed to Open!" << endl;
    }
    in.close();
    return 0;
}

