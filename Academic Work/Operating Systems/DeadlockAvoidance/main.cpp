#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <string.h>
#include <fstream>

using namespace std;

typedef struct ProcessStruct
{
    string name;
    int has; //Instances currently assigned to process
    int max; //Max # of instances that process needs
    int complete; //Determines if process has successfully completed
}Process;

bool isSafe(int available, int size, Process *pArr);

int main(int argc, char** argv) 
{
    const int SIZE = 100; //Size of process array
    ifstream in;
    in.open(argv[1]); //Open file
    string line; //holds entire line
    string s; //used for parsing lines
    stringstream stream; //String Stream
    bool firstLine = true;
    int available; //Saves number of available instances
    Process pArr[SIZE]; //Array of process structs
    int index = 0; //index into Process Array
    int x; //Holds integer values read from file
    bool safe; //Holds result of isSafe function
    
    //Check if file opened properly before running
    if(in.is_open())
    {
        //While not at EOF
        while(getline(in, line))
        {
            stringstream ss;
            
            //If we are on the first line of the project, save # of available
            //instances of the single resource
            if(firstLine)
            {
                ss << dec << line;
                ss >> available;
                firstLine = false;
            }
            else
            {
                //Parse name and assign
                ss << line;
                ss >> pArr[index].name;
                //Parse program's allocated instances of our resource
                ss << dec << line;
                ss >> x;
                pArr[index].has = x;
                //Parse max # of instances that process needs
                ss << dec << line;
                ss >> x;
                pArr[index].max = x;
                //Mark as incomplete
                pArr[index].complete = false;
                //Increment index of our process array, this variable will later
                //serve as number of occupied spaces in array
                index++;
            }
        }
    }
    else //If file opened incorrectly
    {
        cout << "File Failed to Open!" << endl;
    }
    in.close();
    
    //Determine if safe or unsafe
    cout << argv[1] << ": ";
    safe = isSafe(available, index, pArr);
    if(safe)
    {
        cout << " (Safe)" << endl;
    }
    else
    {
        cout << " (Unsafe)" << endl;
    }
    
    return 0;
}

bool isSafe(int available, int size, Process *pArr)
{
    int diff; //holds difference between has and max for each process
    bool success = true; //Used to determine if safe or unsafe, controls loop
    bool found; //Determines if a runnable process was found.
    int done; //Counts number of finished processes
    while(success)
    {
        found = false; //reset found
        for(int i = 0; i < size; i++)
        {
            //If process is complete, skip to next
            if(pArr[i].complete)
            {
                continue;
            }
            //Calculate needed resources
            diff = (pArr[i].max - pArr[i].has);
            if(diff <= available)
            {
                available -= diff; //Subtract resources
                pArr[i].complete = true; //Mark this process as complete
                found = true; //Found a process that can run
                done++; //Increment count of finished processes
                available += pArr[i].max; //Return all resources to available
                cout << pArr[i].name << " "; //Print completed process name
                break; //Break for loop so it restarts at beginning
            }
        }
        if(done == size)
        {
            //safe
            return true;
        }
        if(!found)
        {
            //unsafe
            success = false; //This is just to ensure infinite loop is impossible.
        }
    }
    //Return unsafe
    return false;
}

