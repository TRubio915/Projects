/* 
 * File:   BitMapAllocator.cpp
 * Author: trubio915
 * 
 * Created on October 28, 2020, 11:24 PM
 */

#include "BitMapAllocator.h"

BitMapAllocator::BitMapAllocator(int numPgFrames) 
{
    //Resize memory to number of PFs * size of a PF
    memory.resize(numPgFrames * PAGESIZE);
    
    //Build Bitmap
    uint32_t freeFrames = numPgFrames; //4 Bytes indicating the number of free frames available  
    
    //Initialize first 37 bytes of frame to 0
    for(int i = 0; i < 37; i++)
    {
        memory[i] = 0;
    }
    
    //Change numPgFrames number of bits to 1 to indicate they are available
    for(int i = 1; i <= numPgFrames; i++)
    {
        setBit(i, 1);
    }
}

BitMapAllocator::~BitMapAllocator() 
{}

bool BitMapAllocator::GetFrames(uint32_t count, std::vector<uint32_t> &pageFrames)
{
    //If # of available frames is less than count, return false
    if(count > GetFreeCount())
    {
        return false;
    }
    //while there are remaining page frames to allocate
    while(count > 0)
    {
        int i = 0; //Index variable
        //If bit indicates available page frame, push index onto back of vector,
        //Set this bit to unavailable.
        //Allocate PF and set all bytes to 0
        if(getBit(i) == 1)
        {
            pageFrames.push_back(i);
            setBit(i,0);
            for(int j = i; j < PAGESIZE; j++)
            {
                memory[j] = 0;
            }
            count--; //Decrement count until 0
        }
        i++;
    }
    
    return true;
}

bool BitMapAllocator::FreeFrames(uint32_t count, std::vector<uint32_t> &pageFrames)
{
    int frame; //temp variable to hold frame address when popped off back of vector
    
    //If count is less than number of pageFrames, pop off back of vector
    //and mark corresponding memory as available
    if(count <= pageFrames.size())
    {
        for(int i = count; i > 0; i--)
        {   
            frame = pageFrames.back();
            pageFrames.pop_back();
            setBit(frame, 1);
        }
        return true;
    }
    else
    {
        //Return false without releasing any page frames
        return false;
    }
}

uint32_t BitMapAllocator::GetFreeCount() const
{
    int available;
    for(int i = 0; i < 4; i++)
    {
        available += memory[i];
    }
    return available;
}

std::string BitMapAllocator::GetBitMapString()const
{
    std::string bMap; //Holds a string representation of the bitmap
    std::string byteStr; //Holds each string representation of bytes
    std::stringstream stream; //Used for hex conversion
    
    for(int i = 5; i < 37; i++)
    {
        stream << std::hex << memory[i];
        byteStr = stream.str();
        bMap += " " + byteStr;
    }
    return bMap;
}

//Sets specified bit in bitmap
void BitMapAllocator::setBit(int frameNum, int value)
{
    int index = frameNum / 8;
    int shift = frameNum % 8;
    memory.at(5 + index) = ((memory.at(5 + index) & ~(1 << shift)) | ((value & 1) << shift));
}
    
//Gets specified bit from bitmap
int BitMapAllocator::getBit(int frameNum)
{
    int index = frameNum / 8;
    int shift = frameNum % 8;
    
    return (memory.at(5 + index) >> shift) & 1;
}
