#ifndef BITMAPALLOCATOR_H
#define BITMAPALLOCATOR_H

#include <vector>
#include <iomanip>

class BitMapAllocator {
public:
    BitMapAllocator(int numPgFrames);
    BitMapAllocator(const BitMapAllocator& orig) = delete; //Copy
    BitMapAllocator(BitMapAllocator&& other) = delete; //Move
    BitMapAllocator &operator = (const BitMapAllocator &) = delete; //Copy Assignment
    BitMapAllocator &operator = (BitMapAllocator&&) = delete; //Move Assignment
    virtual ~BitMapAllocator(); //Destructor
    bool GetFrames(uint32_t count, std::vector<uint32_t> &pageFrames); //Gets  available frames for program
    bool FreeFrames(uint32_t count, std::vector<uint32_t> &pageFrames); //Frees specified number of frames
    uint32_t GetFreeCount() const; //Gets number of available frames
    std::string GetBitMapString()const; //Gets string of bitmap
    void setBit(int frameNum, int value); //Sets specific bit to 0 or 1
    int getBit(int frameNum); //Gets value of specific bit
    
private:
    std::vector<uint8_t> memory; //Vector to simulate memory
    const int PAGESIZE = 1024; //Const to hold size of each page frame
};

#endif /* BITMAPALLOCATOR_H */

