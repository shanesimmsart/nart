#include "../../include/nart/core/memory.h"

MemoryArena::MemoryArena(size_t blockSize)  // ~2KB by default
{
    // Aligning to 16 bytes for cache-friendliness
    currentSize = (blockSize + 15) & (~15);
    currentBlockAddress = (uint8_t*)::operator new(currentSize);
}

void* MemoryArena::Allocate(size_t allocSize) {
    allocSize = (allocSize + 15) & (~15);

    if (allocSize > (currentSize - currentOffset)) {
        // Put currentBlock in used blocks
        usedBlocks.push_back(
            std::pair<uint8_t*, size_t>(currentBlockAddress, currentSize));
        currentOffset = 0;
        currentBlockAddress = nullptr;

        // Search for an available block of appropriate size
        for (auto block : availableBlocks) {
            if (allocSize < block.second) {
                currentBlockAddress = block.first;
                currentSize = block.second;
                availableBlocks.erase(std::remove(availableBlocks.begin(),
                                                  availableBlocks.end(), block),
                                      availableBlocks.end());
                break;
            }
        }

        // If no block found, create a new one
        if (!currentBlockAddress) {
            // Defaulting to 2KB
            currentSize = (allocSize + 2047) & (~2047);
            currentBlockAddress = (uint8_t*)::operator new(currentSize);
        }
    }

    uint8_t* allocAddress = currentBlockAddress + currentOffset;
    currentOffset += allocSize;
    return allocAddress;
}

void MemoryArena::Refresh() {
    currentOffset = 0;

    for (auto block : usedBlocks) {
        availableBlocks.push_back(block);
    }

    usedBlocks.clear();
}

MemoryArena::~MemoryArena() {
    ::operator delete(currentBlockAddress);

    for (auto block : usedBlocks) {
        ::operator delete(block.first);
    }

    for (auto block : availableBlocks) {
        ::operator delete(block.first);
    }

    usedBlocks.clear();
    availableBlocks.clear();
}
