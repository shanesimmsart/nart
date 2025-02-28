#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

class MemoryArena {
public:
    MemoryArena(size_t blockSize = 2048);  // ~2KB by default

    void* Allocate(size_t allocSize);

    void Refresh();

    ~MemoryArena();

private:
    uint8_t* currentBlockAddress = nullptr;
    size_t currentOffset = 0;
    size_t currentSize = 0;
    std::vector<std::pair<uint8_t*, size_t>> usedBlocks;
    std::vector<std::pair<uint8_t*, size_t>> availableBlocks;
};
