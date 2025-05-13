#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
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

template<typename T>
class ArenaAllocator {
public:
    using value_type = T;

    ArenaAllocator(MemoryArena* _memoryArena) : memoryArena(_memoryArena) {}

    template <typename U>
    ArenaAllocator(const ArenaAllocator<U>& alloc) : memoryArena(alloc.memoryArena) {}

    T* allocate(size_t n) { 
        return static_cast<T*>(memoryArena->Allocate(n * sizeof(T)));
    }

    // Do nothing
    // (Deallocation is done all at once in MemoryArena dtor)
    void deallocate(T* p, size_t n) {}

    // To allow templated copies access to arena ptr
    friend class ArenaAllocator;

private:
    MemoryArena* memoryArena;
};
