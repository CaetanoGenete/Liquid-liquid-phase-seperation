#ifndef LLPS_ALIGNED_ALLOCATOR_HPP_INCLUDED
#define LLPS_ALIGNED_ALLOCATOR_HPP_INCLUDED

#ifdef LLPS_USE_MKL
#include <memory>

#include "fftw/fftw3.h"

namespace llps {

    template<class Type>
    struct fftw_allocator : public std::allocator<Type>
    {
    public:
        Type* allocate(size_t n)
        {
            return static_cast<Type*>(fftw_malloc(sizeof(Type) * n));
        }
        void deallocate(Type* type, size_t n)
        {
            fftw_free(type);
        }
    };

}

#endif // LLPS_USE_MKL

#endif // !LLPS_ALIGNED_ALLOCATOR_HPP_INCLUDED