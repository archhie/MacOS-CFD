#pragma once

#include <cstdlib>
#include <cstring>
#include <new>
#include <utility>

// Simple aligned 2D field with ghost cells.
template <typename T> struct Field2D {
    int nx = 0, ny = 0, pitch = 0, ngx = 0, ngy = 0;
    T *data = nullptr;

    Field2D() = default;
    ~Field2D() { free(); }

    Field2D(const Field2D &) = delete;
    Field2D &operator=(const Field2D &) = delete;

    Field2D(Field2D &&other) noexcept { *this = std::move(other); }
    Field2D &operator=(Field2D &&other) noexcept {
        if (this != &other) {
            free();
            nx = other.nx;
            ny = other.ny;
            pitch = other.pitch;
            ngx = other.ngx;
            ngy = other.ngy;
            data = other.data;
            other.data = nullptr;
            other.nx = other.ny = other.pitch = other.ngx = other.ngy = 0;
        }
        return *this;
    }

    void allocate(int NX, int NY, int pitch_, int ngx_, int ngy_) {
        free();
        nx = NX;
        ny = NY;
        pitch = pitch_;
        ngx = ngx_;
        ngy = ngy_;
        size_t total = static_cast<size_t>(pitch) * (ny + 2 * ngy);
#if defined(_MSC_VER)
        data = static_cast<T *>(_aligned_malloc(total * sizeof(T), 64));
        if (!data)
            throw std::bad_alloc();
#else
        if (posix_memalign(reinterpret_cast<void **>(&data), 64,
                           total * sizeof(T)) != 0)
            data = nullptr;
        if (!data)
            throw std::bad_alloc();
#endif
        std::memset(data, 0, total * sizeof(T));
    }

    void free() {
#if defined(_MSC_VER)
        if (data)
            _aligned_free(data);
#else
        if (data)
            std::free(data);
#endif
        data = nullptr;
        nx = ny = pitch = ngx = ngy = 0;
    }

    inline T &at_raw(int ii, int jj) { return data[jj * pitch + ii]; }
    inline const T &at_raw(int ii, int jj) const {
        return data[jj * pitch + ii];
    }
};

