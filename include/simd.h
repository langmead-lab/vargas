/**
 * Ravi Gaddipati
 * Jan 10, 2016
 * rgaddip1@jhu.edu
 *
 * @brief
 * SIMD wrapper class for SSE and AVX2.
 * @details
 * Note: element access uses reinterpret_cast<native_t*> by default, and is technically undefined behavior
 * for native_t != [unsigned] char. It seems for GCC the resulting
 * ASM is same for native_t = [unsigned] char using compliant code, but not for Intel compiler.
 * Compliant access option:
 *
 * @code{.cpp}
 * __RG_STRONG_INLINE__
 * native_t &operator[](const int i) {
 *     native_t buf;
 *     std::memcpy(&buf, reinterpret_cast<char>(&v) + (i*sizeof(native_t)), sizeof(native_t));
 *     return buf;
 * }
 * @endcode
 *
 * @file
 */

#ifndef VARGAS_SIMD_H
#define VARGAS_SIMD_H

#include "utils.h"
#include "doctest.h"

#include <type_traits>
#include <x86intrin.h>
#include <cstdint>
#include <memory>
#include <vector>
#include <stdlib.h>



#if !defined(VA_SIMD_USE_SSE) && !defined(VA_SIMD_USE_AVX2) && !defined(VA_SIMD_USE_AVX512)
#error("No SIMD instruction set defined.")
#endif

#ifdef VA_SIMD_USE_AVX512
#ifndef VA_SIMD_USE_AVX2
#define VA_SIMD_USE_AVX2
#endif
#endif

#ifdef VA_SIMD_USE_AVX2
#define VA_MAX_INT8 32
#define VA_MAX_INT16 16
#endif

#ifdef VA_SIMD_USE_SSE
#ifndef VA_MAX_INT8
#define VA_MAX_INT8 16
#endif
#ifndef VA_MAX_INT16
#define VA_MAX_INT16 8
#endif
#endif

namespace vargas {

  /**
   * @brief
   * Allocate memory aligned to a boundary
   * @tparam T Allocator type
   * @tparam A alignment boundary, min 4 for 32 bit, 8 for 64 bit systems.
   */
  template<class T, std::size_t A>
  struct aligned_allocator {
      static_assert(!(A & (A - 1)), "A should be a power of two.");
      // Align on >= 4 byte for 32bit, 8 byte for 64bit
      const unsigned al = A < 2 * sizeof(void *) ? 2 * sizeof(void *) : A;

      using value_type = T;
      using pointer = T *;
      using const_pointer = const T *;
      using reference = T &;
      using const_reference = const T &;
      using size_type = std::size_t;
      using difference_type = std::ptrdiff_t;

      aligned_allocator() = default;
      aligned_allocator(const aligned_allocator &) = default;
      template<class U>
      aligned_allocator(const aligned_allocator<U, A> &) {}
      aligned_allocator &operator=(const aligned_allocator &) = delete;
      ~aligned_allocator() = default;

      template<class U>
      struct rebind {
          typedef aligned_allocator<U, A> other;
      };

      constexpr std::size_t max_size() const {
          return std::numeric_limits<size_t>::max() / sizeof(T);
      }

      // stateless allocator
      bool operator!=(const aligned_allocator &) const { return false; }
      bool operator==(const aligned_allocator &) const { return true; }

      T *allocate(std::size_t n) const {
          if (n == 0) return nullptr;

          // aligned_alloc needs a multiple of al
          if (n % al) n += al - (n % al);
          if (n > max_size()) throw std::length_error("aligned_allocator<T,A>::allocate() - Integer overflow.");

          void *p;
          int err = posix_memalign(&p, al, n * sizeof(T)); // HHPC, MARCC likes this better than aligned_alloc
          if (p == nullptr || err != 0) throw std::bad_alloc();

          return static_cast<T *>(p);
      }

      void deallocate(T *p, std::size_t n) const {
          (void) n; // get rid of compiler warning
          free(p);
      }
  };


  template<typename T, unsigned N>
  struct SIMD {
      static_assert(std::is_same<T, char>::value || std::is_same<T, int16_t>::value, "Invalid T in SIMD<T,N>");

      using native_t = T;

      #ifdef VA_SIMD_USE_SSE
      using simd_t = __m128i;
      #else
      using simd_t = __m256i;
      #endif

      static constexpr unsigned length = N;
      static constexpr unsigned size = sizeof(native_t) * N;

      SIMD() = default;
      SIMD(const native_t o) {
          *this = o;
      }
      SIMD(const simd_t &o) {
          v = o;
      }


      __RG_STRONG_INLINE__ SIMD<T, N> operator==(const SIMD<T, N> &o) const;
      __RG_STRONG_INLINE__ SIMD<T, N> &operator=(const SIMD<T, N>::native_t o);
      __RG_STRONG_INLINE__ SIMD<T, N> operator+(const SIMD<T, N> &o) const;
      __RG_STRONG_INLINE__ SIMD<T, N> operator-(const SIMD<T, N> &o) const;
      __RG_STRONG_INLINE__ SIMD<T, N> operator>(const SIMD<T, N> &o) const;
      __RG_STRONG_INLINE__ SIMD<T, N> operator<(const SIMD<T, N> &o) const;
      __RG_STRONG_INLINE__ SIMD<T, N> operator^(const SIMD<T, N> &o) const;
      __RG_STRONG_INLINE__ SIMD<T, N> operator&(const SIMD<T, N> &o) const;
      __RG_STRONG_INLINE__ SIMD<T, N> operator|(const SIMD<T, N> &o) const;
      __RG_STRONG_INLINE__ SIMD<T, N> and_not(const SIMD<T, N> &o) const;
      __RG_STRONG_INLINE__ bool any() const;

      __RG_STRONG_INLINE__
      native_t &operator[](const int i) {
          return reinterpret_cast<native_t *>(&v)[i];
      };

      __RG_STRONG_INLINE__
      SIMD<T, N> operator!() const {
          // XOR with all ones
          return v ^ (v == v);
      };

      __RG_STRONG_INLINE__
      SIMD<T, N> &operator=(const SIMD<T, N> &o) {
          v = o.v;
          return *this;
      };

      __RG_STRONG_INLINE__
      SIMD<T, N> operator>=(const SIMD<T, N> &o) const {
          return !(*this < o);
      };

      __RG_STRONG_INLINE__
      SIMD<T, N> operator<=(const SIMD<T, N> &o) const {
          return !(*this > o);
      };

      __RG_STRONG_INLINE__
      SIMD<T, N> operator!=(const SIMD<T, N> &o) const {
          return !(*this == o);
      };

      __RG_STRONG_INLINE__
      explicit operator bool() const {
          return any();
      };

      simd_t v;
  };

  // SSE2
  using int8x16 = SIMD<char, 16>;
  using int16x8 = SIMD<int16_t, 8>;
  // AVX2
  using int8x32 = SIMD<char, 32>;
  using int16x16 = SIMD<int16_t, 16>;

  using int8_fast = SIMD<char, VA_MAX_INT8>;
  using int16_fast = SIMD<int16_t, VA_MAX_INT16>;

  /**
   * std::vector with an aligned allocator
   */
  template<typename T>
  using SIMDVector = std::vector<T, aligned_allocator<T, T::size>>;

  /************************************ 128b ************************************/

  #ifdef VA_SIMD_USE_SSE

  template<>
  int8x16 int8x16::operator==(const int8x16 &o) const {
      return _mm_cmpeq_epi8(v, o.v);
  }
  template<>
  int8x16 int8x16::operator^(const int8x16 &o) const {
      // XOR with all ones
      return _mm_xor_si128(v, o.v);
  }
  template<>
  int8x16 &int8x16::operator=(const int8x16::native_t o) {
      v = _mm_set1_epi8(o);
      return *this;
  }
  template<>
  int8x16 int8x16::operator+(const int8x16 &o) const {
      return _mm_adds_epi8(v, o.v);
  }
  template<>
  int8x16 int8x16::operator-(const int8x16 &o) const {
      return _mm_subs_epi8(v, o.v);
  }
  template<>
  int8x16 int8x16::operator>(const int8x16 &o) const {
      return _mm_cmpgt_epi8(v, o.v);
  }
  template<>
  int8x16 int8x16::operator<(const int8x16 &o) const {
      return _mm_cmplt_epi8(v, o.v);
  }
  template<>
  int8x16 int8x16::operator&(const int8x16 &o) const {
      return _mm_and_si128(v, o.v);
  }
  template<>
  int8x16 int8x16::operator|(const int8x16 &o) const {
      return _mm_or_si128(v, o.v);
  }
  template<>
  bool int8x16::any() const {
      return _mm_movemask_epi8(v);
  }
  template<>
  int8x16 int8x16::and_not(const int8x16 &o) const {
      return _mm_andnot_si128(o.v, v);
  }

  __RG_STRONG_INLINE__
  int8x16 max(const int8x16 &a, const int8x16 &b) {
      return _mm_max_epi8(a.v, b.v);
  }
  __RG_STRONG_INLINE__
  int8x16 blend(const int8x16 &mask, const int8x16 &t, const int8x16 &f) {
      return _mm_blendv_epi8(f.v, t.v, mask.v);
  }


  template<>
  int16x8 int16x8::operator==(const int16x8 &o) const {
      return _mm_cmpeq_epi16(v, o.v);
  }
  template<>
  int16x8 int16x8::operator^(const int16x8 &o) const {
      return _mm_xor_si128(v, o.v);
  }
  template<>
  int16x8 &int16x8::operator=(const int16x8::native_t o) {
      v = _mm_set1_epi16(o);
      return *this;
  }
  template<>
  int16x8 int16x8::operator+(const int16x8 &o) const {
      return _mm_adds_epi16(v, o.v);
  }
  template<>
  int16x8 int16x8::operator-(const int16x8 &o) const {
      return _mm_subs_epi16(v, o.v);
  }
  template<>
  int16x8 int16x8::operator>(const int16x8 &o) const {
      return _mm_cmpgt_epi16(v, o.v);
  }
  template<>
  int16x8 int16x8::operator<(const int16x8 &o) const {
      return _mm_cmplt_epi16(v, o.v);
  }
  template<>
  int16x8 int16x8::operator&(const int16x8 &o) const {
      return _mm_and_si128(v, o.v);
  }
  template<>
  int16x8 int16x8::operator|(const int16x8 &o) const {
      return _mm_or_si128(v, o.v);
  }
  template<>
  bool int16x8::any() const {
      return _mm_movemask_epi8(v);
  }
  template<>
  int16x8 int16x8::and_not(const int16x8 &o) const {
      return _mm_andnot_si128(o.v, v);
  }

  __RG_STRONG_INLINE__
  int16x8 max(const int16x8 &a, const int16x8 &b) {
      return _mm_max_epi16(a.v, b.v);
  }
  __RG_STRONG_INLINE__
  int16x8 blend(const int16x8 &mask, const int16x8 &t, const int16x8 &f) {
      return _mm_blendv_epi8(f.v, t.v, mask.v);
  }

  #endif

  /************************************ 256b ************************************/

  #ifdef VA_SIMD_USE_AVX2

  template<> int8x32 int8x32::operator==(const int8x32 &o) const {
      return _mm256_cmpeq_epi8(v, o.v);
  }
  template<> int8x32 int8x32::operator^(const int8x32 &o) const {
      return _mm256_xor_si256(v, o.v);
  }
  template<> int8x32 &int8x32::operator=(const int8x32::native_t o) {
      v = _mm256_set1_epi8(o);
      return *this;
  }
  template<> int8x32 int8x32::operator+(const int8x32 &o) const {
      return _mm256_adds_epi8(v, o.v);
  }
  template<> int8x32 int8x32::operator-(const int8x32 &o) const {
      return _mm256_subs_epi8(v, o.v);
  }
  template<> int8x32 int8x32::operator>(const int8x32 &o) const {
      return _mm256_cmpgt_epi8(v, o.v);
  }
  template<> int8x32 int8x32::operator<(const int8x32 &o) const {
      return _mm256_cmpgt_epi8(o.v, v);
  }
  template<> int8x32 int8x32::operator&(const int8x32 &o) const {
      return _mm256_and_si256(v, o.v);
  }
  template<> int8x32 int8x32::operator|(const int8x32 &o) const {
      return _mm256_or_si256(v, o.v);
  }
  template<> bool int8x32::any() const {
      return _mm256_movemask_epi8(v);
  }
    template <>
  int8x32 int8x32::and_not(const int8x32 &o) const {
      return _mm256_andnot_si256(o.v, v);
  }
  __RG_STRONG_INLINE__
  int8x32 max(const int8x32 &a, const int8x32 &b) {
      return _mm256_max_epi8(a.v, b.v);
  }
  __RG_STRONG_INLINE__
  int8x32 blend(const int8x32 &mask, const int8x32 &t, const int8x32 &f) {
      return _mm256_blendv_epi8(f.v, t.v, mask.v);
  }


  template<> int16x16 int16x16::operator==(const int16x16 &o) const {
      return _mm256_cmpeq_epi16(v, o.v);
  }
  template<> int16x16 int16x16::operator^(const int16x16 &o) const {
      return _mm256_xor_si256(v, o.v);
  }
  template<> int16x16 &int16x16::operator=(const int16x16::native_t o) {
      v = _mm256_set1_epi16(o);
      return *this;
  }
  template<> int16x16 int16x16::operator+(const int16x16 &o) const {
      return _mm256_adds_epi16(v, o.v);
  }
  template<> int16x16 int16x16::operator-(const int16x16 &o) const {
      return _mm256_subs_epi16(v, o.v);
  }
  template<> int16x16 int16x16::operator>(const int16x16 &o) const {
      return _mm256_cmpgt_epi16(v, o.v);
  }
  template<> int16x16 int16x16::operator<(const int16x16 &o) const {
      return _mm256_cmpgt_epi16(o.v, v);
  }
  template<> int16x16 int16x16::operator&(const int16x16 &o) const {
      return _mm256_and_si256(v, o.v);
  }
  template<> int16x16 int16x16::operator|(const int16x16 &o) const {
      return _mm256_or_si256(v, o.v);
  }
  template<> bool int16x16::any() const {
      return _mm256_movemask_epi8(v);
  }
  template <>
  int16x16 int16x16::and_not(const int16x16 &o) const {
      return _mm256_andnot_si256(o.v, v);
  }
  __RG_STRONG_INLINE__
  int16x16 max(const int16x16 &a, const int16x16 &b) {
      return _mm256_max_epi16(a.v, b.v);
  }
  __RG_STRONG_INLINE__
  int16x16 blend(const int16x16 &mask, const int16x16 &t, const int16x16 &f) {
      return _mm256_blendv_epi8(f.v, t.v, mask.v);
  }

  #endif

}


#endif //VARGAS_SIMD_H
