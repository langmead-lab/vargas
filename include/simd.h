//
// Created by gaddra on 1/6/17.
//

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
   * @tparam T
   * @tparam A alignment
   */
  template<class T, std::size_t A>
  struct aligned_allocator {
      static_assert(!(A & (A - 1)), "A should be a power of two.");
      // Align on >= 4 byte for 32bit, 8 byte for 64bit
      const std::size_t al = A < 2 * sizeof(void *) ? 2 * sizeof(void *) : A;

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

      bool operator!=(const aligned_allocator &) const { return false; }
      bool operator==(const aligned_allocator &) const { return true; }

      T *allocate(std::size_t n) const {
          if (n == 0) return nullptr;

          // aligned_alloc needs a multiple of al
          if (n % al) n += al - (n % al);
          if (n > max_size()) throw std::length_error("aligned_allocator<T,A>::allocate() - Integer overflow.");

          void *p = aligned_alloc(al, n * sizeof(T));
          if (p == nullptr) throw std::bad_alloc();

          return static_cast<T *>(p);
      }

      void deallocate(T *p, std::size_t n) const {
          (void) n; // get rid of compiler warning
          free(p);
      }
  };


  template<typename T, size_t N>
  struct SIMD {
      static_assert(N == 8 || N == 16 || N == 32, "Invalid N in SIMD<T,N>");
      static_assert(std::is_same<T, int8_t>::value || std::is_same<T, int16_t>::value, "Invalid T in SIMD<T,N>");

      using native_t = T;
      using simd_t = typename std::conditional<std::is_same<int8_t, T>::value,
                                               typename std::conditional<N == 16, __m128i,
                                                                         typename std::conditional<N == 32,
                                                                                                   __m256i,
                                                                                                   __m512i>::type>::type,
                                               typename std::conditional<N == 8, __m128i,
                                                                         typename std::conditional<N == 16,
                                                                                                   __m256i,
                                                                                                   __m512i>::type>::type>::type;

      static constexpr size_t length = N;
      static constexpr size_t size = sizeof(native_t) * N;

      SIMD() = default;
      SIMD(const native_t o) {
          *this = o;
      }
      SIMD(const simd_t &o) {
          v = o;
      }


      __RG_STRONG_INLINE__ SIMD<T, N> operator==(const SIMD<T, N> &o) const;
      __RG_STRONG_INLINE__ SIMD<T, N> operator!() const;
      __RG_STRONG_INLINE__ SIMD<T, N> &operator=(const SIMD<T, N>::native_t o);
      __RG_STRONG_INLINE__ SIMD<T, N> operator+(const SIMD<T, N> &o) const;
      __RG_STRONG_INLINE__ SIMD<T, N> operator-(const SIMD<T, N> &o) const;
      __RG_STRONG_INLINE__ SIMD<T, N> operator>(const SIMD<T, N> &o) const;
      __RG_STRONG_INLINE__ SIMD<T, N> operator<(const SIMD<T, N> &o) const;
      __RG_STRONG_INLINE__ SIMD<T, N> operator&(const SIMD<T, N> &o) const;
      __RG_STRONG_INLINE__ SIMD<T, N> operator|(const SIMD<T, N> &o) const;
      __RG_STRONG_INLINE__ native_t at(const int i) const;
      __RG_STRONG_INLINE__ void insert(const int i, const native_t e);
      __RG_STRONG_INLINE__ SIMD<T, N> &operator=(const SIMD<T, N> &o) {
          v = o.v;
          return *this;
      };
      __RG_STRONG_INLINE__ SIMD<T, N> operator>=(const SIMD<T, N> &o) const {
          return !(*this < o);
      };
      __RG_STRONG_INLINE__ SIMD<T, N> operator<=(const SIMD<T, N> &o) const {
          return !(*this > o);
      };
      __RG_STRONG_INLINE__ SIMD<T, N> operator!=(const SIMD<T, N> &o) const {
          return !(*this == o);
      };
      __RG_STRONG_INLINE__ bool any() const;
      __RG_STRONG_INLINE__ SIMD<T, N> and_not(const SIMD<T, N> &o) const;

      simd_t v;
  };

  template<typename T>
  __RG_STRONG_INLINE__
  typename T::native_t extract(const size_t i, const T &v) {
      return reinterpret_cast<const typename T::native_t *>(&v.v)[i];
  }

  template<typename T>
  __RG_STRONG_INLINE__
  void insert(typename T::native_t elem, const size_t i, T &v) {
      reinterpret_cast<typename T::native_t *>(&v.v)[i] = elem;
  }

  // SSE2
  using int8x16 = SIMD<int8_t, 16>;
  using int16x8 = SIMD<int16_t, 8>;
  // AVX2
  using int8x32 = SIMD<int8_t, 32>;
  using int16x16 = SIMD<int16_t, 16>;

  template<typename T>
  using SIMDVector = std::vector<T, aligned_allocator<T, T::size>>;

  /************************************ 128b ************************************/

  #ifdef VA_SIMD_USE_SSE

  template<>
  int8x16 int8x16::operator==(const int8x16 &o) const {
      return _mm_cmpeq_epi8(v, o.v);
  }
  template<>
  int8x16 int8x16::operator!() const {
      // XOR with all ones
      return _mm_xor_si128(v, _mm_cmpeq_epi8(v, v));
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
  template<>
  typename int8x16::native_t int8x16::at(const int i) const {
      switch (i) {
          default:
          case 0:
              return _mm_extract_epi8(v, 0);
          case 1:
              return _mm_extract_epi8(v, 1);
          case 2:
              return _mm_extract_epi8(v, 2);
          case 3:
              return _mm_extract_epi8(v, 3);
          case 4:
              return _mm_extract_epi8(v, 4);
          case 5:
              return _mm_extract_epi8(v, 5);
          case 6:
              return _mm_extract_epi8(v, 6);
          case 7:
              return _mm_extract_epi8(v, 7);
          case 8:
              return _mm_extract_epi8(v, 8);
          case 9:
              return _mm_extract_epi8(v, 9);
          case 10:
              return _mm_extract_epi8(v, 10);
          case 11:
              return _mm_extract_epi8(v, 11);
          case 12:
              return _mm_extract_epi8(v, 12);
          case 13:
              return _mm_extract_epi8(v, 13);
          case 14:
              return _mm_extract_epi8(v, 14);
          case 15:
              return _mm_extract_epi8(v, 15);
      }
  }
  template<>
  void int8x16::insert(const int i, typename int8x16::native_t e) {
      switch (i) {
          default:
          case 0:
              v = _mm_insert_epi8(v, e, 0);
              return;
          case 1:
              v = _mm_insert_epi8(v, e, 1);
              return;
          case 2:
              v = _mm_insert_epi8(v, e, 2);
              return;
          case 3:
              v = _mm_insert_epi8(v, e, 3);
              return;
          case 4:
              v = _mm_insert_epi8(v, e, 4);
              return;
          case 5:
              v = _mm_insert_epi8(v, e, 5);
              return;
          case 6:
              v = _mm_insert_epi8(v, e, 6);
              return;
          case 7:
              v = _mm_insert_epi8(v, e, 7);
              return;
          case 8:
              v = _mm_insert_epi8(v, e, 8);
              return;
          case 9:
              v = _mm_insert_epi8(v, e, 9);
              return;
          case 10:
              v = _mm_insert_epi8(v, e, 10);
              return;
          case 11:
              v = _mm_insert_epi8(v, e, 11);
              return;
          case 12:
              v = _mm_insert_epi8(v, e, 12);
              return;
          case 13:
              v = _mm_insert_epi8(v, e, 13);
              return;
          case 14:
              v = _mm_insert_epi8(v, e, 14);
              return;
          case 15:
              v = _mm_insert_epi8(v, e, 15);
              return;
      }
  }

  __RG_STRONG_INLINE__ int8x16 max(const int8x16 &a, const int8x16 &b) {
      return _mm_max_epi8(a.v, b.v);
  }
  __RG_STRONG_INLINE__ int8x16 blend(const int8x16 &mask, const int8x16 &t, const int8x16 &f) {
      return _mm_blendv_epi8(f.v, t.v, mask.v);
  }


  template<>
  int16x8 int16x8::operator==(const int16x8 &o) const {
      return _mm_cmpeq_epi16(v, o.v);
  }
  template<>
  int16x8 int16x8::operator!() const {
      return _mm_xor_si128(v, _mm_cmpeq_epi16(v, v));
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
  template<>
  typename int16x8::native_t int16x8::at(const int i) const {
      switch (i) {
          default:
          case 0:
              return _mm_extract_epi16(v, 0);
          case 1:
              return _mm_extract_epi16(v, 1);
          case 2:
              return _mm_extract_epi16(v, 2);
          case 3:
              return _mm_extract_epi16(v, 3);
          case 4:
              return _mm_extract_epi16(v, 4);
          case 5:
              return _mm_extract_epi16(v, 5);
          case 6:
              return _mm_extract_epi16(v, 6);
          case 7:
              return _mm_extract_epi16(v, 7);
      }
  }
  template<>
  void int16x8::insert(const int i, typename int16x8::native_t e) {
      switch (i) {
          default:
          case 0:
              v = _mm_insert_epi16(v, e, 0);
              return;
          case 1:
              v = _mm_insert_epi16(v, e, 1);
              return;
          case 2:
              v = _mm_insert_epi16(v, e, 2);
              return;
          case 3:
              v = _mm_insert_epi16(v, e, 3);
              return;
          case 4:
              v = _mm_insert_epi16(v, e, 4);
              return;
          case 5:
              v = _mm_insert_epi16(v, e, 5);
              return;
          case 6:
              v = _mm_insert_epi16(v, e, 6);
              return;
          case 7:
              v = _mm_insert_epi16(v, e, 7);
              return;
      }
  }

  __RG_STRONG_INLINE__ int16x8 max(const int16x8 &a, const int16x8 &b) {
      return _mm_max_epi16(a.v, b.v);
  }
  __RG_STRONG_INLINE__ int16x8 blend(const int16x8 &mask, const int16x8 &t, const int16x8 &f) {
      return _mm_blendv_epi8(f.v, t.v, mask.v);
  }

  #endif

  /************************************ 256b ************************************/

  #ifdef VA_SIMD_USE_AVX2

  template<> int8x32 int8x32::operator==(const int8x32 &o) const {
      return _mm256_cmpeq_epi8(v, o.v);
  }
  template<> int8x32 int8x32::operator!() const {
      return _mm256_xor_si256(v, _mm256_cmpeq_epi8(v,v));
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
  template<> typename int8x32::native_t int8x32::at(const int i) const {
      switch (i) {
          default:
          case 0: return _mm256_extract_epi8(v, 0);
          case 1: return _mm256_extract_epi8(v, 1);
          case 2: return _mm256_extract_epi8(v, 2);
          case 3: return _mm256_extract_epi8(v, 3);
          case 4: return _mm256_extract_epi8(v, 4);
          case 5: return _mm256_extract_epi8(v, 5);
          case 6: return _mm256_extract_epi8(v, 6);
          case 7: return _mm256_extract_epi8(v, 7);
          case 8: return _mm256_extract_epi8(v, 8);
          case 9: return _mm256_extract_epi8(v, 9);
          case 10: return _mm256_extract_epi8(v, 10);
          case 11: return _mm256_extract_epi8(v, 11);
          case 12: return _mm256_extract_epi8(v, 12);
          case 13: return _mm256_extract_epi8(v, 13);
          case 14: return _mm256_extract_epi8(v, 14);
          case 15: return _mm256_extract_epi8(v, 15);
          case 16: return _mm256_extract_epi8(v, 16);
          case 17: return _mm256_extract_epi8(v, 17);
          case 18: return _mm256_extract_epi8(v, 18);
          case 19: return _mm256_extract_epi8(v, 19);
          case 20: return _mm256_extract_epi8(v, 20);
          case 21: return _mm256_extract_epi8(v, 21);
          case 22: return _mm256_extract_epi8(v, 22);
          case 23: return _mm256_extract_epi8(v, 23);
          case 24: return _mm256_extract_epi8(v, 24);
          case 25: return _mm256_extract_epi8(v, 25);
          case 26: return _mm256_extract_epi8(v, 26);
          case 27: return _mm256_extract_epi8(v, 27);
          case 28: return _mm256_extract_epi8(v, 28);
          case 29: return _mm256_extract_epi8(v, 29);
          case 30: return _mm256_extract_epi8(v, 30);
          case 31: return _mm256_extract_epi8(v, 31);
      }
  }
    template<>
  void int8x32::insert(const int i, typename int8x32::native_t e) {
      switch(i) {
          default:
          case 0: v = _mm256_insert_epi8(v, e, 0); return;
          case 1: v = _mm256_insert_epi8(v, e, 1); return;
          case 2: v = _mm256_insert_epi8(v, e, 2); return;
          case 3: v = _mm256_insert_epi8(v, e, 3); return;
          case 4: v = _mm256_insert_epi8(v, e, 4); return;
          case 5: v = _mm256_insert_epi8(v, e, 5); return;
          case 6: v = _mm256_insert_epi8(v, e, 6); return;
          case 7: v = _mm256_insert_epi8(v, e, 7); return;
          case 8: v = _mm256_insert_epi8(v, e, 8); return;
          case 9: v = _mm256_insert_epi8(v, e, 9); return;
          case 10: v = _mm256_insert_epi8(v, e, 10); return;
          case 11: v = _mm256_insert_epi8(v, e, 11); return;
          case 12: v = _mm256_insert_epi8(v, e, 12); return;
          case 13: v = _mm256_insert_epi8(v, e, 13); return;
          case 14: v = _mm256_insert_epi8(v, e, 14); return;
          case 15: v = _mm256_insert_epi8(v, e, 15); return;
          case 16: v = _mm256_insert_epi8(v, e, 16); return;
          case 17: v = _mm256_insert_epi8(v, e, 17); return;
          case 18: v = _mm256_insert_epi8(v, e, 18); return;
          case 19: v = _mm256_insert_epi8(v, e, 19); return;
          case 20: v = _mm256_insert_epi8(v, e, 20); return;
          case 21: v = _mm256_insert_epi8(v, e, 21); return;
          case 22: v = _mm256_insert_epi8(v, e, 22); return;
          case 23: v = _mm256_insert_epi8(v, e, 23); return;
          case 24: v = _mm256_insert_epi8(v, e, 24); return;
          case 25: v = _mm256_insert_epi8(v, e, 25); return;
          case 26: v = _mm256_insert_epi8(v, e, 26); return;
          case 27: v = _mm256_insert_epi8(v, e, 27); return;
          case 28: v = _mm256_insert_epi8(v, e, 28); return;
          case 29: v = _mm256_insert_epi8(v, e, 29); return;
          case 30: v = _mm256_insert_epi8(v, e, 30); return;
          case 31: v = _mm256_insert_epi8(v, e, 31); return;
      }
  }
    template <>
  int8x32 int8x32::and_not(const int8x32 &o) const {
      return _mm256_andnot_si256(o.v, v);
  }
  __RG_STRONG_INLINE__ int8x32 max(const int8x32 &a, const int8x32 &b) {
      return _mm256_max_epi8(a.v, b.v);
  }
  __RG_STRONG_INLINE__ int8x32 blend(const int8x32 &mask, const int8x32 &t, const int8x32 &f) {
      return _mm256_blendv_epi8(f.v, t.v, mask.v);
  }


  template<> int16x16 int16x16::operator==(const int16x16 &o) const {
      return _mm256_cmpeq_epi16(v, o.v);
  }
  template<> int16x16 int16x16::operator!() const {
      return _mm256_xor_si256(v, _mm256_cmpeq_epi16(v,v));
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
  template<> typename int16x16::native_t int16x16::at(const int i) const {
      switch (i) {
          default:
          case 0: return _mm256_extract_epi16(v, 0);
          case 1: return _mm256_extract_epi16(v, 1);
          case 2: return _mm256_extract_epi16(v, 2);
          case 3: return _mm256_extract_epi16(v, 3);
          case 4: return _mm256_extract_epi16(v, 4);
          case 5: return _mm256_extract_epi16(v, 5);
          case 6: return _mm256_extract_epi16(v, 6);
          case 7: return _mm256_extract_epi16(v, 7);
          case 8: return _mm256_extract_epi16(v, 8);
          case 9: return _mm256_extract_epi16(v, 9);
          case 10: return _mm256_extract_epi16(v, 10);
          case 11: return _mm256_extract_epi16(v, 11);
          case 12: return _mm256_extract_epi16(v, 12);
          case 13: return _mm256_extract_epi16(v, 13);
          case 14: return _mm256_extract_epi16(v, 14);
          case 15: return _mm256_extract_epi16(v, 15);

      }
  }
    template<>
  void int16x16::insert(const int i, typename int16x16::native_t e) {
      switch(i) {
          default:
          case 0: v = _mm256_insert_epi16(v, e, 0); return;
          case 1: v = _mm256_insert_epi16(v, e, 1); return;
          case 2: v = _mm256_insert_epi16(v, e, 2); return;
          case 3: v = _mm256_insert_epi16(v, e, 3); return;
          case 4: v = _mm256_insert_epi16(v, e, 4); return;
          case 5: v = _mm256_insert_epi16(v, e, 5); return;
          case 6: v = _mm256_insert_epi16(v, e, 6); return;
          case 7: v = _mm256_insert_epi16(v, e, 7); return;
          case 8: v = _mm256_insert_epi16(v, e, 8); return;
          case 9: v = _mm256_insert_epi16(v, e, 9); return;
          case 10: v = _mm256_insert_epi16(v, e, 10); return;
          case 11: v = _mm256_insert_epi16(v, e, 11); return;
          case 12: v = _mm256_insert_epi16(v, e, 12); return;
          case 13: v = _mm256_insert_epi16(v, e, 13); return;
          case 14: v = _mm256_insert_epi16(v, e, 14); return;
          case 15: v = _mm256_insert_epi16(v, e, 15); return;
      }
  }
  template <>
  int16x16 int16x16::and_not(const int16x16 &o) const {
      return _mm256_andnot_si256(o.v, v);
  }
  __RG_STRONG_INLINE__ int16x16 max(const int16x16 &a, const int16x16 &b) {
      return _mm256_max_epi16(a.v, b.v);
  }
  __RG_STRONG_INLINE__ int16x16 blend(const int16x16 &mask, const int16x16 &t, const int16x16 &f) {
      return _mm256_blendv_epi8(f.v, t.v, mask.v);
  }

  #endif

}


#endif //VARGAS_SIMD_H
