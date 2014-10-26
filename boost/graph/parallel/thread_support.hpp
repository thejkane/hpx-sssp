// Copyright 2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nicholas Edmonds
//           Andrew Lumsdaine

#ifndef BOOST_PARALLEL_THREAD_SUPPORT_HPP
#define BOOST_PARALLEL_THREAD_SUPPORT_HPP

#if defined(__APPLE__)
#  include <libkern/OSAtomic.h>
#endif

#include <boost/utility/enable_if.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/bool.hpp>
#include <stdint.h>

namespace boost { namespace parallel {

#if defined(__MACH__) && defined(__APPLE__)

  template <typename T, typename = void> struct atomics_supported: boost::mpl::false_ {};

  // Need a more sophisticated for detecting available atomics on a
  // variety of compilers, but for now just set these appropriately

  // #ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_1
  template <typename T> struct atomics_supported<T, typename boost::enable_if_c<sizeof(T) == 1>::type>: boost::mpl::true_ {};

  // #ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_2
  template <typename T> struct atomics_supported<T, typename boost::enable_if_c<sizeof(T) == 2>::type>: boost::mpl::true_ {};

  // #ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
  template <typename T> struct atomics_supported<T, typename boost::enable_if_c<sizeof(T) == 4>::type>: boost::mpl::true_ {};

  // #ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_8
  template <typename T> struct atomics_supported<T, typename boost::enable_if_c<sizeof(T) == 8>::type>: boost::mpl::true_ {};

  // #ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_16
  // template <typename T> struct atomics_supported<T, typename boost::enable_if_c<sizeof(T) == 16>::type>: boost::mpl::true_ {};

  //
  // TODO: Need atomics for unsigned integral types and floats!
  //


//   inline bool bool_compare_and_swap(volatile int* pa, int oldval, int newval)
//   { return OSAtomicCompareAndSwapInt(oldval, newval, pa); }

   inline bool bool_compare_and_swap(volatile long* pa, long oldval, long newval)
   { return OSAtomicCompareAndSwapLong(oldval, newval, pa); }

   inline bool bool_compare_and_swap(volatile unsigned long* pa, unsigned long oldval, unsigned long newval)
   { return OSAtomicCompareAndSwapLong(static_cast<long>(oldval), static_cast<long>(newval), *(volatile long**)&pa); }

  inline bool bool_compare_and_swap(void* volatile *pa, void* oldval, void* newval)
  { return OSAtomicCompareAndSwapPtr(oldval, newval, pa); }

  inline bool bool_compare_and_swap(int32_t* pa, int32_t oldval, int32_t newval)
  { return OSAtomicCompareAndSwap32(oldval, newval, *(volatile int32_t**)&pa); }

  inline bool bool_compare_and_swap(unsigned int* pa, unsigned int oldval, unsigned int newval)
  { return OSAtomicCompareAndSwap32(static_cast<int32_t>(oldval), static_cast<int32_t>(newval),
				    *(volatile int32_t**)&pa); }

  inline bool bool_compare_and_swap(int64_t* pa, int64_t oldval, int64_t newval)
  { return OSAtomicCompareAndSwap64(oldval, newval, *(volatile int64_t**)&pa); }

  inline bool bool_compare_and_swap(unsigned long long* pa, unsigned long long oldval, 
				    unsigned long long newval)
  { return OSAtomicCompareAndSwap64(static_cast<int64_t>(oldval), static_cast<int64_t>(newval), 
				    *(volatile int64_t**)&pa); }

  inline int32_t val_compare_and_swap(int32_t* pa, int32_t oldval, int32_t newval)
  { 
    if (OSAtomicCompareAndSwap32(oldval, newval, pa))
      return oldval;
    else
      return *pa;
  }

  inline unsigned int val_compare_and_swap(unsigned int* pa, unsigned int oldval, unsigned int newval)
  { 
    if (OSAtomicCompareAndSwap32(static_cast<int32_t>(oldval), static_cast<int32_t>(newval), 
				 *(volatile int32_t**)&pa))
      return oldval;
    else
      return *pa;
  }

  inline int64_t val_compare_and_swap(int64_t* pa, int64_t oldval, int64_t newval)
  { 
    if (OSAtomicCompareAndSwap64(oldval, newval, pa))
      return oldval;
    else
      return *pa;
  }

  inline int32_t fetch_and_add(volatile int32_t* pa, int32_t val)
  { return OSAtomicAdd32(val, pa); }

  inline int64_t fetch_and_add(volatile int64_t* pa, int64_t val)
  { return OSAtomicAdd64(val, pa); }

#else

  template <typename T, typename = void> struct atomics_supported: boost::mpl::false_ {};

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_1
  template <typename T> struct atomics_supported<T, typename boost::enable_if_c<sizeof(T) == 1>::type>: boost::mpl::true_ {};
#endif

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_2
  template <typename T> struct atomics_supported<T, typename boost::enable_if_c<sizeof(T) == 2>::type>: boost::mpl::true_ {};
#endif

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
  template <typename T> struct atomics_supported<T, typename boost::enable_if_c<sizeof(T) == 4>::type>: boost::mpl::true_ {};
#endif

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_8
  template <typename T> struct atomics_supported<T, typename boost::enable_if_c<sizeof(T) == 8>::type>: boost::mpl::true_ {};
#endif

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_16
  template <typename T> struct atomics_supported<T, typename boost::enable_if_c<sizeof(T) == 16>::type>: boost::mpl::true_ {};
#endif

  template <typename T>
  inline bool bool_compare_and_swap(T* pa, T oldval, T newval)
  { BOOST_MPL_ASSERT((atomics_supported<T>)); return __sync_bool_compare_and_swap(pa, oldval, newval); }

  template <typename T>
  inline T val_compare_and_swap(T* pa, T oldval, T newval)
  { BOOST_MPL_ASSERT((atomics_supported<T>)); return __sync_val_compare_and_swap(pa, oldval, newval); }

  template <typename T>
  inline T fetch_and_add(T* pa, T val)
  { BOOST_MPL_ASSERT((atomics_supported<T>)); return __sync_fetch_and_add(pa, val); }

#endif

} } // end namespace boost::parallel

#endif // BOOST_PARALLEL_THREAD_SUPPORT_HPP
