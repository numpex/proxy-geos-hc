#include "commonConfig.hpp"

#ifdef SEM_USE_RAJA

#include "solverRaja.hpp"

#elif defined SEM_USE_OMP

#include "solverOMP.hpp"

#elif defined SEM_USE_KOKKOS

#include "solverKokkos.hpp"

#else

#include "solverSequential.hpp"

#endif

#ifdef SEM_USE_CALIPER

#include <caliper/cali.h>
#include <sys/time.h>
#include <string>
#include <iostream>

/// Mark a function or scope for timing with a given name
#define SEM_CALIPER_MARK_SCOPE( name ) cali::Function __cali_ann ## __LINE__((name))

/// Mark a function for timing using a compiler-provided name
#define SEM_CALIPER_MARK_FUNCTION cali::Function __cali_ann ## __func__( timingHelpers::stripPF( __PRETTY_FUNCTION__ ).c_str())

/// Mark the beginning of timed statement group
#define SEM_CALIPER_MARK_BEGIN( name ) CALI_MARK_BEGIN((name))

/// Mark the end of timed statements group
#define SEM_CALIPER_MARK_END( name ) CALI_MARK_END((name))

/// Mark the beginning of function, only useful when you don't want to or can't mark the whole function.
#define SEM_CALIPER_MARK_FUNCTION_BEGIN CALI_MARK_FUNCTION_BEGIN

/// Mark the end of function, only useful when you don't want to or can't mark the whole function.
#define SEM_CALIPER_MARK_FUNCTION_END CALI_MARK_FUNCTION_END

#else // SEM_USE_CALIPER

/// @cond DO_NOT_DOCUMENT
#define SEM_CALIPER_MARK_SCOPE( name )
#define SEM_CALIPER_MARK_FUNCTION

#define SEM_CALIPER_MARK_BEGIN( name )
#define SEM_CALIPER_MARK_END( name )

#define SEM_CALIPER_MARK_FUNCTION_BEGIN
#define SEM_CALIPER_MARK_FUNCTION_END
/// @endcond

#endif // SEM_USE_CALIPER
