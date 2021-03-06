//*****************************************************************************
// config.h           Compiler language support flags
//
// This file (if named config.h) was generated automatically when running cmake
// from the file config.cmake.h.in.

// Macro for declaring aligned variables.
/* #undef BZ_HAVE_ALIGNEMENT_DIRECTIVE_WINDOWS_STYLE */
#ifdef BZ_HAVE_ALIGNEMENT_DIRECTIVE_WINDOWS_STYLE
    #define BZ_ALIGN_VARIABLE(vartype,varname,alignment) __declspec(align(alignment)) vartype varname;
    #undef BZ_HAVE_ALIGNEMENT_DIRECTIVE_WINDOWS_STYLE
#endif

/* #undef BZ_HAVE_ALIGNEMENT_DIRECTIVE_GCC_STYLE */
#ifdef BZ_HAVE_ALIGNEMENT_DIRECTIVE_GCC_STYLE
    #define BZ_ALIGN_VARIABLE(vartype,varname,alignment) vartype __attribute__ ((aligned (alignment))) varname;
    #undef BZ_HAVE_ALIGNEMENT_DIRECTIVE_GCC_STYLE
#endif

#ifndef BZ_ALIGN_VARIABLE
#define BZ_ALIGN_VARIABLE(vartype,varname,alignment) vartype varname;
#endif

// Enable dimensions with > 2^31 elements (NOT IMPLEMENTED).
/* #undef BZ_FULLY64BIT */

// define if bool is a built-in type
#define BZ_HAVE_BOOL

// define if the Boost library is available
/* #undef BZ_HAVE_BOOST */

// Define to 1 if you have the <boost/mpi.hpp> header file.
/* #undef BZ_HAVE_BOOST_MPI */

// define if the Boost::Serialization library is available
/* #undef BZ_HAVE_BOOST_SERIALIZATION */

// define if the compiler has <climits> header
#define BZ_HAVE_CLIMITS

// define if the compiler has complex<T>
#define BZ_HAVE_COMPLEX

// define if the compiler has standard complex<T> functions
#define BZ_HAVE_COMPLEX_FCNS

// define if the compiler has complex math functions
#define BZ_HAVE_COMPLEX_MATH1

// define if the compiler has more complex math functions
/* #undef BZ_HAVE_COMPLEX_MATH2 */

// define if complex math functions are in namespace std
#define BZ_HAVE_COMPLEX_MATH_IN_NAMESPACE_STD

// define if the compiler supports const_cast<>
#define BZ_HAVE_CONST_CAST

// Define to 1 if you have the <cstring> header file.
#define BZ_HAVE_CSTRING

// define if the compiler supports default template parameters
#define BZ_HAVE_DEFAULT_TEMPLATE_PARAMETERS

// Obsolete ?
// Define to 1 if you have the <dlfcn.h> header file.
/* #undef BZ_HAVE_DLFCN_H */

// define if the compiler supports dynamic_cast<>
#define BZ_HAVE_DYNAMIC_CAST

// define if the compiler handle computations inside an enum
#define BZ_HAVE_ENUM_COMPUTATIONS

// define if the compiler handles (int) casts in enum computations
#define BZ_HAVE_ENUM_COMPUTATIONS_WITH_CAST

// define if the compiler supports exceptions
#define BZ_HAVE_EXCEPTIONS

// define if the compiler supports the explicit keyword
#define BZ_HAVE_EXPLICIT

// define if the compiler supports explicit template function qualification
#define BZ_HAVE_EXPLICIT_TEMPLATE_FUNCTION_QUALIFICATION

// define if the compiler recognizes the full specialization syntax
#define BZ_HAVE_FULL_SPECIALIZATION_SYNTAX

// define if the compiler supports function templates with non-type parameters
#define BZ_HAVE_FUNCTION_NONTYPE_PARAMETERS

// define if the compiler supports IEEE math library
#define BZ_HAVE_IEEE_MATH

// Define to 1 if you have the <inttypes.h> header file.
#define BZ_HAVE_INTTYPES_H

// Obsolete ?
// Define to 1 if you have the `m' library (-lm).
#define BZ_HAVE_LIBM

// Define to 1 if you have the `papi' library (-lpapi).
/* #undef BZ_HAVE_LIBPAPI */

// define if the compiler supports member constants
#define BZ_HAVE_MEMBER_CONSTANTS

// define if the compiler supports member templates
#define BZ_HAVE_MEMBER_TEMPLATES

// define if the compiler supports member templates outside the class
// declaration
#define BZ_HAVE_MEMBER_TEMPLATES_OUTSIDE_CLASS

// define if the compiler supports the mutable keyword
#define BZ_HAVE_MUTABLE

// define if the compiler supports the Numerical C Extensions Group restrict
// keyword
/* #undef BZ_HAVE_NCEG_RESTRICT */

// define if the compiler supports the __restrict__ keyword
#define BZ_HAVE_NCEG_RESTRICT_EGCS

// define if the compiler has numeric_limits<T>
#define BZ_HAVE_NUMERIC_LIMITS

// define if the compiler accepts the old for scoping rules
/* #undef BZ_HAVE_OLD_FOR_SCOPING */

// define if the compiler supports partial ordering
#define BZ_HAVE_PARTIAL_ORDERING

// define if the compiler supports partial specialization
#define BZ_HAVE_PARTIAL_SPECIALIZATION

// define if the compiler supports reinterpret_cast<>
#define BZ_HAVE_REINTERPRET_CAST

// define if the compiler supports Run-Time Type Identification
#define BZ_HAVE_RTTI

// define if the compiler has getrusage() function
#define BZ_HAVE_RUSAGE

// define if the compiler supports static_cast<>
#define BZ_HAVE_STATIC_CAST

// define if the compiler supports ISO C++ standard library
#define BZ_HAVE_STD

// Obsolete ?
// Define to 1 if you have the <stdint.h> header file.
/* #undef BZ_HAVE_STDINT_H */

// Define to 1 if you have the <stdlib.h> header file.
/* #undef BZ_HAVE_STDLIB_H */

// define if the compiler supports Standard Template Library
#define BZ_HAVE_STL

// Define to 1 if you have the <strings.h> header file.
/* #undef BZ_HAVE_STRINGS_H */

// Define to 1 if you have the <string.h> header file.
/* #undef BZ_HAVE_STRING_H */

// define if the compiler supports System V math library
/* #undef BZ_HAVE_SYSTEM_V_MATH */

// Define to 1 if you have the <sys/stat.h> header file.
/* #undef BZ_HAVE_SYS_STAT_H */

// Define to 1 if you have the <sys/types.h> header file.
/* #undef BZ_HAVE_SYS_TYPES_H */

// Define to 1 if you have the <tbb/atomic.h> header file.
/* #undef BZ_HAVE_TBB_ATOMIC_H */

// define if the compiler supports basic templates
#define BZ_HAVE_TEMPLATES

// define if the compiler supports templates as template arguments
#define BZ_HAVE_TEMPLATES_AS_TEMPLATE_ARGUMENTS

// define if the compiler supports use of the template keyword as a qualifier
#define BZ_HAVE_TEMPLATE_KEYWORD_QUALIFIER

// define if the compiler supports template-qualified base class specifiers
#define BZ_HAVE_TEMPLATE_QUALIFIED_BASE_CLASS

// define if the compiler supports template-qualified return types
#define BZ_HAVE_TEMPLATE_QUALIFIED_RETURN_TYPE

// define if the compiler supports function matching with argument types which
// are template scope-qualified
#define BZ_HAVE_TEMPLATE_SCOPED_ARGUMENT_MATCHING

// define if the compiler recognizes typename
#define BZ_HAVE_TYPENAME

// define if the compiler supports the vector type promotion mechanism
#define BZ_HAVE_TYPE_PROMOTION

// Define to 1 if you have the <unistd.h> header file.
/* #undef BZ_HAVE_UNISTD_H */

// define if the compiler supports numeric traits promotions
#define BZ_HAVE_USE_NUMTRAIT

// define if the compiler has valarray<T>
#define BZ_HAVE_VALARRAY

// define if the compiler has isnan function in namespace std
#define BZ_ISNAN_IN_NAMESPACE_STD

// define if the compiler has C math abs(integer types) in namespace std
#define BZ_MATH_ABSINT_IN_NAMESPACE_STD

// define if the compiler has C math functions in namespace std
#define BZ_MATH_FN_IN_NAMESPACE_STD

// Name of package
/* #undef BZ_PACKAGE */

// Define to the address where bug reports for this package should be sent.
/* #undef BZ_PACKAGE_BUGREPORT */

// Define to the full name of this package.
/* #undef BZ_PACKAGE_NAME */

// Define to the full name and version of this package.
#define BZ_PACKAGE_STRING " 1.0"

// Define to the one symbol short name of this package.
/* #undef BZ_PACKAGE_TARNAME */

// Define to the home page for this package.
/* #undef BZ_PACKAGE_URL */

// Define to the version of this package.
/* #undef BZ_PACKAGE_VERSION */

// Pad array lengths to SIMD width.
/* #undef BZ_PAD_ARRAYS */

// Set SIMD instruction width in bytes.
#define BZ_SIMD_WIDTH 1

// Define to 1 if you have the ANSI C header files.
/* #undef BZ_STDC_HEADERS */

// Enable Blitz thread-safety features
/* #undef BZ_THREADSAFE */

// Use TBB atomic types.
/* #undef BZ_THREADSAFE_USE_TBB */

// Specifies whether compiler alignment pragmas should be used.
/* #undef BZ_USE_ALIGNMENT_PRAGMAS */

// Version number of package
/* #undef BZ_VERSION */

// CXX
/* #undef BZ__compiler_name */

// CXXFLAGS
/* #undef BZ__compiler_options */

// date
/* #undef BZ__config_date */

// uname -a
/* #undef BZ__os_name */

// target
/* #undef BZ__platform */
