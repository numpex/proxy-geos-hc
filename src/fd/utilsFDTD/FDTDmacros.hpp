#ifndef FDTDMACROS_HPP
#define FDTDMACROS_HPP

#define POW2(x) ((x)*(x))
#define IDX3(i,j,k) (nz*ny*(i) + nz*(j) + (k))
#define IDX3_l(i,j,k) ((nz+2*lz)*(ny+2*ly)*((i)+lx) + (nz+2*lz)*((j)+ly) + ((k)+lz))
#define IDX3_eta1(i,j,k) ((nz+2)*(ny+2)*((i)+1) + (nz+2)*((j)+1) + ((k)+1))
#define IDX3_eta0(i,j,k) ((nz+2)*(ny+2)*(i) + (nz+2)*(j) + (k))

#ifdef USE_RAJA
// implementation of sqrt via binary search
// copied from https://stackoverflow.com/questions/8622256/in-c11-is-sqrt-defined-as-constexpr
constexpr size_t sqrt_helper(size_t n, size_t lo, size_t hi)
{
  return (lo == hi)
           ? lo // search complete
           : ((n / ((lo + hi + 1) / 2) < ((lo + hi + 1) / 2))
                ? sqrt_helper(n, lo, ((lo + hi + 1) / 2)-1) // search lower half
                : sqrt_helper(n, ((lo + hi + 1) / 2), hi)); // search upper half
}
// constexpr integer sqrt
constexpr size_t sqrt(size_t n)
{               
  return sqrt_helper(n, 0, n/2 + 1); 
}
// implementation of lesser_of_squarest_factor_pair via linear search
constexpr size_t lesser_of_squarest_factor_pair_helper(size_t n, size_t guess)
{
  return ((n / guess) * guess == n)
           ? guess // search complete, guess is a factor
           : lesser_of_squarest_factor_pair_helper(n, guess - 1); // continue searching
}
// constexpr return the lesser of the most square pair of factors of n
// ex. 12 has pairs of factors (1, 12) (2, 6) *(3, 4)* and returns 3
constexpr size_t lesser_of_squarest_factor_pair(size_t n)
{
  return (n == 0)
      ? 0 // return 0 in the 0 case
      : lesser_of_squarest_factor_pair_helper(n, sqrt(n));
}
// constexpr return the greater of the most square pair of factors of n
// ex. 12 has pairs of factors (1, 12) (2, 6) *(3, 4)* and returns 4
constexpr size_t greater_of_squarest_factor_pair(size_t n)
{
  return (n == 0)
      ? 0 // return 0 in the 0 case
      : n / lesser_of_squarest_factor_pair_helper(n, sqrt(n));
}

//#define block_size (512)
#define block_size (128)
#define z_block_sz (32)
#define y_block_sz (greater_of_squarest_factor_pair(block_size/z_block_sz))
#define x_block_sz (lesser_of_squarest_factor_pair(block_size/z_block_sz))

#ifdef ENABLE_CUDA
#define RAJANestedLoop(x3,y3,z3,x4,y4,z4)\
     RAJA::TypedRangeSegment<int> KRange(z3, z4);\
     RAJA::TypedRangeSegment<int> JRange(y3, y4);\
     RAJA::TypedRangeSegment<int> IRange(x3, x4);\
     using EXEC_POL =\
     RAJA::KernelPolicy<\
        RAJA::statement::CudaKernelFixedAsync<x_block_sz*y_block_sz*z_block_sz,\
          RAJA::statement::For<0, RAJA::cuda_global_size_z_direct<x_block_sz>, \
            RAJA::statement::For<1, RAJA::cuda_global_size_y_direct<y_block_sz>,\
              RAJA::statement::For<2, RAJA::cuda_global_size_x_direct<z_block_sz>,\
               RAJA::statement::Lambda<0,RAJA::Segs<0,1,2>>\
             >\
           >\
         >\
       >\
     >;\
     RAJA::kernel<EXEC_POL>( RAJA::make_tuple(IRange, JRange, KRange), [=] __device__ (int i, int j, int k) 
#endif


#ifdef ENABLE_HIP
#define RAJANestedLoop(x3,y3,z3,x4,y4,z4)\
     RAJA::TypedRangeSegment<int> KRange(z3, z4);\
     RAJA::TypedRangeSegment<int> JRange(y3, y4);\
     RAJA::TypedRangeSegment<int> IRange(x3, x4);\
     using EXEC_POL =\
     RAJA::KernelPolicy<\
        RAJA::statement::HipKernelFixedAsync<x_block_sz*y_block_sz*z_block_sz,\
          RAJA::statement::For<0, RAJA::hip_global_size_z_direct<x_block_sz>, \
            RAJA::statement::For<1, RAJA::hip_global_size_y_direct<y_block_sz>,\
              RAJA::statement::For<2, RAJA::hip_global_size_x_direct<z_block_sz>,\
               RAJA::statement::Lambda<0,RAJA::Segs<0,1,2>>\
             >\
           >\
         >\
       >\
     >;\
     RAJA::kernel<EXEC_POL>( RAJA::make_tuple(IRange, JRange, KRange), [=] __device__ (int i, int j, int k) 
#endif
#endif

#if defined (USE_RAJA)
  #define LOOP3DHEAD RAJANestedLoop(x3,y3,z3,x4,y4,z4) {
#elif defined (USE_KOKKOS)
  #define LOOP3DHEAD Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({z3,x3,y3},{z4,x4,y4}),KOKKOS_LAMBDA(int k,int i,int j) {
#else
   #define LOOP3DHEAD\
      for (int i = x3; i < x4; ++i)\
      {\
          for (int j = y3; j < y4; ++j)\
          {\
              for (int k = z3; k < z4; ++k)\
              {
#endif

#if defined(USE_RAJA) || defined(USE_KOKKOS)
  #define LOOP3DEND   });
  #define VECTORVIEW vectorRealView const &
#else
  #define LOOP3DEND }}}
  #define VECTORVIEW vectorRealView &
#endif

#endif //FDTDMACROS_HPP

