#pragma once

#include <Eigen/Core>

#include <csCore/BasicTypes.hpp>
#include <csLie/Types.hpp>

namespace csLie {
  
  /**
   * @brief transform a S2 unit vector into a SO3 rotation matrix with x aligned to the unit vector
   * @tparam S scalar type
   * @param v the unit vector
   * @return a rotation matrix
   */
  template <class S> SO3<S> s2toZAlignedSO3(S2<S> const & v);

  // given a s2 algebra vector delta (2D vector)
  // return the S2 exp, i.e., a unitary 3D vector
  template <class S> S2<S> s2Exp(S2Alg<S> const & delta);

  // given 
  // - a unitary 3D vector v (in S2)
  // - a s2 algebra vector delta (2D vector)
  // return the oplus(v, delta)
  // equivalent to s2toZAlignedSO3(v) * s2Exp(delta)
  // for readability rewritten as
  // R(V) * exp(delta) 
  // note: repeating oplus in loops might introduce numerical errors, renormalize the vector when needed.
  template <class S> S2<S> s2OPlus(S2<S> const& v, S2Alg<S> const& delta);

  // the typical operation `oplus` given v : S2 and delta : S2Alg is defined as
  // v_p = s2toZAlignedSO3(v) * s2Exp(delta)
  // for readability rewritten as
  // v_p = R(V) * exp(delta) 
  // in an optimization framework, is often need to evaluate the jacobian wrt the delta evaluated at zero
  //  d v_p  |
  // ------- |         
  // d delta |delta=0
  // this is equal to the last two column of R(v), which are returned by this function
  template <class S> S2OPlusJac<S> s2_d_OPlus_d_delta_0(S2<S> const& v);

}