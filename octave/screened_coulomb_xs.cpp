#include <octave/oct.h>

#include "xs.h"

double exec_func(double x, double y)
{  
  return (std::isnan(x) || std::isnan(y)) ? std::numeric_limits<double>::quiet_NaN() : xs_quad<Screening::ZBL>::crossSection(x,y);
}

DEFUN_DLD (screened_coulomb_xs, args, nargout, "Center mass scattering angle for screened Coulomb potential")
{
  if (args.length () != 2)
    print_usage ();

  NDArray A = args(0).array_value ();
  NDArray B = args(1).array_value ();

  octave_idx_type na = A.numel ();
  octave_idx_type nb = B.numel ();

  if (na==nb) {
    NDArray C(A);
    for (octave_idx_type i=0; i<na; i++)
      C(i) = exec_func(A(i),B(i));
      return octave_value( C );
  }
  if (na==1) {
    NDArray C(B);
    for (octave_idx_type i=0; i<nb; i++)
      C(i) = exec_func(A(0),B(i));
      return octave_value( C );
  }
  if (nb==1) {
    NDArray C(A);
    for (octave_idx_type i=0; i<na; i++)
      C(i) = exec_func(A(i),B(0));
      return octave_value( C );
  }

  error("not compatible input dims");

  return octave_value ();
}