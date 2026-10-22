/*
  Copyright (C) 2026 Vissarion Fisikopoulos

  This file is part of the ovolesti package for GNU Octave.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 3 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, see <http://www.gnu.org/licenses/>.
*/

#include <octave/oct.h>
#include <string>
#include <list>
#include <cmath>
#include <new>
#include <stdexcept>

using namespace std;

// Let Ctrl-C interrupt volesti's random walks: they call this hook once per
// step, and octave_quit throws octave::interrupt_exception when an interrupt
// is pending.  It must be defined before any volesti header is included.
#define VOLESTI_CHECK_INTERRUPT() octave_quit ()

#include <Eigen/Eigen>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "sampling/random_point_generators.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "random_walks/uniform_cdhr_walk.hpp"
#include "generators/boost_random_number_generator.hpp"


// ---------------------------------------------------------------------------
//  uniform_sampling (extracted from volesti sampling/sampling.hpp to avoid
//  pulling in CRHMC and ODE-solver dependency chains we don't need)
// ---------------------------------------------------------------------------

/// Generate N uniform samples from a polytope using a random walk.
/// PushBackWalkPolicy is provided by volume/sampling_policies.hpp.
template <typename WalkTypePolicy,
          typename PointList,
          typename Polytope,
          typename RandomNumberGenerator,
          typename Point>
void
uniform_sampling (PointList &randPoints,
                  Polytope &P,
                  RandomNumberGenerator &rng,
                  const unsigned int &walk_len,
                  const unsigned int &rnum,
                  const Point &starting_point,
                  unsigned int const& nburns)
{
  typedef typename WalkTypePolicy::template Walk
          <Polytope, RandomNumberGenerator> walk;

  PushBackWalkPolicy push_back_policy;
  Point p = starting_point;

  typedef RandomPointGenerator<walk> RandomPointGenerator;
  if (nburns > 0)
    {
      RandomPointGenerator::apply (P, p, nburns, walk_len, randPoints,
                                   push_back_policy, rng);
      randPoints.clear ();
    }
  RandomPointGenerator::apply (P, p, rnum, walk_len, randPoints,
                               push_back_policy, rng);
}


typedef Cartesian<double>                 Kernel;
typedef Kernel::Point                     Point;
typedef HPolytope<Point>                  HPoly;
typedef BoostRandomNumberGenerator
        <boost::mt19937, double>          RNG;


// ---------------------------------------------------------------------------
//  Helpers
// ---------------------------------------------------------------------------

// Convert an octave (m × d) double matrix to an Eigen dense matrix.
static Eigen::MatrixXd
oct_to_eigen_matrix (const Matrix &m)
{
  octave_idx_type rows = m.rows ();
  octave_idx_type cols = m.columns ();
  Eigen::MatrixXd E (rows, cols);
  for (octave_idx_type j = 0; j < cols; ++j)
    for (octave_idx_type i = 0; i < rows; ++i)
      E (i, j) = m(i, j);
  return E;
}

// Convert an octave column vector to an Eigen vector.
static Eigen::VectorXd
oct_to_eigen_vector (const Matrix &v)
{
  octave_idx_type n = v.numel ();
  Eigen::VectorXd E (n);
  for (octave_idx_type i = 0; i < n; ++i)
    E (i) = v(i);
  return E;
}

// Convert an Eigen vector to a Point (volesti).
static Point
eigen_to_point (const Eigen::VectorXd &v)
{
  return Point (v);
}

// Convert a list of volesti Points to a (d × N) octave matrix.
static Matrix
point_list_to_matrix (const list<Point> &randPoints)
{
  if (randPoints.empty ())
    return Matrix (0, 0);

  unsigned int d = randPoints.front ().dimension ();
  Matrix M (d, randPoints.size ());
  octave_idx_type col = 0;
  for (const auto &pt : randPoints)
    {
      for (unsigned int i = 0; i < d; ++i)
        M(i, col) = pt[i];
      ++col;
    }
  return M;
}


// volesti reports a failed inner-ball computation with a negative radius.
// HPolytope.m rejects empty, flat and unbounded polytopes first, so this only
// fires on numerically hard inputs.
static void
compute_inner_ball (HPoly &P, const char *scope)
{
  if (P.ComputeInnerBall ().second < 0)
    error ("%s: unable to compute an interior point of the polytope.", scope);
}


// ---------------------------------------------------------------------------
//  Method dispatchers
// ---------------------------------------------------------------------------

static octave_value_list
do_sample (const octave_value_list &args)
{
  // --- unpack arguments ---
  // args(0): method, args(1): A, args(2): b, args(3): walk_len,
  // args(4): num_points, args(5): starting_point, args(6): seed
  Matrix A_oct = args(1).matrix_value ();
  Matrix b_oct = args(2).matrix_value ();
  unsigned int walk_len = args(3).uint_value ();
  unsigned int num_points = args(4).uint_value ();
  double seed = args(6).double_value ();

  // Check for empty starting point (args(5))
  bool have_start = false;
  Eigen::VectorXd start_vec;
  if (args(5).is_matrix_type ())
    {
      Matrix start_oct = args(5).matrix_value ();
      if (start_oct.numel () > 0)
        {
          start_vec = oct_to_eigen_vector (start_oct);
          have_start = true;
        }
    }

  // --- construct HPolytope ---
  Eigen::MatrixXd Ae = oct_to_eigen_matrix (A_oct);
  Eigen::VectorXd be = oct_to_eigen_vector (b_oct);
  unsigned int d = static_cast<unsigned int> (Ae.cols ());
  HPoly P (d, Ae, be);

  // --- interior point ---
  // HPolytope.m has checked the size and finiteness of STARTING_POINT, so
  // volesti can only reject it for not being strictly inside.
  Point start_pt (d);
  if (have_start)
    {
      start_pt = eigen_to_point (start_vec);
      try
        {
          P.set_interior_point (start_pt);
        }
      catch (const std::invalid_argument &)
        {
          error ("HPolytope.sample: STARTING_POINT must lie strictly "
                 "inside the polytope.");
        }
    }
  else
    {
      compute_inner_ball (P, "HPolytope.sample");
      start_pt = P.InnerBall ().first;
    }

  // --- random number generator ---
  RNG rng (d);
  if (seed >= 0)
    rng.set_seed (static_cast<unsigned> (seed));

  // --- sample ---
  list<Point> randPoints;
  uniform_sampling<CDHRWalk> (randPoints, P, rng,
                              walk_len, num_points, start_pt, 0);

  // --- return ---
  octave_value_list ret;
  ret(0) = octave_value (point_list_to_matrix (randPoints));
  return ret;
}

static octave_value_list
do_volume (const octave_value_list &args)
{
  // --- unpack arguments ---
  Matrix A_oct = args(1).matrix_value ();
  Matrix b_oct = args(2).matrix_value ();
  unsigned int walk_len = args(3).uint_value ();
  double error_tol = args(4).double_value ();
  double seed = args(5).double_value ();

  // --- construct HPolytope ---
  Eigen::MatrixXd Ae = oct_to_eigen_matrix (A_oct);
  Eigen::VectorXd be = oct_to_eigen_vector (b_oct);
  unsigned int d = static_cast<unsigned int> (Ae.cols ());
  HPoly P (d, Ae, be);

  // --- inner ball ---
  compute_inner_ball (P, "HPolytope.volume");

  // --- random number generator ---
  RNG rng (d);
  if (seed >= 0)
    rng.set_seed (static_cast<unsigned> (seed));

  // --- compute volume ---
  double vol = volume_sequence_of_balls<CDHRWalk> (P, rng, error_tol,
                                                   walk_len, 1);
  if (! (vol > 0) || ! std::isfinite (vol))
    error ("HPolytope.volume: volume estimation failed.");

  // --- return ---
  octave_value_list ret;
  ret(0) = octave_value (vol);
  return ret;
}


// ---------------------------------------------------------------------------
//  Main entry point
// ---------------------------------------------------------------------------

DEFUN_DLD (__volesti_HPolytope__, args, nargout,
           "-*- texinfo -*-\n\
@deftypefn {ovolesti} {} __volesti_HPolytope__ (@var{method}, @var{A}, \
@var{b}, @dots{})\n\
\n\
Internal oct-file binding for the @qcode{HPolytope} class.\n\
Dispatches on @var{method} to call the appropriate volesti C++ function.\n\
\n\
Users should call the class methods @code{HPolytope.sample} and\n\
@code{HPolytope.volume} instead of calling this function directly.\n\
\n\
@seealso{HPolytope}\n\
@end deftypefn")
{
  octave_value_list retval;

  if (args.length () < 2)
    error ("__volesti_HPolytope__: at least METHOD and A, b are required.");

  string method = args(0).string_value ();

  if (method == "sample")
    {
      if (args.length () != 7)
        error ("__volesti_HPolytope__: 'sample' requires A, b, WALK_LEN, "
               "NUM_POINTS, STARTING_POINT, SEED.");
    }
  else if (method == "volume")
    {
      if (args.length () != 6)
        error ("__volesti_HPolytope__: 'volume' requires A, b, WALK_LEN, "
               "ERROR, SEED.");
    }
  else
    error ("__volesti_HPolytope__: unknown method '%s'.", method.c_str ());

  // Translate volesti's C++ exceptions into Octave errors.  Octave's own
  // exceptions (error, Ctrl-C, exit) and bad_alloc pass through: the first
  // derives from std::runtime_error and the others from std::exception, so
  // they must be rethrown before the generic handler.
  try
    {
      if (method == "sample")
        retval = do_sample (args);
      else
        retval = do_volume (args);
    }
  catch (const octave::execution_exception &)
    {
      throw;
    }
  catch (const octave::interrupt_exception &)
    {
      throw;
    }
  catch (const octave::exit_exception &)
    {
      throw;
    }
  catch (const std::bad_alloc &)
    {
      throw;
    }
  catch (const std::exception &e)
    {
      error ("HPolytope.%s: %s.", method.c_str (), e.what ());
    }

  return retval;
}
