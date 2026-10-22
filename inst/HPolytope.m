## Copyright (C) 2026 Vissarion Fisikopoulos
##
## This file is part of the ovolesti package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftp {ovolesti} HPolytope
##
## An H-representation convex polytope defined by a set of linear inequalities
## @math{@var{A} * x <= @var{b}}.
##
## An @qcode{HPolytope} object holds the constraint matrix @var{A} (each row is
## a halfspace normal) and the right-hand side vector @var{b}.  Methods
## @code{sample} and @code{volume} delegate to the volesti C++ library via
## compiled oct-file bindings.
##
## @code{@var{P} = HPolytope (@var{A}, @var{b})} constructs an H-polytope from
## the constraint matrix @var{A} (m × d) and vector @var{b} (m × 1).
##
## Properties:
##
## @deftp {HPolytope} {property} A
## Constraint matrix (m × d).  Each row is the normal of a halfspace.
## @end deftp
##
## @deftp {HPolytope} {property} b
## Right-hand side vector (m × 1).
## @end deftp
##
## @seealso{HPolytope.sample, HPolytope.volume}
## @end deftp
classdef HPolytope

  properties (SetAccess = private)
    A = []
    b = []
  endproperties

  methods

    ## -*- texinfo -*-
    ## @deftypefn {HPolytope} {@var{P} =} HPolytope (@var{A}, @var{b})
    ##
    ## Construct an H-polytope from a set of linear inequalities.
    ##
    ## @code{@var{P} = HPolytope (@var{A}, @var{b})} returns an @qcode{HPolytope}
    ## object representing the intersection of halfspaces
    ## @math{@var{A} * x <= @var{b}}.
    ##
    ## @var{A} must be a numeric matrix (m constraints × d dimension).
    ## @var{b} must be a numeric column vector with m elements.
    ##
    ## @end deftypefn
    function this = HPolytope (A, b)

      ## Input validation
      if (nargin != 2)
        error ("HPolytope: exactly two arguments (A, b) are required.");
      endif
      if (! isnumeric (A) || ! ismatrix (A))
        error ("HPolytope: A must be a numeric matrix.");
      endif
      if (isempty (A) || rows (A) == 0 || columns (A) == 0)
        error ("HPolytope: A must be non-empty.");
      endif
      if (! isnumeric (b) || ! isvector (b) || isempty (b))
        error ("HPolytope: b must be a non-empty numeric vector.");
      endif
      if (! isreal (A))
        error ("HPolytope: A must be real.");
      endif
      if (! isreal (b))
        error ("HPolytope: b must be real.");
      endif

      ## Coerce b to column vector
      b = b(:);
      m = size (A, 1);

      if (any (isnan (A(:))) || any (isinf (A(:))))
        error ("HPolytope: A must not contain NaN or Inf.");
      endif
      if (any (isnan (b)) || any (isinf (b)))
        error ("HPolytope: b must not contain NaN or Inf.");
      endif
      if (length (b) != m)
        error ("HPolytope: b must have the same number of rows as A.");
      endif

      this.A = A;
      this.b = b;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {HPolytope} {@var{pts} =} sample (@var{P}, @var{walk_len}, @var{num_points})
    ## @deftypefnx {HPolytope} {@var{pts} =} sample (@var{P}, @var{walk_len}, @
    ##                                    @var{num_points}, @var{starting_point})
    ## @deftypefnx {HPolytope} {@var{pts} =} sample (@var{P}, @var{walk_len}, @
    ##                                    @var{num_points}, @var{starting_point}, @
    ##                                    @var{seed})
    ##
    ## Draw uniform samples from the polytope.
    ##
    ## @code{@var{pts} = sample (@var{P}, @var{walk_len}, @var{num_points})}
    ## generates @var{num_points} uniformly distributed points from
    ## @var{P}.  The random walk takes @var{walk_len} steps between each sample.
    ##
    ## @var{walk_len} and @var{num_points} must be integers between 1 and
    ## @math{2^32 - 1}.
    ##
    ## @var{starting_point} is a real vector with d elements that lies strictly
    ## inside the polytope.  If omitted or empty, the Chebyshev center (the
    ## center of the largest inscribed ball) is computed automatically by an
    ## interior-point method.
    ##
    ## @var{seed} is an integer between 0 and @math{2^32 - 1} that seeds the
    ## random number generator, so that repeated calls return the same points.
    ## If omitted or empty, the generator is seeded from the clock.
    ##
    ## @var{P} must be bounded and have a nonempty interior; otherwise an error
    ## is raised.
    ##
    ## The return value @var{pts} is a d × @var{num_points} matrix where each
    ## column is a sample point.
    ##
    ## @seealso{HPolytope.volume}
    ## @end deftypefn
    function pts = sample (this, walk_len, num_points, starting_point, seed)

      ## Input validation
      if (nargin < 3)
        error ("HPolytope.sample: WALK_LEN and NUM_POINTS are required.");
      endif
      if (! isCount (walk_len))
        error (strcat ("HPolytope.sample: WALK_LEN must be an integer", ...
                       " between 1 and 2^32 - 1."));
      endif
      if (! isCount (num_points))
        error (strcat ("HPolytope.sample: NUM_POINTS must be an integer", ...
                       " between 1 and 2^32 - 1."));
      endif
      if (nargin < 5)
        seed = [];
      endif
      seed = checkSeed (seed, "HPolytope.sample");
      checkPolytope (this.A, this.b, "HPolytope.sample");

      ## Validate the starting point against the polytope
      if (nargin < 4 || isempty (starting_point))
        starting_point = [];
      else
        d = columns (this.A);
        if (! isnumeric (starting_point) || ! isreal (starting_point) ...
            || ! isvector (starting_point) || numel (starting_point) != d)
          error (strcat ("HPolytope.sample: STARTING_POINT must be a real", ...
                         " vector with %d elements."), d);
        endif
        starting_point = double (starting_point(:));
        if (! all (isfinite (starting_point)))
          error (strcat ("HPolytope.sample: STARTING_POINT must not", ...
                         " contain NaN or Inf."));
        endif
      endif

      pts = __volesti_HPolytope__ ("sample", this.A, this.b, walk_len, ...
                                   num_points, starting_point, seed);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {HPolytope} {@var{vol} =} volume (@var{P})
    ## @deftypefnx {HPolytope} {@var{vol} =} volume (@var{P}, @var{walk_len})
    ## @deftypefnx {HPolytope} {@var{vol} =} volume (@var{P}, @var{walk_len}, @
    ##                                    @var{error})
    ## @deftypefnx {HPolytope} {@var{vol} =} volume (@var{P}, @var{walk_len}, @
    ##                                    @var{error}, @var{seed})
    ##
    ## Estimate the volume of the polytope.
    ##
    ## @code{@var{vol} = volume (@var{P})} returns an estimate of the volume of
    ## @var{P} using the sequence-of-balls algorithm with default parameters.
    ##
    ## @var{walk_len} is the number of steps per random-walk iteration
    ## (default: 1).
    ##
    ## @var{error} is the target relative error (default: 1.0).
    ##
    ## @var{seed} is an integer between 0 and @math{2^32 - 1} that seeds the
    ## random number generator, so that repeated calls return the same
    ## estimate.  If omitted or empty, the generator is seeded from the clock.
    ##
    ## @var{P} must be bounded and have a nonempty interior; otherwise an error
    ## is raised.
    ##
    ## @seealso{HPolytope.sample}
    ## @end deftypefn
    function vol = volume (this, walk_len, error_tol, seed)

      if (nargin < 2 || isempty (walk_len))
        walk_len = 1;
      endif
      if (nargin < 3 || isempty (error_tol))
        error_tol = 1.0;
      endif
      if (nargin < 4)
        seed = [];
      endif

      if (! isCount (walk_len))
        error (strcat ("HPolytope.volume: WALK_LEN must be an integer", ...
                       " between 1 and 2^32 - 1."));
      endif
      if (! isnumeric (error_tol) || ! isreal (error_tol) ...
          || ! isscalar (error_tol) || ! isfinite (error_tol) || error_tol <= 0)
        error ("HPolytope.volume: ERROR must be a positive finite scalar.");
      endif
      seed = checkSeed (seed, "HPolytope.volume");
      checkPolytope (this.A, this.b, "HPolytope.volume");

      vol = __volesti_HPolytope__ ("volume", this.A, this.b, walk_len, ...
                                   error_tol, seed);

    endfunction

  endmethods

endclassdef

## True if X is a real integer scalar that fits in a C++ unsigned int and is
## at least 1 (used for WALK_LEN and NUM_POINTS).
function tf = isCount (x)
  tf = isnumeric (x) && isreal (x) && isscalar (x) ...
       && x >= 1 && x <= double (intmax ('uint32')) && x == fix (x);
endfunction

## Validate SEED and map it to the oct-file convention (-1 means unseeded).
function seed = checkSeed (seed, scope)
  if (isempty (seed))
    seed = -1;
  elseif (! isnumeric (seed) || ! isreal (seed) || ! isscalar (seed) ...
          || ! (seed >= 0 && seed <= double (intmax ('uint32'))) ...
          || seed != fix (seed))
    error ("%s: SEED must be an integer between 0 and 2^32 - 1.", scope);
  else
    seed = double (seed);
  endif
endfunction

## Raise an error unless the polytope {x : A*x <= b} is bounded and has a
## nonempty interior.  volesti assumes both, and aborts or returns garbage
## otherwise.  Two small LPs decide this:
##  - the Chebyshev ball, max r s.t. A*x + norm (A(i,:))*r <= b, r >= 0, is
##    infeasible for an empty polytope and has r = 0 for a flat one;
##  - by Stiemke's lemma, A*y <= 0 has no solution with A*y != 0 if and only
##    if A'*l = 0 for some l > 0, so with rank (A) = d the polytope is bounded.
function checkPolytope (A, b, scope)
  A = full (double (A));
  b = full (double (b));
  [m, d] = size (A);
  param = struct ('msglev', 0, 'presol', 0);
  GLP_NOFEAS = 4;
  GLP_OPT = 5;
  GLP_UNBND = 6;
  ## Chebyshev ball
  [~, r, errnum, extra] = glpk ([zeros(d, 1); 1], [A, sqrt(sumsq (A, 2))], ...
                                b, [-Inf(d, 1); 0], [], repmat ('U', 1, m), ...
                                repmat ('C', 1, d + 1), -1, param);
  if (errnum == 0 && extra.status == GLP_NOFEAS)
    error ("%s: the polytope is empty.", scope);
  elseif (errnum == 0 && extra.status == GLP_UNBND)
    error ("%s: the polytope is unbounded.", scope);
  elseif (errnum != 0)
    error ("%s: unable to check the polytope (glpk error %d).", scope, errnum);
  elseif (extra.status != GLP_OPT)
    error ("%s: unable to check the polytope (glpk status %d).", scope, ...
           extra.status);
  endif
  if (r <= 1e-8)
    error ("%s: the polytope has an empty interior.", scope);
  endif
  ## Boundedness
  [~, ~, errnum, extra] = glpk (zeros (m, 1), A', zeros (d, 1), ones (m, 1), ...
                                [], repmat ('S', 1, d), repmat ('C', 1, m), ...
                                1, param);
  if (errnum != 0)
    error ("%s: unable to check the polytope (glpk error %d).", scope, errnum);
  elseif (extra.status != GLP_OPT || rank (A) < d)
    error ("%s: the polytope is unbounded.", scope);
  endif
endfunction
