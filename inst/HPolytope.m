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
## @seealso{HPolytope.sample, HPolytope.volume, HPolytope.is_in}
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
    ## @var{starting_point} is a d × 1 interior point of the polytope.  If
    ## omitted or empty, the Chebyshev center (the center of the largest
    ## inscribed ball) is computed automatically by an interior-point method.
    ##
    ## @var{seed} is an optional integer seed for the random number generator.
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
      if (! isnumeric (walk_len) || ! isscalar (walk_len) ...
          || walk_len < 1 || walk_len != round (walk_len))
        error ("HPolytope.sample: WALK_LEN must be a positive integer.");
      endif
      if (! isnumeric (num_points) || ! isscalar (num_points) ...
          || num_points < 1 || num_points != round (num_points))
        error ("HPolytope.sample: NUM_POINTS must be a positive integer.");
      endif

      if (nargin < 4 || isempty (starting_point))
        starting_point = [];
      endif
      if (nargin < 5 || isempty (seed))
        seed = -1;
      endif

      pts = __volesti_HPolytope__ ("sample", this.A, this.b, ...
                                   walk_len, num_points, ...
                                   starting_point, seed);

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
    ## @var{seed} is an optional integer seed for the random number generator.
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
      if (nargin < 4 || isempty (seed))
        seed = -1;
      endif

      if (! isnumeric (walk_len) || ! isscalar (walk_len) ...
          || walk_len < 1 || walk_len != round (walk_len))
        error ("HPolytope.volume: WALK_LEN must be a positive integer.");
      endif
      if (! isnumeric (error_tol) || ! isscalar (error_tol) || error_tol <= 0)
        error ("HPolytope.volume: ERROR must be a positive scalar.");
      endif

      vol = __volesti_HPolytope__ ("volume", this.A, this.b, ...
                                   walk_len, error_tol, seed);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {HPolytope} {@var{res} =} is_in (@var{P}, @var{x})
    ## @deftypefnx {HPolytope} {@var{res} =} is_in (@var{P}, @var{x}, @var{tol})
    ##
    ## Test whether point(s) lie inside the H-polytope.
    ##
    ## @code{@var{res} = is_in (@var{P}, @var{x})} checks if the point @var{x}
    ## (a @math{d \times 1} column vector or a @math{d \times N} matrix of @var{N}
    ## points) satisfies the polytope linear constraints @math{@var{A} * x <= @var{b}}.
    ##
    ## @var{tol} is an optional non-negative numeric tolerance (default: 0.0).
    ## A point @var{x} is considered inside if @math{@var{A} * x <= @var{b} + @var{tol}}.
    ##
    ## Returns a logical scalar if @var{x} is a single point, or a @math{1 \times N}
    ## logical row vector if @var{x} contains multiple points.
    ##
    ## @seealso{HPolytope.sample, HPolytope.volume}
    ## @end deftypefn
    function res = is_in (this, x, tol)

      ## Input validation
      if (nargin < 2)
        error ("HPolytope.is_in: X is required.");
      endif
      if (! isnumeric (x) || ! ismatrix (x) || isempty (x))
        error ("HPolytope.is_in: X must be a non-empty numeric vector or matrix.");
      endif
      if (! isreal (x))
        error ("HPolytope.is_in: X must be real.");
      endif
      if (any (isnan (x(:))) || any (isinf (x(:))))
        error ("HPolytope.is_in: X must not contain NaN or Inf.");
      endif

      d = columns (this.A);

      if (rows (x) != d)
        if (isvector (x))
          if (d == 1)
            x = x.';
          elseif (length (x) == d)
            x = x(:);
          else
            error ("HPolytope.is_in: dimension of X (%d) does not match polytope dimension (%d).", ...
                   rows (x), d);
          endif
        else
          error ("HPolytope.is_in: dimension of X (%d) does not match polytope dimension (%d).", ...
                 rows (x), d);
        endif
      endif

      if (nargin < 3 || isempty (tol))
        tol = 0.0;
      endif
      if (! isnumeric (tol) || ! isscalar (tol) || tol < 0)
        error ("HPolytope.is_in: TOL must be a non-negative numeric scalar.");
      endif
      if (! isreal (tol))
        error ("HPolytope.is_in: TOL must be a real scalar.");
      endif
      if (isnan (tol) || isinf (tol))
        error ("HPolytope.is_in: TOL must not be NaN or Inf.");
      endif

      res = __volesti_HPolytope__ ("is_in", this.A, this.b, x, tol);

    endfunction

  endmethods

endclassdef
