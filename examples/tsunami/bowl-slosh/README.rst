
.. _geoclaw_examples_tsunami_bowl-slosh:

Sloshing water in a parabolic bowl 
==================================

Waves in a parabolic bowl with a flat surface sloshing around.
An exact analytic solution is known in which the surface stays flat.

To create the topo file before running the code::

    make topo

In this code, :math:`x` and :math:`y` are in meters (coordinate_system=1 
in `setrun.py`).

Topography: :math:`B(x,y) = h_0((x^2 + y^2)/a^2 -1)`,

Depth: :math:`h(x,y,t) = \max\left(0,~~ (\sigma h_0/a^2)(2x\cos(\omega t) + 2y\sin(\omega t) -
\sigma) - B(x,y)\right)`

Velocities:  :math:`u(x,y,t) = -\sigma \omega \sin(\omega t),\qquad
v(x,y,t) = \sigma \omega \cos(\omega t).`

where :math:`\omega = \sqrt{2gh_0} / a`.

The period of oscillation is  :math:`T = 2\pi / \omega`.

The following parameters are currently hardwired several places:

:math:`a = 1, ~~\sigma = 0.5, ~~h = 0.1,~~g = 9.81` 

This should be cleaned up: better to put them in a setprob.data file that
is read in where needed.

References
----------

* W. C. Thacker, Some exact solutions to the nonlinear shallow water wave equations,
  J. Fluid Mech. 107 (1981), 499-508.

* J.M. Gallardo, C. Pares, and M. Castro, On a well-balanced high-order
  finite volume scheme for shallow water equations with topography and dry
  areas, J. Comput. Phys. 227(2007) 574-601.

* Y. Xing, X. Zhang and C.-W. Shu, Positivity preserving high order well
  balanced discontinuous Galerkin methods for the shallow water equations ,
  Advances in Water Resources  33 (2010), pp. 1476-1493. 

This test problem has been used in several other papers too.

