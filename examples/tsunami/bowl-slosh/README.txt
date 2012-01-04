begin_html  [use: doc/doc.css] [use:jsMath]
<!--   For a more readable version of this file, execute
                  unix>  make htmls
       in this directory and then point your browser to README.html 
     --------------------------------------------------------------  -->

<h2>
GeoClaw Sample Code
</h2>

Waves in a parabolic bowl with a flat surface sloshing around.
An exact analytic solution is known in which the surface stays flat.

In this code, $x$ and $y$ are in meters (iccordsys=1 in [code:setrun.py]).

Topography: $B(x,y) = h_0((x^2 + y^2)/a^2 -1)$,

Depth: $h(x,y,t) = \max\left(0,~~ (\sigma h_0/a^2)(2x\cos(\omega t) + 2y\sin(\omega t) -
\sigma) - B(x,y)\right)$

Velocities:  $u(x,y,t) = -\sigma \omega \sin(\omega t),\qquad
v(x,y,t) = \sigma \omega \cos(\omega t).$

where $\omega = \sqrt{2gh_0} / a$.

The period of oscillation is  $T = 2\pi / \omega$.

The following parameters are currently hardwired several places:

$a = 1, ~~\sigma = 0.5, ~~h = 0.1,~~g = 9.81$ 

This should be cleaned up: better to put them in a setprob.data file that
is read in where needed.

<h3>References</h3>

<ul>
<li> W. C. Thacker, Some exact solutions to the nonlinear shallow water wave equations,
J. Fluid Mech. 107 (1981), 499-508.
<li> J.M. Gallardo, C. Pares, and M. Castro, On a well-balanced high-order
finite volume scheme for shallow water equations with topography and dry
areas, J. Comput. Phys. 227(2007) 574-601.
<li> Y. Xing, X. Zhang and C.-W. Shu, Positivity preserving high order well
balanced discontinuous Galerkin methods for the shallow water equations ,
submitted to Advances in Water Resources. 
[http://www.dam.brown.edu/scicomp/scg-media/report_files/BrownSC-2010-10.pdf preprint]
</ul>
<p>
This test problem has been used in several other papers too.

To run the code, see [link: #instructions Instructions]

<h4>
Plots of results
</h4>
After running this code and creating plots via "make .plots", you should be
able to view the plots in [link: _plots/_PlotIndex.html].


<h4>
Fortran files
</h4>


<dl>
<dt>[code: Makefile]
<dd> Determines which version of fortran files
are used when compiling the code with make and specifies where output and
plots should be directed.  Type "make .help" at the Unix prompt for options.



</dl>

<h4>
Python files
</h4>
<dl>

<dt>[code: maketopo.py]
<dd> Used to create topo file and qinit data file.

<dt>[code: setrun.py]
<dd> This file contains a function that 
specifies what run-time parameters will be used.

<dt>[code: setplot.py]
<dd> This file contains a function that 
specifies what plots will be done and
sets various plotting parameters. 

</dl>


<h4>
Data files
</h4>

The .data files are automatically generated using the information in 
[code: setrun.py].


[name:  instructions]
<h4>
Instructions
</h4>


To make topo and qinit data files:
{{{
  $ make topo
}}}

To make all data files, edit setrun.py and then
{{{
  $ make .data
}}}

To run code:

You may need to type
{{{
  $ make new 
}}}
to make sure the modules are accessible.

Then run the code with
{{{
  $ make .output
}}}

To plot results, either generate html pages via:
{{{
  $ make .plots
}}}
or view interactively using ipython and Iplotclaw.

All of this can be done with:

{{{
  $ source make_all.py
}}}

end_html

