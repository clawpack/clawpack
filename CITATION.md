# Citing Clawpack

If you use Clawpack in publications, please cite the following:

    @misc{clawpack,
          title={Clawpack software}, 
          author={Clawpack Development Team}, 
          url={http://www.clawpack.org}, 
          note={Version 5.1},
          year={2014}}

If you are using a version of the software other than the one above, please make
sure to cite that instead.  Please also cite at least one of the following 
regarding the algorithms used in Clawpack:


#### Basic algorithms in 1D and 2D

    R. J. LeVeque, 1997. Wave propagation algorithms for multi-dimensional 
    hyperbolic systems. J. Comput. Phys. 131, 327–353.


    @article{rjl:wpalg,
        Author = {R. J. LeVeque},
        Title = {Wave propagation algorithms for multi-dimensional hyperbolic 
                 systems},
        Journal = {J. Comput. Phys.},
        Pages = {327--353},
        Volume = {131},
        Year = {1997}
    }

    R. J. LeVeque. Finite Volume Methods for Hyperbolic Problems. Cambridge 
    University Press, Cambridge, UK, 2002.

    @book{LeVeque-FVMHP,
          Author = {R. J. LeVeque},
          Title = {Finite Volume Methods for Hyperbolic Problems},
          Publisher = {Cambridge University Press},
          Year = {2002},
          Url = {http://www.clawpack.org/book.html}
    }

#### 3D algorithms
    
    J. O. Langseth and R. J. LeVeque. 2000. A wave-propagation method for 
    three-dimensional hyperbolic conservation laws. J. Comput. Phys. 165, 
    126–166.
    
    @article{LangsethLeVeque00,
             Author = {J. O. Langseth and R. J. LeVeque},
             Title = {A wave-propagation method for three-dimensional hyperbolic
                      conservation laws},
             Journal = {J. Comput. Phys.},
             Pages = {126--166},
             Volume = {165},
             Year = {2000}
    }

#### Adaptive Mesh Refinement (AMR)

    M. J. Berger and R. J. LeVeque. 1998. Adaptive Mesh Refinement using Wave-Propagation Algorithms for Hyperbolic Systems. SIAM J. Numer. Anal. 35, 2298–2316.

    @article{BergerLeVeque98,
             Author = {M. J. Berger and R. J. LeVeque},
             Journal = {SIAM J. Numer. Anal.},
             Pages = {2298--2316},
             Title = {Adaptive Mesh Refinement using Wave-Propagation Algorithms 
                      for Hyperbolic Systems},
             Volume = {35},
             Year = {1998}
    }
 
#### F-wave Algorithms

    D. S. Bale, R. J. LeVeque, S. Mitran, and J. A. Rossmanith. A wave-propagation 
    method for conservation laws with spatially varying flux functions, SIAM J. 
    Sci. Comput 24 (2002), 955-978.

    @article{BaleLevMitRoss02,
        Author = {D. Bale and R. J. LeVeque and S. Mitran and J. A. Rossmanith},
        Title = {A wave-propagation method for conservation laws and balance laws
                 with spatially varying flux functions},
        Journal = {SIAM J. Sci. Comput.},
        Pages = {955--978},
        Volume = {24},
        Year = {2002}
    }

#### GeoClaw

    M. J. Berger, D. L. George, R. J. LeVeque and K. M. Mandli, The GeoClaw 
    software for depth-averaged flows with adaptive refinement, Advances in Water 
    Resources 34 (2011), pp. 1195-1206.

    @article{BergerGeorgeLeVequeMandli11,
             Author = {M. J. Berger and D. L. George and R. J. LeVeque and K. T.  
             Mandli},
             Journal = {Adv. Water Res.},
             Pages = {1195-1206},
             Title = {The {GeoClaw} software for depth-averaged flows with 
                      adaptive refinement},
             Volume = {34},
             Year = {2011},
             Url = {\url{www.clawpack.org/links/papers/awr11}}
    }

    R. J. LeVeque, D. L. George, and M. J. Berger, 2011, Tsunami modelling with 
    adaptively refined finite volume methods, Acta Numerica, pp. 211-289.

    @article{mjb-dg-rjl:actanum2011,
             Author = {R.J. LeVeque  and D. L. George and M. J. Berger},
             Title = {Adaptive Mesh Refinement Techniques for Tsunamis and Other
                     Geophysical Flows Over Topography},
             Journal = {Acta Numerica},
             Pages = {211-289},
             Year = {2011}
    }

#### PyClaw

Please change the version number and year to the version you have used.

    @misc{pyclaw,
          title={PyClaw software}, 
          url={http://www.pyclaw.org}, 
          author={Mandli, Kyle T. and Ketcheson, David I. and others}, 
          note={Version 5.1}
          year={2014}
    }

    David I. Ketcheson, Kyle T. Mandli, Aron J. Ahmadia, Amal Alghamdi, Manuel 
    Quezada de Luna, Matteo Parsani, Matthew G. Knepley, and Matthew Emmett, 
    2012, PyClaw: Accessible, Extensible, Scalable Tools for Wave Propagation 
    Problems, SIAM Journal on Scientific Computing, 34(4):C210-C231

    @article{pyclaw-sisc,
             Author = {Ketcheson, David I. and Mandli, Kyle T. and Ahmadia, Aron 
             J. and Alghamdi, Amal and {Quezada de Luna}, Manuel and Parsani, 
             Matteo and Knepley, Matthew G. and Emmett, Matthew},
             Journal = {SIAM Journal on Scientific Computing},
             Month = nov,
             Number = {4},
             Pages = {C210--C231},
             Title = {{PyClaw: Accessible, Extensible, Scalable Tools for Wave    
             Propagation Problems}},
             Volume = {34},
             Year = {2012}
    }

#### SharpClaw (High Order WENO)

    D. I. Ketcheson, Matteo Parsani, and R J LeVeque, 2013, High-order Wave 
    Propagation Algorithms for Hyperbolic Systems, SIAM Journal on Scientific 
    Computing, 35(1):A351-A377 (2013)

    @article{Ketcheson2011,
             Author = {Ketcheson, David I. and Parsani, Matteo and LeVeque,
             Randall J.},
             Journal = {SIAM Journal on Scientific Computing},
             Number = {1},
             Pages = {A351--A377},
             Title = {{High-order Wave Propagation Algorithms for Hyperbolic  
                       Systems}},
             Volume = {35},
             Year = {2013}
    }
