        program test

        IMPLICIT COMPLEX (X)
        DOUBLE PRECISION b, c, s, cp, ff, res


	XI = (0, 1)
	XR = (1, 0)

        b=6.0
        
        s=7000**2
        pi = 4.0*atan(1.0)
        c = .167
        cp = .748
        XSs = s**c/(log(s)**cp)+s**c*exp(-XI*Pi*c)/((log(s)-XI*Pi)**cp)
        FF = 1.031925910*b*besk1(0.577*b)-.2119478974*besk0(0.577*b)
     $  +3.678319512*b*besk1(1.719*b)-17.76215332*besk0(1.719*b)
     $  +19.36687661*besk0(1.858*b)
	
        res = realpart(XR - exp(-XSs*FF))
        
        WRITE(*, *) res
        
          
         end program test
          