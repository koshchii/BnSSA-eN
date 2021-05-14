! Modified by Yongguang Liang to include the R1998 fitting parameters
! 
!File R1990.FORTRN.                                                       
!Reference:  L.W.Whitlow, SLAC-Report-357,                                      
!            Ph.D. Thesis, Stanford University,                                 
!            March 1990.                                                        
!For details see file HELP.DOCUMENT.                                            
                                                                                
!Program contains 135 lines of Fortran code, of 72 characters each, with        
!one subroutine.                                                                
                                                                                
                                                                                
      SUBROUTINE R1998(X,Q2,R,DR,GOODFIT)                                       
                                                                                
! Model for R, based on a fit to world R measurements. Fit performed by         
! program RFIT8 in pseudo-gaussian variable: log(1+.5R).  For details           
! see Reference.                                                                
!                                                                               
! Three models are used, each model has three free parameters.  The             
! functional forms of the models are phenomenological and somewhat              
! contrived.  Each model fits the data very well, and the average of            
! the fits is returned.  The standard deviation of the fit values is            
! used to estimate the systematic uncertainty due to model dependence.          
!                                                                               
! Statistical uncertainties due to fluctuations in measured values have         
! have been studied extensively.  A parametrization of the statistical          
! uncertainty of R1990 is presented in FUNCTION DR1990.                         
!                                                                               
! The three model forms are given by:                                           
!                                                                               
!     R_A = A(1)/LOG(Q2/.04)*FAC + A(2)/[Q2ã4+A(3)ã4]ã.25 ;                     
!     R_B = B(1)/LOG(Q2/.04)*FAC + B(2)/Q2 + B(3)/(Q2**2+.3**2) ;               
!     R_C = C(1)/LOG(Q2/.04)*FAC + C(2)/[(Q2-Q2thr)ã2+C(3)ã2]ã.5 ,              
!                               ...where Q2thr = 5(1-X)ã5 ;                     
!           where FAC = 1+12[Q2/(1+Q2)][.125ã2/(.125ã2+xã2)] gives the          
!           x-dependence of the logarithmic part in order to match Rqcd         
!           at high Q2.                                                         
!                                                                               
! Each model fits very well.  As each model has its own strong points           
! and drawbacks, R1990 returns the average of the models.  The                  
! chisquare for each fit (124 degrees of freedom) are:                          
!     R_A: 110,    R_B: 110,    R_C: 114,    R1990(=avg): 108                   
!                                                                               
! This subroutine returns reasonable values for R for all x and for all         
! Q2 greater than or equal to .3 GeV.                                           
!                                                                               
! The uncertainty in R originates in three sources:                             
!                                                                               
!     D1 = uncertainty in R due to statistical fluctuations of the data         
!          and is parameterized in FUNCTION DR1990, for details see             
!          Reference.                                                           
!                                                                               
!     D2 = uncertainty in R due to possible model dependence, approxi-          
!          mated by the variance between the models.                            
!                                                                               
!     D3 = uncertainty in R due to possible epsilon dependent errors            
!          in the radiative corrections, taken to be +/- .025.  See             
!          theses (mine or Dasu's) for details.                                 
!                                                                               
! and the total error is returned by the program:                               
!                                                                               
!     DR = is the total uncertainty in R, DR = sqrt(D1ã2+D2ã2+D3ã2).            
!          DR is my best estimate of how well we have measured R.  At           
!          high Q2, where R is small, DR is typically larger than R.  If        
!          you have faith in QCD, then, since R1990 = Rqcd at high Q2,          
!          you might wish to assume DR = 0 at very high Q2.                     
!                                                                               
! NOTE:    In many applications, for example the extraction of F2 from          
!          measured cross section, you do not want the full error in R          
!          given by DR.  Rather, you will want to use only the D1 and D2        
!          contributions, and the D3 contribution from radiative                
!          corrections propogates complexely into F2.  For more informa-        
!          tion, see the documentation to dFRC in HELP.DOCUMENT, or             
!          for explicite detail, see Reference.                                 
!                                                                               
!    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        
      IMPLICIT NONE                                                             
      REAL*8 FAC,RLOG,Q2THR,R_A,R_B,R_C,R, D1,D2,D3,DR,DR1990,X,Q2,W2                
      REAL*8 A(3), B(3), C(3), MP,MP2
      LOGICAL GOODFIT                                                           
!    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        
      DATA A/ .06723, .46714, 1.89794 /                                  
      DATA B/ .06347, .57468, -.35342 /                                   
      DATA C/ .05992, .50885, 2.10807 /

      MP = .9382727
      MP2 = MP*MP                                                                        

      FAC   = 1+12.*(Q2/(1.+Q2))*(.125**2/(X**2+.125**2))                       
      RLOG  = FAC/LOG(Q2/.04)!   <--- we use natural logarithms only!           
c      Q2thr = 5.*(1.-X)**5                                                      
                                                                                
c      R_A   = A(1)*RLOG + A(2)/SQRT(SQRT(Q2**4+A(3)**4))                        
c      R_B   = B(1)*RLOG + B(2)/Q2 + B(3)/(Q2**2+.3**2)                          
c      R_C   = C(1)*RLOG + C(2)/SQRT((Q2-Q2thr)**2+C(3)**2)                      


       W2 = Q2/X-Q2+MP2

c       IF (W2 .GT. 3.9 ) THEN 
       q2thr = 12.3708*x-43.1043*x**2  !! r1998
     & +41.7415*x**3

      r_a = 0.0485*rlog
     &+0.5470/sqrt(sqrt(q2**4+2.0621**4))
     &*(1.-0.3804*x+0.5090*x**2)*x**(-0.0285) 
      r_b = 0.0481*rlog
     &+(0.6114/q2-0.3509/(q2**2+.3**2))
     &*(1.-0.4611*x+0.7172*x**2)*x**(-0.0317)                    
      r_c = 0.0577*rlog
     & + 0.4644/sqrt((q2-q2thr)**2+1.8288**2)
c      ELSE
c       q2thr = -5.8119*x+0.9259*x**2  !! r1990
c     & +9.4095*x**3

c      r_a = 0.2839*rlog
c     &+0.0789/sqrt(sqrt(q2**4+0.0**4))
c     &*(1.+32.2507*x-28.5938*x**2)*x**(2.8009) 
c      r_b = 0.2518*rlog
c     &+(0.0110/q2-0.0069/(q2**2+.3**2))
c     &*(1.+57.0509*x-1.4183*x**2)*x**(-0.0934)                    
c      r_c = 0.0621*rlog
c     & + 1.4431/sqrt((q2-q2thr)**2+7.1946**2)
c       write(*,*) r_a, r_b, r_c 
c      ENDIF
      
      R     = (R_A+R_B+R_C)/3.                                                  
                                                                                
      D1    = DR1990(X,Q2)                                                      
      D2    = SQRT(((R_A-R)**2+(R_B-R)**2+(R_C-R)**2)/2.)                       
      D3    = .023*(1.+.5*R)                                                    
              IF (Q2.LT.1.OR.X.LT..1) D3 = 1.5*D3                               
      DR    = SQRT(D1**2+D2**2+D3**2)                                           
                                                                                
      GOODFIT = .TRUE.                                                          
      IF (Q2.LT..3) GOODFIT = .FALSE.                                           
      RETURN                                                                    
      END                                                                       
                                                                                
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        
                                                                                
      FUNCTION DR1990(X,Q2)                                                     
                                                                                
! Parameterizes the uncertainty in R1990 due to the statistical                 
! fluctuations in the data.  Values reflect an average of the R-values          
! about a neighborhood of the specific (x,Q2) value.  That neighborhood         
! is of size [+/-.05] in x, and [+/-33%] in Q2.  For details, see               
! Reference.                                                                    
!                                                                               
! This subroutine is accurate over all (x,Q2), not only the SLAC deep           
! inelastic range.  Where there is no data, for example in the resonance        
! region, it returns a realistic uncertainty, extrapolated from the deep        
! inelastic region (suitably enlarged).  We similarly estimate the              
! uncertainty at very large Q2 by extrapolating from the highest Q2             
! measurments.  For extremely large Q2, R is expected to fall to zero,          
! so the uncertainty in R should not continue to grow.  For this reason         
! DR1990 uses the value at 64 GeV for all larger Q2.                            
!                                                                               
! XHIGH accounts for the rapidly diminishing statistical accuracy for           
! x>.8, and does not contribute for smaller x.                                  
                                                                                
                                                                                
      IMPLICIT NONE                                                             
      REAL*8 U(10,10),DR1990,QMAX,Q,S,A,XLOW,XHIGH,X,Q2                           
                                                                                
                                                                                
      QMAX = 64.                                                                
                                                                                
      Q = MIN(Q2,QMAX)                                                          
      S = .006+.03*X**2                                                         
      A = MAX(.05,8.33*X-.66)                                                   
                                                                                
      XLOW  = .020+ABS(S*LOG(Q/A))                                              
      XHIGH = .1*MAX(.1,X)**20/(.86**20+MAX(.1,X)**20)                          
                                                                                
      DR1990 = SQRT(XLOW**2+XHIGH**2)                                           
      RETURN                                                                    
      END                                                                       
!                                                                               
!End of file R1990.FORTRN.  135 Fortran lines.                                  








