C  Important Notice:
C   These algorithms are modifications and based on the Software Newuoa, authored by M. J. D. Powell,
C   to minimize sum of squares by taking advantage of the problem structure.
C    min_{x \in R^n }  F(x) := Sum_{i=1}^{mv}  v_err_i(x)^2
C   where v_err(x) : R^n \to R^{mv} is a vector function.
C   This subroutine seeks the least value of sum of the squres of the components of v_err(x)
C   by combing trust region method and Levenberg-Marquardt method 
C
C   References:
C
C   1.  M. J. D. Powell, The NEWUOA software for unconstrained optimization without derivatives,
C       DAMTP 2004/ NA 05
C   2.  H. Zhang, A. R. CONN, AND K. SCHEINBERG, A DERIVATIVE-FREE ALGORITHM FOR THE LEAST-SQUARE
C       MINIMIZATION, technical report, 2009
C
C      -----------------------------------------------------------------
C      | This program is free software; you can redistribute it and/or  |
C      |modify it under the terms of the GNU General Public License as  |
C      |published by the Free Software Foundation; either version 2 of  |
C      |the License, or (at your option) any later version.             |
C      |This program is distributed in the hope that it will be useful, |
C      |but WITHOUT ANY WARRANTY; without even the implied warranty of  |
C      |MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   |
C      |GNU General Public License for more details.                    |
C      |                                                                |
C      |You should have received a copy of the GNU General Public       |
C      |License along with this program; if not, write to the Free      |
C      |Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, |
C      |MA  02110-1301  USA                                             |
C      -----------------------------------------------------------------|

      SUBROUTINE NEWUOA_H(N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W,mv)
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION X(*),W(*)
C
C     N must be set to the number of variables and must be at least two.
C     mv must be set to the lengh of the vector function v_err(x):  R^n \to R^{mv}. 
C     The maximum number variables in this codes: nmax =  100
C     The maximum lengh of the vector function v_err(x): mmax = 400
C     If n > 100 or m > 400, the parameter nmax and mmax need to be creased in
C     subroutine  NEWUOB_H and TRSAPP_H
C
C     NPT is the number of interpolation conditions. Its value must be in the
C     interval [N+2,2N+1].  Recommended: NPT = 2*N+1
C
C     Initial values of the variables must be set in X(1),X(2),...,X(N). They
C     will be changed to the values that give the least calculated F(x) = Sum_{i=1}^{mv} v_err_i(x)^2.
C
C     RHOBEG and RHOEND must be set to the initial and final values of a trust
C     region radius, so both must be positive with RHOEND<=RHOBEG. Typically
C     RHOBEG should be about one tenth of the greatest expected change to a
C     variable, and RHOEND should indicate the accuracy that is required in
C     the final values of the variables. Default: RHOBEG = 1.0, RHOEND = 10^{-8} 
C
C     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
C     amount of printing. Specifically, there is no output if IPRINT=0 and
C     there is output only at the return if IPRINT=1. Otherwise, each new
C     value of RHO is printed, with the best vector of variables so far and
C     the corresponding value of the objective function. Further, each new
C     value of F with its variables are output if IPRINT=3.
C
C     MAXFUN must be set to an upper bound on the number of calls of subroutine
C     dfovec(n, mv, x, v_err) which provides the values of the vector function v_err(x).
C     Here: n, mv, x \in R^n are input, v_err \in R^{mv} are output.
C     Default:  MAXFUN= 100(n+1), i.e 100 (simplex) gradients for reasonable accuracy.
C               MAXFUN= infinity, to let the algorithm explore the lowest function value  
C                       as much as it could.
C
C     The array W will be used for working space. Its length must be at least
C     (NPT+11)*(NPT+N)+N*(3*N+11)/2 
C
C     SUBROUTINE dfovec(n, mv, x, v_err) must be provided by the user. It must provide
C     the values of the vector function v_err(x) : R^n to R^{mv} at the variables X(1),X(2),...,X(N).
C
      NP=N+1
      NPTM=NPT-NP
      IF (NPT .LT. N+2 .OR. NPT .GT. 2*N+1) THEN
          PRINT 10
   10     FORMAT (/4X,'Return from NEWUOA because NPT is not in',
     1      '[N+2, 2N+1]')
          GO TO 20
      END IF
      NDIM=NPT+N
      IXB=1
      IXO=IXB+N
      IXN=IXO+N
      IXP=IXN+N
      IGQ=IXP+N*NPT
      IHQ=IGQ+N
      IPQ=IHQ+(N*NP)/2
      IBMAT=IPQ+NPT
      IZMAT=IBMAT+NDIM*N
      ID=IZMAT+NPT*NPTM
      IVL=ID+N
      IW=IVL+NDIM
C
C     The above settings provide a partition of W for subroutine NEWUOB_H.
C
      CALL NEWUOB_H (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W(IXB),
     1  W(IXO),W(IXN),W(IXP),W(IGQ),W(IHQ),W(IPQ),W(IBMAT),W(IZMAT),
     2  NDIM,W(ID),W(IVL),W(IW),mv)
   20 RETURN
      END

