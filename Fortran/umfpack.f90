! ======================================================================
! UMFPACK Fortran interface via the mUMFPACK module
! ======================================================================
! Version 1.0 (Apr 23, 2014) for UMFPACK version 5.6.2
! ======================================================================
! Compile with any Fortran compiler with support of iso_c_binding module
! and link with the UMFPACK C library:
! e.g., gfortran umfpack.f90 my_code.f90 -lumfpack
! ======================================================================
module mUMFPACK
! ======================================================================
use iso_c_binding
implicit none

! private size constants

integer,parameter,private :: i4=4, & ! size of default integer
                             i8=8, & ! size of long integer
                    ip=c_intptr_t, & ! size of pointers used in basic Fortran wrappers
                             r4=4, & ! size of single precision real/complex
                             r8=8    ! size of double precision real/complex

! default pointers to UMFPACK Symbolic and Numeric objects

type(c_ptr),private :: iSymbolic=c_null_ptr,iNumeric=c_null_ptr

! ======================================================================
! UMFPACK constants
! ======================================================================

! size of Info and Control arrays

integer,parameter :: UMFPACK_INFO=90, &
                     UMFPACK_CONTROL=20

! Version, copyright, and license

character(30),parameter :: UMFPACK_VERSION="UMFPACK V5.6.2 (Apr 25, 2013)"
character(79),parameter :: UMFPACK_COPYRIGHT="UMFPACK:  Copyright (c) 2005-2012 by Timothy A. Davis.  All Rights Reserved."
character(12),parameter :: UMFPACK_DATE="Apr 25, 2013"
integer,parameter :: UMFPACK_MAIN_VERSION=5, &
                     UMFPACK_SUB_VERSION=6, &
                     UMFPACK_SUBSUB_VERSION=2
! UMFPACK_VER_CODE(main,sub) ((main) * 1000 + (sub))
! UMFPACK_VER UMFPACK_VER_CODE(UMFPACK_MAIN_VERSION,UMFPACK_SUB_VERSION)

! contents of Info

enum,bind(c) 
enumerator :: & 
! returned by all routines that use Info:
UMFPACK_STATUS=0, &          ! /* UMFPACK_OK, or other result */
UMFPACK_NROW=1, &            ! /* n_row input value */
UMFPACK_NCOL=16, &           ! /* n_col input value */
UMFPACK_NZ=2, &              ! /* # of entries in A */
! computed in UMFPACK_*symbolic and UMFPACK_numeric:
UMFPACK_SIZE_OF_UNIT=3, &    ! /* sizeof (Unit) */
! computed in UMFPACK_*symbolic:
UMFPACK_SIZE_OF_INT=4, &     ! /* sizeof (int) */
UMFPACK_SIZE_OF_LONG=5, &    ! /* sizeof (SuiteSparse_long) */
UMFPACK_SIZE_OF_POINTER=6, & ! /* sizeof (void *) */
UMFPACK_SIZE_OF_ENTRY=7, &   ! /* sizeof (Entry), real or complex */
UMFPACK_NDENSE_ROW=8, &      ! /* number of dense rows */
UMFPACK_NEMPTY_ROW=9, &      ! /* number of empty rows */
UMFPACK_NDENSE_COL=10, &     ! /* number of dense rows */
UMFPACK_NEMPTY_COL=11, &     ! /* number of empty rows */
UMFPACK_SYMBOLIC_DEFRAG=12, & ! /* # of memory compactions */
UMFPACK_SYMBOLIC_PEAK_MEMORY=13, & ! /* memory used by symbolic analysis */
UMFPACK_SYMBOLIC_SIZE=14, &  ! /* size of Symbolic object, in Units */
UMFPACK_SYMBOLIC_TIME=15, &  ! /* time (sec.) for symbolic analysis */
UMFPACK_SYMBOLIC_WALLTIME=17, & ! /* wall clock time for sym. analysis */
UMFPACK_STRATEGY_USED=18, &  ! /* strategy used: sym, unsym */
UMFPACK_ORDERING_USED=19, &  ! /* ordering used: colamd, amd, given */
UMFPACK_QFIXED=31, &         ! /* whether Q is fixed or refined */
UMFPACK_DIAG_PREFERRED=32, & ! /* whether diagonal pivoting attempted*/
UMFPACK_PATTERN_SYMMETRY=33, & ! /* symmetry of pattern of S */
UMFPACK_NZ_A_PLUS_AT=34, &   ! /* nnz (S+S'), excl. diagonal */
UMFPACK_NZDIAG=35, &         ! /* nnz (diag (S)) */
! AMD statistics, computed in UMFPACK_*symbolic:
UMFPACK_SYMMETRIC_LUNZ=36, & ! /* nz in L+U, if AMD ordering used */
UMFPACK_SYMMETRIC_FLOPS=37, & ! /* flops for LU, if AMD ordering used */
UMFPACK_SYMMETRIC_NDENSE=38, & ! /* # of "dense" rows/cols in S+S' */
UMFPACK_SYMMETRIC_DMAX=39, & ! /* max nz in cols of L, for AMD */
! 51:55 unused
! statistics for singleton pruning
UMFPACK_COL_SINGLETONS=56, & ! /* # of column singletons */
UMFPACK_ROW_SINGLETONS=57, & ! /* # of row singletons */
UMFPACK_N2=58, &             ! /* size of S */
UMFPACK_S_SYMMETRIC=59, &    ! /* 1 if S square and symmetricly perm.*/
! estimates computed in UMFPACK_*symbolic:
UMFPACK_NUMERIC_SIZE_ESTIMATE=20, & ! /* final size of Numeric->Memory */
UMFPACK_PEAK_MEMORY_ESTIMATE=21, & ! /* for symbolic & numeric */
UMFPACK_FLOPS_ESTIMATE=22, & ! /* flop count */
UMFPACK_LNZ_ESTIMATE=23, &   ! /* nz in L, incl. diagonal */
UMFPACK_UNZ_ESTIMATE=24, &   ! /* nz in U, incl. diagonal */
UMFPACK_VARIABLE_INIT_ESTIMATE=25, & ! /* initial size of Numeric->Memory*/
UMFPACK_VARIABLE_PEAK_ESTIMATE=26, & ! /* peak size of Numeric->Memory */
UMFPACK_VARIABLE_FINAL_ESTIMATE=27, & ! /* final size of Numeric->Memory */
UMFPACK_MAX_FRONT_SIZE_ESTIMATE=28, & ! /* max frontal matrix size */
UMFPACK_MAX_FRONT_NROWS_ESTIMATE=29, & ! /* max # rows in any front */
UMFPACK_MAX_FRONT_NCOLS_ESTIMATE=30, & ! /* max # columns in any front */
! exact values, (estimates shown above) computed in UMFPACK_numeric:
UMFPACK_NUMERIC_SIZE=40, &   ! /* final size of Numeric->Memory */
UMFPACK_PEAK_MEMORY=41, &    ! /* for symbolic & numeric */
UMFPACK_FLOPS=42, &          ! /* flop count */
UMFPACK_LNZ=43, &            ! /* nz in L, incl. diagonal */
UMFPACK_UNZ=44, &            ! /* nz in U, incl. diagonal */
UMFPACK_VARIABLE_INIT=45, &  ! /* initial size of Numeric->Memory*/
UMFPACK_VARIABLE_PEAK=46, &  ! /* peak size of Numeric->Memory */
UMFPACK_VARIABLE_FINAL=47, & ! /* final size of Numeric->Memory */
UMFPACK_MAX_FRONT_SIZE=48, & ! /* max frontal matrix size */
UMFPACK_MAX_FRONT_NROWS=49, & ! /* max # rows in any front */
UMFPACK_MAX_FRONT_NCOLS=50, & ! /* max # columns in any front */
! computed in UMFPACK_numeric:
UMFPACK_NUMERIC_DEFRAG=60, & ! /* # of garbage collections */
UMFPACK_NUMERIC_REALLOC=61, & ! /* # of memory reallocations */
UMFPACK_NUMERIC_COSTLY_REALLOC=62, & ! /* # of costlly memory realloc's */
UMFPACK_COMPRESSED_PATTERN=63, & ! /* # of integers in LU pattern */
UMFPACK_LU_ENTRIES=64, &     ! /* # of reals in LU factors */
UMFPACK_NUMERIC_TIME=65, &   ! /* numeric factorization time */
UMFPACK_UDIAG_NZ=66, &       ! /* nz on diagonal of U */
UMFPACK_RCOND=67, &          ! /* est. reciprocal condition # */
UMFPACK_WAS_SCALED=68, &     ! /* none, max row, or sum row */
UMFPACK_RSMIN=69, &          ! /* min (max row) or min (sum row) */
UMFPACK_RSMAX=70, &          ! /* max (max row) or max (sum row) */
UMFPACK_UMIN=71, &           ! /* min abs diagonal entry of U */
UMFPACK_UMAX=72, &           ! /* max abs diagonal entry of U */
UMFPACK_ALLOC_INIT_USED=73, & ! /* alloc_init parameter used */
UMFPACK_FORCED_UPDATES=74, & ! /* # of forced updates */
UMFPACK_NUMERIC_WALLTIME=75, & ! /* numeric wall clock time */
UMFPACK_NOFF_DIAG=76, &      ! /* number of off-diagonal pivots */
UMFPACK_ALL_LNZ=77, &        ! /* nz in L, if no dropped entries */
UMFPACK_ALL_UNZ=78, &        ! /* nz in U, if no dropped entries */
UMFPACK_NZDROPPED=79, &      ! /* # of dropped small entries */
! computed in UMFPACK_solve:
UMFPACK_IR_TAKEN=80, &       ! /* # of iterative refinement steps taken */
UMFPACK_IR_ATTEMPTED=81, &   ! /* # of iter. refinement steps attempted */
UMFPACK_OMEGA1=82, &         ! /* omega1, sparse backward error estimate */
UMFPACK_OMEGA2=83, &         ! /* omega2, sparse backward error estimate */
UMFPACK_SOLVE_FLOPS=84, &    ! /* flop count for solve */
UMFPACK_SOLVE_TIME=85, &     ! /* solve time (seconds) */
UMFPACK_SOLVE_WALLTIME=86    ! /* solve time (wall clock, seconds) */
! Info(87,88,89) unused
! Unused parts of Info may be used in future versions of UMFPACK.
end enum

! contents of Control

enum,bind(c) 
enumerator :: & 
! used in all UMFPACK_report_* routines:
UMFPACK_PRL=0, &             ! /* print level */
! used in UMFPACK_*symbolic only:
UMFPACK_DENSE_ROW=1, &       ! /* dense row parameter */
UMFPACK_DENSE_COL=2, &       ! /* dense col parameter */
UMFPACK_BLOCK_SIZE=4, &      ! /* BLAS-3 block size */
UMFPACK_STRATEGY=5, &        ! /* auto, symmetric, or unsym. */
UMFPACK_ORDERING=10, &       ! /* ordering method to use */
UMFPACK_FIXQ=13, &           ! /* -1: no fixQ, 0: default, 1: fixQ */
UMFPACK_AMD_DENSE=14, &      ! /* for AMD ordering */
UMFPACK_AGGRESSIVE=19, &     ! /* whether or not to use aggressive */
UMFPACK_SINGLETONS=11, &     ! /* singleton filter on if true */
! used in UMFPACK_*numeric only:
UMFPACK_PIVOT_TOLERANCE=3, & ! /* threshold partial pivoting setting */
UMFPACK_ALLOC_INIT=6, &      ! /* initial allocation ratio */
UMFPACK_SYM_PIVOT_TOLERANCE=15, & ! /* threshold, only for diag. entries */
UMFPACK_SCALE=16, &          ! /* what row scaling to do */
UMFPACK_FRONT_ALLOC_INIT=17, & ! /* frontal matrix allocation ratio */
UMFPACK_DROPTOL=18, &        ! /* drop tolerance for entries in L,U */
! used in UMFPACK_*solve only:
UMFPACK_IRSTEP=7, &          ! /* max # of iterative refinements */
! compile-time settings - Control(8:11) cannot be changed at run time:
UMFPACK_COMPILED_WITH_BLAS=8 ! /* uses the BLAS */
end enum                     ! /* 9,12: unused */

! Control(UMFPACK_STRATEGY) is one of the following:
enum,bind(c) ; enumerator :: & 
UMFPACK_STRATEGY_AUTO=0, &   ! /* use sym. or unsym. strategy */
UMFPACK_STRATEGY_UNSYMMETRIC, & ! /* COLAMD(A), coletree postorder, not prefer diag*/
UMFPACK_STRATEGY_OBSOLETE, & ! /* 2-by-2 is no longer available */
UMFPACK_STRATEGY_SYMMETRIC   ! /* AMD(A+A'), no coletree postorder, prefer diagonal */
end enum

! Control(UMFPACK_SCALE) is one of the following:
enum,bind(c) ; enumerator :: & 
UMFPACK_SCALE_NONE=0, &      ! /* no scaling */
UMFPACK_SCALE_SUM, &         ! /* default: divide each row by sum (abs (row))*/
UMFPACK_SCALE_MAX            ! /* divide each row by max (abs (row)) */
end enum

! Control(UMFPACK_ORDERING) and Info(UMFPACK_ORDERING_USED) are one of:
enum,bind(c) ; enumerator :: & 
UMFPACK_ORDERING_CHOLMOD=0, & ! /* use CHOLMOD (AMD/COLAMD then METIS)*/
UMFPACK_ORDERING_AMD, &      ! /* use AMD/COLAMD */
UMFPACK_ORDERING_GIVEN, &    ! /* user-provided Qinit */
UMFPACK_ORDERING_METIS, &    ! /* use METIS */
UMFPACK_ORDERING_BEST, &     ! /* try many orderings, pick best */
UMFPACK_ORDERING_NONE, &     ! /* natural ordering */
UMFPACK_ORDERING_USER        ! /* user-provided function */
! /* AMD/COLAMD means: use AMD for symmetric strategy, COLAMD for unsymmetric */
end enum

! status codes

enum,bind(c) ; enumerator :: & 
UMFPACK_OK=0, &
! /* status > 0 means a warning, but the method was successful anyway. */
! /* A Symbolic or Numeric object was still created. */
UMFPACK_WARNING_singular_matrix=1, &
! /* The following warnings were added in umfpack_*_get_determinant */
UMFPACK_WARNING_determinant_underflow=2, &
UMFPACK_WARNING_determinant_overflow=3, &
! /* status < 0 means an error, and the method was not successful. */
! /* No Symbolic of Numeric object was created. */
UMFPACK_ERROR_out_of_memory=-1, &
UMFPACK_ERROR_invalid_Numeric_object=-3, &
UMFPACK_ERROR_invalid_Symbolic_object=-4, &
UMFPACK_ERROR_argument_missing=-5, &
UMFPACK_ERROR_n_nonpositive=-6, &
UMFPACK_ERROR_invalid_matrix=-8, &
UMFPACK_ERROR_different_pattern=-11, &
UMFPACK_ERROR_invalid_system=-13, &
UMFPACK_ERROR_invalid_permutation=-15, &
UMFPACK_ERROR_internal_error=-911, & ! /* yes, call me if you get this! */
UMFPACK_ERROR_file_IO=-17, &
UMFPACK_ERROR_ordering_failed=-18
end enum

! solve codes
! /* Solve the system ( )x=b, where ( ) is defined below.  "t" refers to the */
! /* linear algebraic transpose (complex conjugate if A is complex), or the (') */
! /* operator in MATLAB.  "at" refers to the array transpose, or the (.') */
! /* operator in MATLAB. */

enum,bind(c) ; enumerator :: & 
UMFPACK_A=0, &      ! /* Ax=b    */
UMFPACK_At=1, &     ! /* A'x=b   */
UMFPACK_Aat=2, &    ! /* A.'x=b  */
UMFPACK_Pt_L=3, &   ! /* P'Lx=b  */
UMFPACK_L=4, &      ! /* Lx=b    */
UMFPACK_Lt_P=5, &   ! /* L'Px=b  */
UMFPACK_Lat_P=6, &  ! /* L.'Px=b */
UMFPACK_Lt=7, &     ! /* L'x=b   */
UMFPACK_Lat=8, &    ! /* L.'x=b  */
UMFPACK_U_Qt=9, &   ! /* UQ'x=b  */
UMFPACK_U=10, &     ! /* Ux=b    */
UMFPACK_Q_Ut=11, &  ! /* QU'x=b  */
UMFPACK_Q_Uat=12, & ! /* QU.'x=b */
UMFPACK_Ut=13, &    ! /* U'x=b   */
UMFPACK_Uat=14      ! /* U.'x=b  */
end enum

! ======================================================================
! Full UMFPACK interface: interfaces to C functions
! ======================================================================

interface

! int umfpack_di_symbolic(int n_row,int n_col,const int Ap [ ],const int Ai [ ],const double Ax [ ],
! void **Symbolic,const double Control [UMFPACK_CONTROL],double Info [UMFPACK_INFO]) ;

integer(c_int) function c_umfpack_di_symbolic(n_row,n_col,Ap,Ai,Ax,Symbolic,Control,Info) bind(c,name='umfpack_di_symbolic')
import c_int,c_ptr
integer(c_int),value :: n_row,n_col
type(c_ptr),value,intent(in) :: Ap,Ai
type(c_ptr),value,intent(in) :: Ax
type(c_ptr) :: Symbolic
type(c_ptr),value :: Control,Info
end function

integer(c_int) function c_umfpack_zi_symbolic(n_row,n_col,Ap,Ai,Ax,Az,Symbolic,Control,Info) bind(c,name='umfpack_zi_symbolic')
import c_int,c_ptr
integer(c_int),value :: n_row,n_col
type(c_ptr),value,intent(in) :: Ap,Ai
type(c_ptr),value,intent(in) :: Ax,Az
type(c_ptr) :: Symbolic
type(c_ptr),value :: Control,Info
end function

! int umfpack_di_numeric(const int Ap [ ],const int Ai [ ],const double Ax [ ],
! void *Symbolic,void **Numeric,const double Control [UMFPACK_CONTROL],double Info [UMFPACK_INFO]) ;

integer(c_int) function c_umfpack_di_numeric(Ap,Ai,Ax,Symbolic,Numeric,Control,Info) bind(c,name='umfpack_di_numeric')
import c_int,c_ptr
type(c_ptr),value,intent(in) :: Ap,Ai
type(c_ptr),value,intent(in) :: Ax
type(c_ptr),value,intent(in) :: Symbolic
type(c_ptr) :: Numeric
type(c_ptr),value :: Control,Info
end function

integer(c_int) function c_umfpack_zi_numeric(Ap,Ai,Ax,Az,Symbolic,Numeric,Control,Info) bind(c,name='umfpack_zi_numeric')
import c_int,c_ptr
type(c_ptr),value,intent(in) :: Ap,Ai
type(c_ptr),value,intent(in) :: Ax,Az
type(c_ptr),value,intent(in) :: Symbolic
type(c_ptr) :: Numeric
type(c_ptr),value :: Control,Info
end function

! int umfpack_di_solve(int sys,const int Ap [ ],const int Ai [ ],const double Ax [ ],double X [ ],const double B [ ],
! void *Numeric,const double Control [UMFPACK_CONTROL],double Info [UMFPACK_INFO]) ;

integer(c_int) function c_umfpack_di_solve(sys,Ap,Ai,Ax,X,B,Numeric,Control,Info) bind(c,name='umfpack_di_solve')
import c_int,c_ptr
integer(c_int),value :: sys
type(c_ptr),value,intent(in) :: Ap,Ai
type(c_ptr),value,intent(in) :: Ax
type(c_ptr),value :: X
type(c_ptr),value,intent(in) :: B
type(c_ptr),value,intent(in) :: Numeric
type(c_ptr),value :: Control,Info
end function

integer(c_int) function c_umfpack_zi_solve(sys,Ap,Ai,Ax,Az,Xx,Xz,Bx,Bz,Numeric,Control,Info) bind(c,name='umfpack_zi_solve')
import c_int,c_ptr
integer(c_int),value :: sys
type(c_ptr),value,intent(in) :: Ap,Ai
type(c_ptr),value,intent(in) :: Ax,Az
type(c_ptr),value :: Xx,Xz
type(c_ptr),value,intent(in) :: Bx,Bz
type(c_ptr),value,intent(in) :: Numeric
type(c_ptr),value :: Control,Info
end function

! void umfpack_di_free_symbolic(void **Symbolic) ;

subroutine c_umfpack_di_free_symbolic(Symbolic) bind(c,name='umfpack_di_free_symbolic')
import c_ptr
type(c_ptr) :: Symbolic
end subroutine

subroutine c_umfpack_zi_free_symbolic(Symbolic) bind(c,name='umfpack_zi_free_symbolic')
import c_ptr
type(c_ptr) :: Symbolic
end subroutine

! void umfpack_di_free_numeric(void **Numeric) ;

subroutine c_umfpack_di_free_numeric(Numeric) bind(c,name='umfpack_di_free_numeric')
import c_ptr
type(c_ptr) :: Numeric
end subroutine

subroutine c_umfpack_zi_free_numeric(Numeric) bind(c,name='umfpack_zi_free_numeric')
import c_ptr
type(c_ptr) :: Numeric
end subroutine

! void umfpack_di_defaults(double Control [UMFPACK_CONTROL]) ;

subroutine c_umfpack_di_defaults(Control) bind(c,name='umfpack_di_defaults')
import c_ptr
type(c_ptr),value :: Control
end subroutine

subroutine c_umfpack_zi_defaults(Control) bind(c,name='umfpack_zi_defaults')
import c_ptr
type(c_ptr),value :: Control
end subroutine

! int umfpack_di_scale(double X [ ],const double B [ ],void *Numeric) ;

integer(c_int) function c_umfpack_di_scale(X,B,Numeric) bind(c,name='umfpack_di_scale')
import c_int,c_ptr
type(c_ptr),value :: X
type(c_ptr),value,intent(in) :: B
type(c_ptr),value,intent(in) :: Numeric
end function

integer(c_int) function c_umfpack_zi_scale(Xx,Xz,Bx,Bz,Numeric) bind(c,name='umfpack_zi_scale')
import c_int,c_ptr
type(c_ptr),value :: Xx,Xz
type(c_ptr),value,intent(in) :: Bx,Bz
type(c_ptr),value,intent(in) :: Numeric
end function

! int umfpack_di_save_numeric(void *Numeric,char *filename) ;

integer(c_int) function c_umfpack_di_save_numeric(Numeric,filename) bind(c,name='umfpack_di_save_numeric')
import c_int,c_ptr,c_char
type(c_ptr),value :: Numeric
character(1,c_char) :: filename(*)
end function

integer(c_int) function c_umfpack_zi_save_numeric(Numeric,filename) bind(c,name='umfpack_zi_save_numeric')
import c_int,c_ptr,c_char
type(c_ptr),value :: Numeric
character(1,c_char) :: filename(*)
end function

! int umfpack_di_save_symbolic(void *Symbolic,char *filename) ;

integer(c_int) function c_umfpack_di_save_symbolic(Symbolic,filename) bind(c,name='umfpack_di_save_symbolic')
import c_int,c_ptr,c_char
type(c_ptr),value :: Symbolic
character(1,c_char) :: filename(*)
end function

integer(c_int) function c_umfpack_zi_save_symbolic(Symbolic,filename) bind(c,name='umfpack_zi_save_symbolic')
import c_int,c_ptr,c_char
type(c_ptr),value :: Symbolic
character(1,c_char) :: filename(*)
end function

! int umfpack_di_load_numeric(void **Numeric,char *filename) ;

integer(c_int) function c_umfpack_di_load_numeric(Numeric,filename) bind(c,name='umfpack_di_load_numeric')
import c_int,c_ptr,c_char
type(c_ptr) :: Numeric
character(1,c_char) :: filename(*)
end function

integer(c_int) function c_umfpack_zi_load_numeric(Numeric,filename) bind(c,name='umfpack_zi_load_numeric')
import c_int,c_ptr,c_char
type(c_ptr) :: Numeric
character(1,c_char) :: filename(*)
end function

! int umfpack_di_load_symbolic(void **Symbolic,char *filename) ;

integer(c_int) function c_umfpack_di_load_symbolic(Symbolic,filename) bind(c,name='umfpack_di_load_symbolic')
import c_int,c_ptr,c_char
type(c_ptr) :: Symbolic
character(1,c_char) :: filename(*)
end function

integer(c_int) function c_umfpack_zi_load_symbolic(Symbolic,filename) bind(c,name='umfpack_zi_load_symbolic')
import c_int,c_ptr,c_char
type(c_ptr) :: Symbolic
character(1,c_char) :: filename(*)
end function

! void umfpack_di_report_status(const double Control [UMFPACK_CONTROL],int status) ;

subroutine c_umfpack_di_report_status(Control,status) bind(c,name='umfpack_di_report_status')
import c_ptr,c_int
type(c_ptr),value :: Control
integer(c_int),value :: status
end subroutine

subroutine c_umfpack_zi_report_status(Control,status) bind(c,name='umfpack_zi_report_status')
import c_ptr,c_int
type(c_ptr),value :: Control
integer(c_int),value :: status
end subroutine

! void umfpack_di_report_control(const double Control [UMFPACK_CONTROL]) ;

subroutine c_umfpack_di_report_control(Control) bind(c,name='umfpack_di_report_control')
import c_ptr
type(c_ptr),value :: Control
end subroutine

subroutine c_umfpack_zi_report_control(Control) bind(c,name='umfpack_zi_report_control')
import c_ptr
type(c_ptr),value :: Control
end subroutine

! void umfpack_di_report_info(const double Control [UMFPACK_CONTROL],const double Info [UMFPACK_INFO]) ;

subroutine c_umfpack_di_report_info(Control,Info) bind(c,name='umfpack_di_report_info')
import c_ptr
type(c_ptr),value :: Control,Info
end subroutine

subroutine c_umfpack_zi_report_info(Control,Info) bind(c,name='umfpack_zi_report_info')
import c_ptr
type(c_ptr),value :: Control,Info
end subroutine

! int umfpack_di_report_numeric(void *Numeric,const double Control [UMFPACK_CONTROL]) ;

integer(c_int) function c_umfpack_di_report_numeric(Numeric,Control) bind(c,name='umfpack_di_report_numeric')
import c_int,c_ptr
type(c_ptr),value :: Numeric,Control
end function

integer(c_int) function c_umfpack_zi_report_numeric(Numeric,Control) bind(c,name='umfpack_zi_report_numeric')
import c_int,c_ptr
type(c_ptr),value :: Numeric,Control
end function

! int umfpack_di_report_symbolic(void *Symbolic,const double Control [UMFPACK_CONTROL]) ;

integer(c_int) function  c_umfpack_di_report_symbolic(Symbolic,Control) bind(c,name='umfpack_di_report_symbolic')
import c_int,c_ptr
type(c_ptr),value :: Symbolic,Control
end function

integer(c_int) function  c_umfpack_zi_report_symbolic(Symbolic,Control) bind(c,name='umfpack_zi_report_symbolic')
import c_int,c_ptr
type(c_ptr),value :: Symbolic,Control
end function

end interface

! ======================================================================
! Full UMFPACK interface: overloaded names
! ======================================================================

interface umfpack_zi_symbolic
module procedure umfpack_zi_symbolic,umfpack_ci_symbolic
end interface
interface s_umfpack_zi_symbolic
module procedure s_umfpack_zi_symbolic,s_umfpack_ci_symbolic
end interface
interface umfpack_symbolic
module procedure umfpack_di_symbolic,umfpack_zi_symbolic,umfpack_ci_symbolic
end interface
interface s_umfpack_symbolic
module procedure s_umfpack_di_symbolic,s_umfpack_zi_symbolic,s_umfpack_ci_symbolic
end interface

interface umfpack_zi_numeric
module procedure umfpack_zi_numeric,umfpack_ci_numeric
end interface
interface s_umfpack_zi_numeric
module procedure s_umfpack_zi_numeric,s_umfpack_ci_numeric
end interface
interface umfpack_numeric
module procedure umfpack_di_numeric,umfpack_zi_numeric,umfpack_ci_numeric
end interface
interface s_umfpack_numeric
module procedure s_umfpack_di_numeric,s_umfpack_zi_numeric,s_umfpack_ci_numeric
end interface

interface umfpack_zi_solve
module procedure umfpack_zi_solve,umfpack_ci_solve
end interface
interface s_umfpack_zi_solve
module procedure s_umfpack_zi_solve,s_umfpack_ci_solve
end interface
interface umfpack_solve
module procedure umfpack_di_solve,umfpack_zi_solve,umfpack_ci_solve
end interface
interface s_umfpack_solve
module procedure s_umfpack_di_solve,s_umfpack_zi_solve,s_umfpack_ci_solve
end interface

interface s_umfpack_free_symbolic
module procedure umfpack_free_symbolic
end interface
interface s_umfpack_di_free_symbolic
module procedure umfpack_di_free_symbolic
end interface
interface s_umfpack_zi_free_symbolic
module procedure umfpack_zi_free_symbolic
end interface

interface s_umfpack_free_numeric
module procedure umfpack_free_numeric
end interface
interface s_umfpack_di_free_numeric
module procedure umfpack_di_free_numeric
end interface
interface s_umfpack_zi_free_numeric
module procedure umfpack_zi_free_numeric
end interface

interface s_umfpack_defaults
module procedure umfpack_defaults
end interface
interface s_umfpack_di_defaults
module procedure umfpack_di_defaults
end interface
interface s_umfpack_zi_defaults
module procedure umfpack_zi_defaults
end interface

! a conflict with the constant UMFPACK_SCALE
interface umfpack_zi_scale
module procedure umfpack_zi_scale,umfpack_ci_scale
end interface
interface umfpack_scale_function
module procedure umfpack_di_scale,umfpack_zi_scale,umfpack_ci_scale
end interface

interface s_umfpack_zi_scale
module procedure s_umfpack_zi_scale,s_umfpack_ci_scale
end interface
interface s_umfpack_scale
module procedure s_umfpack_di_scale,s_umfpack_zi_scale,s_umfpack_ci_scale
end interface

interface s_umfpack_report_control
module procedure umfpack_report_control
end interface
interface s_umfpack_di_report_control
module procedure umfpack_di_report_control
end interface
interface s_umfpack_zi_report_control
module procedure umfpack_zi_report_control
end interface

interface s_umfpack_report_info
module procedure umfpack_report_info
end interface
interface s_umfpack_di_report_info
module procedure umfpack_di_report_info
end interface
interface s_umfpack_zi_report_info
module procedure umfpack_zi_report_info
end interface

! ======================================================================
! Basic UMFPACK interface: overloaded names
! ======================================================================

interface umf4csym
module procedure umf4csym,umf4csym_ip
end interface
interface umf4zsym
module procedure umf4zsym,umf4zsym_ip,umf4csym,umf4csym_ip
end interface
interface umf4sym
module procedure umf4sym,umf4sym_ip,umf4zsym,umf4zsym_ip,umf4csym,umf4csym_ip
end interface

interface umf4cnum
module procedure umf4cnum,umf4cnum_ip
end interface
interface umf4znum
module procedure umf4znum,umf4znum_ip,umf4cnum,umf4cnum_ip
end interface
interface umf4num
module procedure umf4num,umf4num_ip,umf4znum,umf4znum_ip,umf4cnum,umf4cnum_ip
end interface

interface umf4csolr
module procedure umf4csolr,umf4csolr_ip
end interface
interface umf4zsolr
module procedure umf4zsolr,umf4zsolr_ip,umf4csolr,umf4csolr_ip
end interface
interface umf4solr
module procedure umf4solr,umf4solr_ip,umf4zsolr,umf4zsolr_ip,umf4csolr,umf4csolr_ip
end interface

interface umf4csol
module procedure umf4csol,umf4csol_ip
end interface
interface umf4zsol
module procedure umf4zsol,umf4zsol_ip,umf4csol,umf4csol_ip
end interface
interface umf4sol
module procedure umf4sol,umf4sol_ip,umf4zsol,umf4zsol_ip,umf4csol,umf4csol_ip
end interface

interface umf4cscal
module procedure umf4cscal,umf4cscal_ip
end interface
interface umf4zscal
module procedure umf4zscal,umf4zscal_ip,umf4cscal,umf4cscal_ip
end interface
interface umf4scal
module procedure umf4scal,umf4scal_ip,umf4zscal,umf4zscal_ip,umf4cscal,umf4cscal_ip
end interface

interface umf4cfnum
module procedure umf4cfnum,umf4cfnum_ip
end interface
interface umf4zfnum
module procedure umf4zfnum,umf4zfnum_ip
end interface
interface umf4fnum
module procedure umf4fnum,umf4fnum_ip
end interface

interface umf4cfsym
module procedure umf4cfsym,umf4cfsym_ip
end interface
interface umf4zfsym
module procedure umf4zfsym,umf4zfsym_ip
end interface
interface umf4fsym
module procedure umf4fsym,umf4fsym_ip
end interface

interface umf4csnum
module procedure umf4csnum,umf4csnum_ip
end interface
interface umf4zsnum
module procedure umf4zsnum,umf4zsnum_ip
end interface
interface umf4snum
module procedure umf4snum,umf4snum_ip
end interface

interface umf4cssym
module procedure umf4cssym,umf4cssym_ip
end interface
interface umf4zssym
module procedure umf4zssym,umf4zssym_ip
end interface
interface umf4ssym
module procedure umf4ssym,umf4ssym_ip
end interface

interface umf4clnum
module procedure umf4clnum,umf4clnum_ip
end interface
interface umf4zlnum
module procedure umf4zlnum,umf4zlnum_ip
end interface
interface umf4lnum
module procedure umf4lnum,umf4lnum_ip
end interface

interface umf4clsym
module procedure umf4clsym,umf4clsym_ip
end interface
interface umf4zlsym
module procedure umf4zlsym,umf4zlsym_ip
end interface
interface umf4lsym
module procedure umf4lsym,umf4lsym_ip
end interface

! ======================================================================
! Defined operator .umfpack.
! ======================================================================
type tCSC_di
integer,allocatable  :: Ap(:)
integer,allocatable  :: Ai(:)
real(r8),allocatable :: Ax(:)
end type

type tCSC_zi
integer,allocatable  :: Ap(:)
integer,allocatable  :: Ai(:)
real(r8),allocatable :: Ax(:)
real(r8),allocatable :: Az(:)
end type

type tCSC_ci
integer,allocatable     :: Ap(:)
integer,allocatable     :: Ai(:)
complex(r8),allocatable :: Ax(:)
end type

type tCSR_di
integer,allocatable  :: Ap(:)
integer,allocatable  :: Ai(:)
real(r8),allocatable :: Ax(:)
end type

type tCSR_zi
integer,allocatable  :: Ap(:)
integer,allocatable  :: Ai(:)
real(r8),allocatable :: Ax(:)
real(r8),allocatable :: Az(:)
end type

type tCSR_ci
integer,allocatable     :: Ap(:)
integer,allocatable     :: Ai(:)
complex(r8),allocatable :: Ax(:)
end type

type tVec_zi
real(r8),allocatable :: x(:)
real(r8),allocatable :: z(:)
end type

type pCSC_di
integer,pointer  :: Ap(:)
integer,pointer  :: Ai(:)
real(r8),pointer :: Ax(:)
end type

type pCSC_zi
integer,pointer  :: Ap(:)
integer,pointer  :: Ai(:)
real(r8),pointer :: Ax(:)
real(r8),pointer :: Az(:)
end type

type pCSC_ci
integer,pointer  :: Ap(:)
integer,pointer  :: Ai(:)
complex(r8),pointer :: Ax(:)
end type

type pCSR_di
integer,pointer  :: Ap(:)
integer,pointer  :: Ai(:)
real(r8),pointer :: Ax(:)
end type

type pCSR_zi
integer,pointer  :: Ap(:)
integer,pointer  :: Ai(:)
real(r8),pointer :: Ax(:)
real(r8),pointer :: Az(:)
end type

type pCSR_ci
integer,pointer  :: Ap(:)
integer,pointer  :: Ai(:)
complex(r8),pointer :: Ax(:)
end type

type pVec_zi
real(r8),pointer :: x(:)
real(r8),pointer :: z(:)
end type

! overloaded structure constructors (not implemented in g95)
! interface pCSC_di
! module procedure make_CSC_di
! end interface
! interface pCSC_zi
! module procedure make_CSC_zi
! end interface
! interface pCSC_ci
! module procedure make_CSC_ci
! end interface
interface pCSC
module procedure make_CSC_di,make_CSC_zi,make_CSC_ci
end interface
interface pCSR
module procedure make_CSR_di,make_CSR_zi,make_CSR_ci
end interface
! interface pVec_zi
! module procedure make_Vec_zi
! end interface
interface pVec
module procedure make_Vec_zi
end interface

interface operator(.umfpack.)
module procedure umfpack_di_operator_CSC,umfpack_zi_operator_CSC,umfpack_ci_operator_CSC, &
                 umfpack_di_operator_CSR,umfpack_zi_operator_CSR,umfpack_ci_operator_CSR, &
                 umfpack_di_operator_pCSC,umfpack_zi_operator_pCSC,umfpack_ci_operator_pCSC, &
                 umfpack_di_operator_pCSR,umfpack_zi_operator_pCSR,umfpack_ci_operator_pCSR
end interface

contains

! ======================================================================
! Full UMFPACK interface: Fortran wrappers
! ======================================================================

integer function umfpack_di_symbolic(n_row,n_col,Ap,Ai,Ax,Symbolic,Control,Info)
integer,intent(in) :: n_row,n_col
integer,target,intent(in) :: Ap(*),Ai(*)
real(r8),target,optional,intent(in) :: Ax(*)
type(c_ptr),optional :: Symbolic
real(r8),target,optional :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
integer :: c_n_row,c_n_col
type(c_ptr) :: c_Ap,c_Ai,c_Ax,c_Symbolic,c_Control,c_Info
c_n_row=n_row
c_n_col=n_col
c_Ap=c_loc(Ap)
c_Ai=c_loc(Ai)
if (present(Ax)) then ; c_Ax=c_loc(Ax) ; else ; c_Ax=c_null_ptr ; endif
if (present(Symbolic)) then ; c_Symbolic=Symbolic ; else ; c_Symbolic=iSymbolic ; endif
if (present(Control)) then ; c_Control=c_loc(Control) ; else ; c_Control=c_null_ptr ; endif
if (present(Info)) then ; c_Info=c_loc(Info) ; else ; c_Info=c_null_ptr ; endif
umfpack_di_symbolic=c_umfpack_di_symbolic(c_n_row,c_n_col,c_Ap,c_Ai,c_Ax,c_Symbolic,c_Control,c_Info)
if (present(Symbolic)) then ; Symbolic=c_Symbolic ; else ; iSymbolic=c_Symbolic ; endif
end function

subroutine s_umfpack_di_symbolic(n_row,n_col,Ap,Ai,Ax,Symbolic,Control,Info,status)
integer,intent(in) :: n_row,n_col
integer,intent(in) :: Ap(*),Ai(*)
real(r8),optional,intent(in) :: Ax(*)
type(c_ptr),optional :: Symbolic
real(r8),optional :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
integer,optional :: status
integer :: c_status
c_status=umfpack_di_symbolic(n_row,n_col,Ap,Ai,Ax,Symbolic,Control,Info)
if (present(status)) status=c_status
end subroutine

integer function umfpack_zi_symbolic(n_row,n_col,Ap,Ai,Ax,Az,Symbolic,Control,Info)
integer,intent(in) :: n_row,n_col
integer,target,intent(in) :: Ap(*),Ai(*)
real(r8),target,optional,intent(in) :: Ax(*)
real(r8),target,intent(in) :: Az(*)
type(c_ptr),optional :: Symbolic
real(r8),target,optional :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
integer :: c_n_row,c_n_col
type(c_ptr) :: c_Ap,c_Ai,c_Ax,c_Az,c_Symbolic,c_Control,c_Info
c_n_row=n_row
c_n_col=n_col
c_Ap=c_loc(Ap)
c_Ai=c_loc(Ai)
if (present(Ax)) then ; c_Ax=c_loc(Ax) ; else ; c_Ax=c_null_ptr ; endif
c_Az=c_loc(Az)
if (present(Symbolic)) then ; c_Symbolic=Symbolic ; else ; c_Symbolic=iSymbolic ; endif
if (present(Control)) then ; c_Control=c_loc(Control) ; else ; c_Control=c_null_ptr ; endif
if (present(Info)) then ; c_Info=c_loc(Info) ; else ; c_Info=c_null_ptr ; endif
umfpack_zi_symbolic=c_umfpack_zi_symbolic(c_n_row,c_n_col,c_Ap,c_Ai,c_Ax,c_Az,c_Symbolic,c_Control,c_Info)
if (present(Symbolic)) then ; Symbolic=c_Symbolic ; else ; iSymbolic=c_Symbolic ; endif
end function

subroutine s_umfpack_zi_symbolic(n_row,n_col,Ap,Ai,Ax,Az,Symbolic,Control,Info,status)
integer,intent(in) :: n_row,n_col
integer,intent(in) :: Ap(*),Ai(*)
real(r8),optional,intent(in) :: Ax(*)
real(r8),intent(in) :: Az(*)
type(c_ptr),optional :: Symbolic
real(r8),optional :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
integer,optional :: status
integer :: c_status
c_status=umfpack_zi_symbolic(n_row,n_col,Ap,Ai,Ax,Az,Symbolic,Control,Info)
if (present(status)) status=c_status
end subroutine

integer function umfpack_ci_symbolic(n_row,n_col,Ap,Ai,Ax,Symbolic,Control,Info)
integer,intent(in) :: n_row,n_col
integer,target,intent(in) :: Ap(*),Ai(*)
complex(r8),target,intent(in) :: Ax(*)
type(c_ptr),optional :: Symbolic
real(r8),target,optional :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
integer :: c_n_row,c_n_col
type(c_ptr) :: c_Ap,c_Ai,c_Ax,c_Az,c_Symbolic,c_Control,c_Info
c_n_row=n_row
c_n_col=n_col
c_Ap=c_loc(Ap)
c_Ai=c_loc(Ai)
c_Ax=c_loc(Ax)
c_Az=c_null_ptr
if (present(Symbolic)) then ; c_Symbolic=Symbolic ; else ; c_Symbolic=iSymbolic ; endif
if (present(Control)) then ; c_Control=c_loc(Control) ; else ; c_Control=c_null_ptr ; endif
if (present(Info)) then ; c_Info=c_loc(Info) ; else ; c_Info=c_null_ptr ; endif
umfpack_ci_symbolic=c_umfpack_zi_symbolic(c_n_row,c_n_col,c_Ap,c_Ai,c_Ax,c_Az,c_Symbolic,c_Control,c_Info)
if (present(Symbolic)) then ; Symbolic=c_Symbolic ; else ; iSymbolic=c_Symbolic ; endif
end function

subroutine s_umfpack_ci_symbolic(n_row,n_col,Ap,Ai,Ax,Symbolic,Control,Info,status)
integer,intent(in) :: n_row,n_col
integer,intent(in) :: Ap(*),Ai(*)
complex(r8),intent(in) :: Ax(*)
type(c_ptr),optional :: Symbolic
real(r8),optional :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
integer,optional :: status
integer :: c_status
c_status=umfpack_ci_symbolic(n_row,n_col,Ap,Ai,Ax,Symbolic,Control,Info)
if (present(status)) status=c_status
end subroutine

integer function umfpack_di_numeric(Ap,Ai,Ax,Symbolic,Numeric,Control,Info)
integer,target,intent(in) :: Ap(*),Ai(*)
real(r8),target,intent(in) :: Ax(*)
type(c_ptr),optional,intent(in) :: Symbolic
type(c_ptr),optional :: Numeric
real(r8),target,optional :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
type(c_ptr) :: c_Ap,c_Ai,c_Ax,c_Symbolic,c_Numeric,c_Control,c_Info
c_Ap=c_loc(Ap)
c_Ai=c_loc(Ai)
c_Ax=c_loc(Ax)
if (present(Symbolic)) then ; c_Symbolic=Symbolic ; else ; c_Symbolic=iSymbolic ; endif
if (present(Numeric)) then ; c_Numeric=Numeric ; else ; c_Numeric=iNumeric ; endif
if (present(Control)) then ; c_Control=c_loc(Control) ; else ; c_Control=c_null_ptr ; endif
if (present(Info)) then ; c_Info=c_loc(Info) ; else ; c_Info=c_null_ptr ; endif
umfpack_di_numeric=c_umfpack_di_numeric(c_Ap,c_Ai,c_Ax,c_Symbolic,c_Numeric,c_Control,c_Info)
if (present(Numeric)) then ; Numeric=c_Numeric ; else ; iNumeric=c_Numeric ; endif
end function

subroutine s_umfpack_di_numeric(Ap,Ai,Ax,Symbolic,Numeric,Control,Info,status)
integer,intent(in) :: Ap(*),Ai(*)
real(r8),intent(in) :: Ax(*)
type(c_ptr),optional,intent(in) :: Symbolic
type(c_ptr),optional :: Numeric
real(r8),optional :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
integer,optional :: status
integer :: c_status
c_status=umfpack_di_numeric(Ap,Ai,Ax,Symbolic,Numeric,Control,Info)
if (present(status)) status=c_status
end subroutine

integer function umfpack_zi_numeric(Ap,Ai,Ax,Az,Symbolic,Numeric,Control,Info)
integer,target,intent(in) :: Ap(*),Ai(*)
real(r8),target,intent(in) :: Ax(*),Az(*)
type(c_ptr),optional,intent(in) :: Symbolic
type(c_ptr),optional :: Numeric
real(r8),target,optional :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
type(c_ptr) :: c_Ap,c_Ai,c_Ax,c_Az,c_Symbolic,c_Numeric,c_Control,c_Info
c_Ap=c_loc(Ap)
c_Ai=c_loc(Ai)
c_Ax=c_loc(Ax)
c_Az=c_loc(Az)
if (present(Symbolic)) then ; c_Symbolic=Symbolic ; else ; c_Symbolic=iSymbolic ; endif
if (present(Numeric)) then ; c_Numeric=Numeric ; else ; c_Numeric=iNumeric ; endif
if (present(Control)) then ; c_Control=c_loc(Control) ; else ; c_Control=c_null_ptr ; endif
if (present(Info)) then ; c_Info=c_loc(Info) ; else ; c_Info=c_null_ptr ; endif
umfpack_zi_numeric=c_umfpack_zi_numeric(c_Ap,c_Ai,c_Ax,c_Az,c_Symbolic,c_Numeric,c_Control,c_Info)
if (present(Numeric)) then ; Numeric=c_Numeric ; else ; iNumeric=c_Numeric ; endif
end function

subroutine s_umfpack_zi_numeric(Ap,Ai,Ax,Az,Symbolic,Numeric,Control,Info,status)
integer,intent(in) :: Ap(*),Ai(*)
real(r8),intent(in) :: Ax(*),Az(*)
type(c_ptr),optional,intent(in) :: Symbolic
type(c_ptr),optional :: Numeric
real(r8),optional :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
integer,optional :: status
integer :: c_status
c_status=umfpack_zi_numeric(Ap,Ai,Ax,Az,Symbolic,Numeric,Control,Info)
if (present(status)) status=c_status
end subroutine

integer function umfpack_ci_numeric(Ap,Ai,Ax,Symbolic,Numeric,Control,Info)
integer,target,intent(in) :: Ap(*),Ai(*)
complex(r8),target,intent(in) :: Ax(*)
type(c_ptr),optional,intent(in) :: Symbolic
type(c_ptr),optional :: Numeric
real(r8),target,optional :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
type(c_ptr) :: c_Ap,c_Ai,c_Ax,c_Az,c_Symbolic,c_Numeric,c_Control,c_Info
c_Ap=c_loc(Ap)
c_Ai=c_loc(Ai)
c_Ax=c_loc(Ax)
c_Az=c_null_ptr
if (present(Symbolic)) then ; c_Symbolic=Symbolic ; else ; c_Symbolic=iSymbolic ; endif
if (present(Numeric)) then ; c_Numeric=Numeric ; else ; c_Numeric=iNumeric ; endif
if (present(Control)) then ; c_Control=c_loc(Control) ; else ; c_Control=c_null_ptr ; endif
if (present(Info)) then ; c_Info=c_loc(Info) ; else ; c_Info=c_null_ptr ; endif
umfpack_ci_numeric=c_umfpack_zi_numeric(c_Ap,c_Ai,c_Ax,c_Az,c_Symbolic,c_Numeric,c_Control,c_Info)
if (present(Numeric)) then ; Numeric=c_Numeric ; else ; iNumeric=c_Numeric ; endif
end function

subroutine s_umfpack_ci_numeric(Ap,Ai,Ax,Symbolic,Numeric,Control,Info,status)
integer,intent(in) :: Ap(*),Ai(*)
complex(r8),intent(in) :: Ax(*)
type(c_ptr),optional,intent(in) :: Symbolic
type(c_ptr),optional :: Numeric
real(r8),optional :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
integer,optional :: status
integer :: c_status
c_status=umfpack_ci_numeric(Ap,Ai,Ax,Symbolic,Numeric,Control,Info)
if (present(status)) status=c_status
end subroutine

integer function umfpack_di_solve(sys,Ap,Ai,Ax,X,B,Numeric,Control,Info)
integer,optional,intent(in) :: sys
integer,target,optional,intent(in) :: Ap(*),Ai(*)
real(r8),target,optional,intent(in) :: Ax(*)
real(r8),target :: X(*)
real(r8),target,intent(in) :: B(*)
type(c_ptr),optional,intent(in) :: Numeric
real(r8),target,optional :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
integer :: c_sys
type(c_ptr) :: c_Ap,c_Ai,c_Ax,c_X,c_B,c_Numeric,c_Control,c_Info
if (present(sys)) then ; c_sys=sys ; else ; c_sys=UMFPACK_A ; endif
if (present(Ap)) then ; c_Ap=c_loc(Ap) ; else ; c_Ap=c_null_ptr ; endif
if (present(Ai)) then ; c_Ai=c_loc(Ai) ; else ; c_Ai=c_null_ptr ; endif
if (present(Ax)) then ; c_Ax=c_loc(Ax) ; else ; c_Ax=c_null_ptr ; endif
c_X=c_loc(X)
c_B=c_loc(B)
if (present(Numeric)) then ; c_Numeric=Numeric ; else ; c_Numeric=iNumeric ; endif
if (present(Control)) then ; c_Control=c_loc(Control) ; else ; c_Control=c_null_ptr ; endif
if (present(Info)) then ; c_Info=c_loc(Info) ; else ; c_Info=c_null_ptr ; endif
umfpack_di_solve=c_umfpack_di_solve(c_sys,c_Ap,c_Ai,c_Ax,c_X,c_B,c_Numeric,c_Control,c_Info)
end function

subroutine s_umfpack_di_solve(sys,Ap,Ai,Ax,X,B,Numeric,Control,Info,status)
integer,optional,intent(in) :: sys
integer,optional,intent(in) :: Ap(*),Ai(*)
real(r8),optional,intent(in) :: Ax(*)
real(r8),intent(out) :: X(*)
real(r8),intent(in) :: B(*)
type(c_ptr),optional,intent(in) :: Numeric
real(r8),optional :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
integer,optional :: status
integer :: c_status
c_status=umfpack_di_solve(sys,Ap,Ai,Ax,X,B,Numeric,Control,Info)
if (present(status)) status=c_status
end subroutine

integer function umfpack_zi_solve(sys,Ap,Ai,Ax,Az,Xx,Xz,Bx,Bz,Numeric,Control,Info)
integer,optional,intent(in) :: sys
integer,target,optional,intent(in) :: Ap(*),Ai(*)
real(r8),target,optional,intent(in) :: Ax(*),Az(*)
real(r8),target :: Xx(*),Xz(*)
real(r8),target,intent(in) :: Bx(*),Bz(*)
type(c_ptr),optional,intent(in) :: Numeric
real(r8),target,optional :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
integer :: c_sys
type(c_ptr) :: c_Ap,c_Ai,c_Ax,c_Az,c_Xx,c_Xz,c_Bx,c_Bz,c_Numeric,c_Control,c_Info
if (present(sys)) then ; c_sys=sys ; else ; c_sys=UMFPACK_A ; endif
if (present(Ap)) then ; c_Ap=c_loc(Ap) ; else ; c_Ap=c_null_ptr ; endif
if (present(Ai)) then ; c_Ai=c_loc(Ai) ; else ; c_Ai=c_null_ptr ; endif
if (present(Ax)) then ; c_Ax=c_loc(Ax) ; else ; c_Ax=c_null_ptr ; endif
if (present(Az)) then ; c_Az=c_loc(Az) ; else ; c_Az=c_null_ptr ; endif
c_Xx=c_loc(Xx)
c_Xz=c_loc(Xz)
c_Bx=c_loc(Bx)
c_Bz=c_loc(Bz)
if (present(Numeric)) then ; c_Numeric=Numeric ; else ; c_Numeric=iNumeric ; endif
if (present(Control)) then ; c_Control=c_loc(Control) ; else ; c_Control=c_null_ptr ; endif
if (present(Info)) then ; c_Info=c_loc(Info) ; else ; c_Info=c_null_ptr ; endif
umfpack_zi_solve=c_umfpack_zi_solve(c_sys,c_Ap,c_Ai,c_Ax,c_Az,c_Xx,c_Xz,c_Bx,c_Bz,c_Numeric,c_Control,c_Info)
end function

subroutine s_umfpack_zi_solve(sys,Ap,Ai,Ax,Az,Xx,Xz,Bx,Bz,Numeric,Control,Info,status)
integer,optional,intent(in) :: sys
integer,optional,intent(in) :: Ap(*),Ai(*)
real(r8),optional,intent(in) :: Ax(*),Az(*)
real(r8),intent(out) :: Xx(*),Xz(*)
real(r8),intent(in) :: Bx(*),Bz(*)
type(c_ptr),optional,intent(in) :: Numeric
real(r8),optional :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
integer,optional :: status
integer :: c_status
c_status=umfpack_zi_solve(sys,Ap,Ai,Ax,Az,Xx,Xz,Bx,Bz,Numeric,Control,Info)
if (present(status)) status=c_status
end subroutine

integer function umfpack_ci_solve(sys,Ap,Ai,Ax,Xx,Bx,Numeric,Control,Info)
integer,optional,intent(in) :: sys
integer,target,optional,intent(in) :: Ap(*),Ai(*)
complex(r8),target,optional,intent(in) :: Ax(*)
complex(r8),target :: Xx(*)
complex(r8),target,intent(in) :: Bx(*)
type(c_ptr),optional,intent(in) :: Numeric
real(r8),target,optional :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
integer :: c_sys
type(c_ptr) :: c_Ap,c_Ai,c_Ax,c_Az,c_Xx,c_Xz,c_Bx,c_Bz,c_Numeric,c_Control,c_Info
if (present(sys)) then ; c_sys=sys ; else ; c_sys=UMFPACK_A ; endif
if (present(Ap)) then ; c_Ap=c_loc(Ap) ; else ; c_Ap=c_null_ptr ; endif
if (present(Ai)) then ; c_Ai=c_loc(Ai) ; else ; c_Ai=c_null_ptr ; endif
if (present(Ax)) then ; c_Ax=c_loc(Ax) ; else ; c_Ax=c_null_ptr ; endif
c_Az=c_null_ptr
c_Xx=c_loc(Xx)
c_Xz=c_null_ptr
c_Bx=c_loc(Bx)
c_Bz=c_null_ptr
if (present(Numeric)) then ; c_Numeric=Numeric ; else ; c_Numeric=iNumeric ; endif
if (present(Control)) then ; c_Control=c_loc(Control) ; else ; c_Control=c_null_ptr ; endif
if (present(Info)) then ; c_Info=c_loc(Info) ; else ; c_Info=c_null_ptr ; endif
umfpack_ci_solve=c_umfpack_zi_solve(c_sys,c_Ap,c_Ai,c_Ax,c_Az,c_Xx,c_Xz,c_Bx,c_Bz,c_Numeric,c_Control,c_Info)
end function

subroutine s_umfpack_ci_solve(sys,Ap,Ai,Ax,Xx,Bx,Numeric,Control,Info,status)
integer,optional,intent(in) :: sys
integer,optional,intent(in) :: Ap(*),Ai(*)
complex(r8),optional,intent(in) :: Ax(*)
complex(r8),intent(out) :: Xx(*)
complex(r8),intent(in) :: Bx(*)
type(c_ptr),optional,intent(in) :: Numeric
real(r8),optional :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
integer,optional :: status
integer :: c_status
c_status=umfpack_ci_solve(sys,Ap,Ai,Ax,Xx,Bx,Numeric,Control,Info)
if (present(status)) status=c_status
end subroutine

subroutine umfpack_di_free_symbolic(Symbolic)
type(c_ptr),optional :: Symbolic
type(c_ptr) :: c_Symbolic
if (present(Symbolic)) then ; c_Symbolic=Symbolic ; else ; c_Symbolic=iSymbolic ; endif
call c_umfpack_di_free_symbolic(c_Symbolic)
if (present(Symbolic)) then ; Symbolic=c_Symbolic ; else ; iSymbolic=c_Symbolic ; endif
end subroutine

subroutine umfpack_zi_free_symbolic(Symbolic)
type(c_ptr),optional :: Symbolic
type(c_ptr) :: c_Symbolic
if (present(Symbolic)) then ; c_Symbolic=Symbolic ; else ; c_Symbolic=iSymbolic ; endif
call c_umfpack_zi_free_symbolic(c_Symbolic)
if (present(Symbolic)) then ; Symbolic=c_Symbolic ; else ; iSymbolic=c_Symbolic ; endif
end subroutine

subroutine umfpack_free_symbolic(Symbolic,version)
type(c_ptr),optional :: Symbolic
character(*),optional :: version
type(c_ptr) :: c_Symbolic
character(2) :: c_version
if (present(Symbolic)) then ; c_Symbolic=Symbolic ; else ; c_Symbolic=iSymbolic ; endif
if (present(version)) then ; c_version=version ; else ; c_version="di" ; endif
select case (c_version)
case ("zi") ; call umfpack_zi_free_symbolic(c_Symbolic)
case ("ci") ; call umfpack_zi_free_symbolic(c_Symbolic)
case default ; call umfpack_di_free_symbolic(c_Symbolic)
end select
end subroutine

subroutine umfpack_di_free_numeric(Numeric)
type(c_ptr),optional :: Numeric
type(c_ptr) :: c_Numeric
if (present(Numeric)) then ; c_Numeric=Numeric ; else ; c_Numeric=iNumeric ; endif
call c_umfpack_di_free_numeric(c_Numeric)
if (present(Numeric)) then ; Numeric=c_Numeric ; else ; iNumeric=c_Numeric ; endif
end subroutine

subroutine umfpack_zi_free_numeric(Numeric)
type(c_ptr),optional :: Numeric
type(c_ptr) :: c_Numeric
if (present(Numeric)) then ; c_Numeric=Numeric ; else ; c_Numeric=iNumeric ; endif
call c_umfpack_zi_free_numeric(c_Numeric)
if (present(Numeric)) then ; Numeric=c_Numeric ; else ; iNumeric=c_Numeric ; endif
end subroutine

subroutine umfpack_free_numeric(Numeric,version)
type(c_ptr),optional :: Numeric
character(*),optional :: version
type(c_ptr) :: c_Numeric
character(2) :: c_version
if (present(Numeric)) then ; c_Numeric=Numeric ; else ; c_Numeric=iNumeric ; endif
if (present(version)) then ; c_version=version ; else ; c_version="di" ; endif
select case (c_version)
case ("zi") ; call umfpack_zi_free_numeric(c_Numeric)
case ("ci") ; call umfpack_zi_free_numeric(c_Numeric)
case default ; call umfpack_di_free_numeric(c_Numeric)
end select
end subroutine

integer function umfpack_di_scale(X,B,Numeric)
real(r8),target :: X(*)
real(r8),target,intent(in) :: B(*)
type(c_ptr),optional,intent(in) :: Numeric
type(c_ptr) :: c_X,c_B,c_Numeric
c_X=c_loc(X)
c_B=c_loc(B)
if (present(Numeric)) then ; c_Numeric=Numeric ; else ; c_Numeric=iNumeric ; endif
umfpack_di_scale=c_umfpack_di_scale(c_X,c_B,c_Numeric)
end function

subroutine s_umfpack_di_scale(X,B,Numeric,status)
real(r8),intent(out) :: X(*)
real(r8),intent(in) :: B(*)
type(c_ptr),optional,intent(in) :: Numeric
integer,optional :: status
integer :: c_status
c_status=umfpack_di_scale(X,B,Numeric)
if (present(status)) status=c_status
end subroutine

integer function umfpack_zi_scale(Xx,Xz,Bx,Bz,Numeric)
real(r8),target :: Xx(*),Xz(*)
real(r8),target,intent(in) :: Bx(*),Bz(*)
type(c_ptr),optional,intent(in) :: Numeric
type(c_ptr) :: c_Xx,c_Xz,c_Bx,c_Bz,c_Numeric
c_Xx=c_loc(Xx)
c_Xz=c_loc(Xz)
c_Bx=c_loc(Bx)
c_Bz=c_loc(Bz)
if (present(Numeric)) then ; c_Numeric=Numeric ; else ; c_Numeric=iNumeric ; endif
umfpack_zi_scale=c_umfpack_zi_scale(c_Xx,c_Xz,c_Bx,c_Bz,c_Numeric)
end function

subroutine s_umfpack_zi_scale(Xx,Xz,Bx,Bz,Numeric,status)
real(r8),intent(out) :: Xx(*),Xz(*)
real(r8),intent(in) :: Bx(*),Bz(*)
type(c_ptr),optional,intent(in) :: Numeric
integer,optional :: status
integer :: c_status
c_status=umfpack_zi_scale(Xx,Xz,Bx,Bz,Numeric)
if (present(status)) status=c_status
end subroutine

integer function umfpack_ci_scale(Xx,Bx,Numeric)
complex(r8),target :: Xx(*)
complex(r8),target,intent(in) :: Bx(*)
type(c_ptr),optional,intent(in) :: Numeric
type(c_ptr) :: c_Xx,c_Xz,c_Bx,c_Bz,c_Numeric
c_Xx=c_loc(Xx)
c_Xz=c_null_ptr
c_Bx=c_loc(Bx)
c_Bz=c_null_ptr
if (present(Numeric)) then ; c_Numeric=Numeric ; else ; c_Numeric=iNumeric ; endif
umfpack_ci_scale=c_umfpack_zi_scale(c_Xx,c_Xz,c_Bx,c_Bz,c_Numeric)
end function

subroutine s_umfpack_ci_scale(Xx,Bx,Numeric,status)
complex(r8),intent(out) :: Xx(*)
complex(r8),intent(in) :: Bx(*)
type(c_ptr),optional,intent(in) :: Numeric
integer,optional :: status
integer :: c_status
c_status=umfpack_ci_scale(Xx,Bx,Numeric)
if (present(status)) status=c_status
end subroutine

subroutine umfpack_di_defaults(Control)
real(r8),target,optional :: Control(0:UMFPACK_CONTROL-1)
type(c_ptr) :: c_Control
if (present(Control)) then ; c_Control=c_loc(Control) ; else ; c_Control=c_null_ptr ; endif
call c_umfpack_di_defaults(c_Control)
end subroutine

subroutine umfpack_zi_defaults(Control)
real(r8),target,optional :: Control(0:UMFPACK_CONTROL-1)
type(c_ptr) :: c_Control
if (present(Control)) then ; c_Control=c_loc(Control) ; else ; c_Control=c_null_ptr ; endif
call c_umfpack_zi_defaults(c_Control)
end subroutine

subroutine umfpack_defaults(Control,version)
real(r8),optional :: Control(0:UMFPACK_CONTROL-1)
character(*),optional :: version
character(2) :: c_version
if (present(version)) then ; c_version=version ; else ; c_version="di" ; endif
select case (c_version)
case ("zi") ; call umfpack_zi_defaults(Control)
case ("ci") ; call umfpack_zi_defaults(Control)
case default ; call umfpack_di_defaults(Control)
end select
end subroutine

subroutine umfpack_di_report_control(Control)
use iso_fortran_env,only : output_unit
real(r8),target,optional :: Control(0:UMFPACK_CONTROL-1)
type(c_ptr) :: c_Control
logical :: isOpened
if (present(Control)) then ; c_Control=c_loc(Control) ; else ; c_Control=c_null_ptr ; endif
inquire (output_unit,opened=isOpened)
if (isOpened) flush (output_unit)
call c_umfpack_di_report_control(c_Control)
if (isOpened) flush (output_unit)
end subroutine

subroutine umfpack_zi_report_control(Control)
use iso_fortran_env,only : output_unit
real(r8),target,optional :: Control(0:UMFPACK_CONTROL-1)
type(c_ptr) :: c_Control
logical :: isOpened
if (present(Control)) then ; c_Control=c_loc(Control) ; else ; c_Control=c_null_ptr ; endif
inquire (output_unit,opened=isOpened)
if (isOpened) flush (output_unit)
call c_umfpack_zi_report_control(c_Control)
if (isOpened) flush (output_unit)
end subroutine

subroutine umfpack_report_control(Control,version)
real(r8),optional :: Control(0:UMFPACK_CONTROL-1)
character(*),optional :: version
character(2) :: c_version
if (present(version)) then ; c_version=version ; else ; c_version="di" ; endif
select case (c_version)
case ("zi") ; call umfpack_zi_report_control(Control)
case ("ci") ; call umfpack_zi_report_control(Control)
case default ; call umfpack_di_report_control(Control)
end select
end subroutine

subroutine umfpack_di_report_info(Control,Info)
use iso_fortran_env,only : output_unit
real(r8),target,optional :: Control(0:UMFPACK_CONTROL-1)
real(r8),target :: Info(0:UMFPACK_INFO-1)
type(c_ptr) :: c_Control,c_Info
logical :: isOpened
if (present(Control)) then ; c_Control=c_loc(Control) ; else ; c_Control=c_null_ptr ; endif
c_Info=c_loc(Info)
inquire (output_unit,opened=isOpened)
if (isOpened) flush (output_unit)
call c_umfpack_di_report_info(c_Control,c_Info)
if (isOpened) flush (output_unit)
end subroutine

subroutine umfpack_zi_report_info(Control,Info)
use iso_fortran_env,only : output_unit
real(r8),target,optional :: Control(0:UMFPACK_CONTROL-1)
real(r8),target :: Info(0:UMFPACK_INFO-1)
type(c_ptr) :: c_Control,c_Info
logical :: isOpened
if (present(Control)) then ; c_Control=c_loc(Control) ; else ; c_Control=c_null_ptr ; endif
c_Info=c_loc(Info)
inquire (output_unit,opened=isOpened)
if (isOpened) flush (output_unit)
call c_umfpack_zi_report_info(c_Control,c_Info)
if (isOpened) flush (output_unit)
end subroutine

subroutine umfpack_report_info(Control,Info,version)
real(r8),optional :: Control(0:UMFPACK_CONTROL-1)
real(r8) :: Info(0:UMFPACK_INFO-1)
character(*),optional :: version
character(2) :: c_version
if (present(version)) then ; c_version=version ; else ; c_version="di" ; endif
select case (c_version)
case ("zi") ; call umfpack_zi_report_info(Control,Info)
case ("ci") ; call umfpack_zi_report_info(Control,Info)
case default ; call umfpack_di_report_info(Control,Info)
end select
end subroutine

integer function umfpack_di_save_numeric(Numeric,filename)
type(c_ptr),intent(in) :: Numeric
character(*),intent(in) :: filename
umfpack_di_save_numeric=c_umfpack_di_save_numeric(Numeric,trim(filename)//c_null_char)
end function

subroutine s_umfpack_di_save_numeric(Numeric,filename,status)
type(c_ptr),intent(in) :: Numeric
character(*),intent(in) :: filename
integer,optional :: status
integer :: c_status
c_status=umfpack_di_save_numeric(Numeric,filename)
if (present(status)) status=c_status
end subroutine

integer function umfpack_zi_save_numeric(Numeric,filename)
type(c_ptr),intent(in) :: Numeric
character(*),intent(in) :: filename
umfpack_zi_save_numeric=c_umfpack_zi_save_numeric(Numeric,trim(filename)//c_null_char)
end function

subroutine s_umfpack_zi_save_numeric(Numeric,filename,status)
type(c_ptr),intent(in) :: Numeric
character(*),intent(in) :: filename
integer,optional :: status
integer :: c_status
c_status=umfpack_zi_save_numeric(Numeric,filename)
if (present(status)) status=c_status
end subroutine

integer function umfpack_save_numeric(Numeric,filename,version)
type(c_ptr),intent(in) :: Numeric
character(*),intent(in) :: filename
character(*),optional :: version
character(2) :: c_version
if (present(version)) then ; c_version=version ; else ; c_version="di" ; endif
select case (c_version)
case ("zi") ; umfpack_save_numeric=umfpack_zi_save_numeric(Numeric,filename)
case ("ci") ; umfpack_save_numeric=umfpack_zi_save_numeric(Numeric,filename)
case default ; umfpack_save_numeric=umfpack_di_save_numeric(Numeric,filename)
end select
end function

subroutine s_umfpack_save_numeric(Numeric,filename,status,version)
type(c_ptr),intent(in) :: Numeric
character(*),intent(in) :: filename
integer,optional :: status
character(*),optional :: version
integer :: c_status
c_status=umfpack_save_numeric(Numeric,filename,version)
if (present(status)) status=c_status
end subroutine

integer function umfpack_di_save_symbolic(Symbolic,filename)
type(c_ptr),intent(in) :: Symbolic
character(*),intent(in) :: filename
umfpack_di_save_symbolic=c_umfpack_di_save_symbolic(Symbolic,trim(filename)//c_null_char)
end function

subroutine s_umfpack_di_save_symbolic(Symbolic,filename,status)
type(c_ptr),intent(in) :: Symbolic
character(*),intent(in) :: filename
integer,optional :: status
integer :: c_status
c_status=umfpack_di_save_symbolic(Symbolic,filename)
if (present(status)) status=c_status
end subroutine

integer function umfpack_zi_save_symbolic(Symbolic,filename)
type(c_ptr),intent(in) :: Symbolic
character(*),intent(in) :: filename
umfpack_zi_save_symbolic=c_umfpack_zi_save_symbolic(Symbolic,trim(filename)//c_null_char)
end function

subroutine s_umfpack_zi_save_symbolic(Symbolic,filename,status)
type(c_ptr),intent(in) :: Symbolic
character(*),intent(in) :: filename
integer,optional :: status
integer :: c_status
c_status=umfpack_zi_save_symbolic(Symbolic,filename)
if (present(status)) status=c_status
end subroutine

integer function umfpack_save_symbolic(Symbolic,filename,version)
type(c_ptr),intent(in) :: Symbolic
character(*),intent(in) :: filename
character(*),optional :: version
character(2) :: c_version
if (present(version)) then ; c_version=version ; else ; c_version="di" ; endif
select case (c_version)
case ("zi") ; umfpack_save_symbolic=umfpack_zi_save_symbolic(Symbolic,filename)
case ("ci") ; umfpack_save_symbolic=umfpack_zi_save_symbolic(Symbolic,filename)
case default ; umfpack_save_symbolic=umfpack_di_save_symbolic(Symbolic,filename)
end select
end function

subroutine s_umfpack_save_symbolic(Symbolic,filename,status,version)
type(c_ptr),intent(in) :: Symbolic
character(*),intent(in) :: filename
integer,optional :: status
character(*),optional :: version
integer :: c_status
c_status=umfpack_save_symbolic(Symbolic,filename,version)
if (present(status)) status=c_status
end subroutine

integer function umfpack_di_load_numeric(Numeric,filename)
type(c_ptr),intent(in) :: Numeric
character(*),intent(in) :: filename
umfpack_di_load_numeric=c_umfpack_di_load_numeric(Numeric,trim(filename)//c_null_char)
end function

subroutine s_umfpack_di_load_numeric(Numeric,filename,status)
type(c_ptr),intent(in) :: Numeric
character(*),intent(in) :: filename
integer,optional :: status
integer :: c_status
c_status=umfpack_di_load_numeric(Numeric,filename)
if (present(status)) status=c_status
end subroutine

integer function umfpack_zi_load_numeric(Numeric,filename)
type(c_ptr),intent(in) :: Numeric
character(*),intent(in) :: filename
umfpack_zi_load_numeric=c_umfpack_zi_load_numeric(Numeric,trim(filename)//c_null_char)
end function

subroutine s_umfpack_zi_load_numeric(Numeric,filename,status)
type(c_ptr),intent(in) :: Numeric
character(*),intent(in) :: filename
integer,optional :: status
integer :: c_status
c_status=umfpack_zi_load_numeric(Numeric,filename)
if (present(status)) status=c_status
end subroutine

integer function umfpack_load_numeric(Numeric,filename,version)
type(c_ptr),intent(in) :: Numeric
character(*),intent(in) :: filename
character(*),optional :: version
character(2) :: c_version
if (present(version)) then ; c_version=version ; else ; c_version="di" ; endif
select case (c_version)
case ("zi") ; umfpack_load_numeric=umfpack_zi_load_numeric(Numeric,filename)
case ("ci") ; umfpack_load_numeric=umfpack_zi_load_numeric(Numeric,filename)
case default ; umfpack_load_numeric=umfpack_di_load_numeric(Numeric,filename)
end select
end function

subroutine s_umfpack_load_numeric(Numeric,filename,status,version)
type(c_ptr),intent(in) :: Numeric
character(*),intent(in) :: filename
integer,optional :: status
character(*),optional :: version
integer :: c_status
c_status=umfpack_load_numeric(Numeric,filename,version)
if (present(status)) status=c_status
end subroutine

integer function umfpack_di_load_symbolic(Symbolic,filename)
type(c_ptr),intent(in) :: Symbolic
character(*),intent(in) :: filename
umfpack_di_load_symbolic=c_umfpack_di_load_symbolic(Symbolic,trim(filename)//c_null_char)
end function

subroutine s_umfpack_di_load_symbolic(Symbolic,filename,status)
type(c_ptr),intent(in) :: Symbolic
character(*),intent(in) :: filename
integer,optional :: status
integer :: c_status
c_status=umfpack_di_load_symbolic(Symbolic,filename)
if (present(status)) status=c_status
end subroutine

integer function umfpack_zi_load_symbolic(Symbolic,filename)
type(c_ptr),intent(in) :: Symbolic
character(*),intent(in) :: filename
umfpack_zi_load_symbolic=c_umfpack_zi_load_symbolic(Symbolic,trim(filename)//c_null_char)
end function

subroutine s_umfpack_zi_load_symbolic(Symbolic,filename,status)
type(c_ptr),intent(in) :: Symbolic
character(*),intent(in) :: filename
integer,optional :: status
integer :: c_status
c_status=umfpack_zi_load_symbolic(Symbolic,filename)
if (present(status)) status=c_status
end subroutine

integer function umfpack_load_symbolic(Symbolic,filename,version)
type(c_ptr),intent(in) :: Symbolic
character(*),intent(in) :: filename
character(*),optional :: version
character(2) :: c_version
if (present(version)) then ; c_version=version ; else ; c_version="di" ; endif
select case (c_version)
case ("zi") ; umfpack_load_symbolic=umfpack_zi_load_symbolic(Symbolic,filename)
case ("ci") ; umfpack_load_symbolic=umfpack_zi_load_symbolic(Symbolic,filename)
case default ; umfpack_load_symbolic=umfpack_di_load_symbolic(Symbolic,filename)
end select
end function

subroutine s_umfpack_load_symbolic(Symbolic,filename,status,version)
type(c_ptr),intent(in) :: Symbolic
character(*),intent(in) :: filename
integer,optional :: status
character(*),optional :: version
integer :: c_status
c_status=umfpack_load_symbolic(Symbolic,filename,version)
if (present(status)) status=c_status
end subroutine

! ======================================================================
! Basic UMFPACK interface: Fortran wrappers
! ======================================================================

subroutine umf4def(Control)
real(r8) :: Control(0:UMFPACK_CONTROL-1)
call umfpack_defaults(Control)
end subroutine

subroutine umf4zdef(Control)
real(r8) :: Control(0:UMFPACK_CONTROL-1)
call umfpack_defaults(Control,"zi")
end subroutine

subroutine umf4cdef(Control)
real(r8) :: Control(0:UMFPACK_CONTROL-1)
call umfpack_defaults(Control,"zi")
end subroutine

subroutine umf4pcon(Control)
real(r8) :: Control(0:UMFPACK_CONTROL-1)
call umfpack_report_control(Control)
end subroutine

subroutine umf4zpcon(Control)
real(r8) :: Control(0:UMFPACK_CONTROL-1)
call umfpack_report_control(Control,"zi")
end subroutine

subroutine umf4cpcon(Control)
real(r8) :: Control(0:UMFPACK_CONTROL-1)
call umfpack_report_control(Control,"zi")
end subroutine

subroutine umf4sym(m,n,Ap,Ai,Ax,Symbolic,Control,Info)
integer,intent(in) :: m,n
integer,intent(in) :: Ap(*),Ai(*)
real(r8),intent(in) :: Ax(*)
type(c_ptr) :: Symbolic
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
call s_umfpack_symbolic(m,n,Ap,Ai,Ax,Symbolic,Control,Info)
end subroutine

subroutine umf4zsym(m,n,Ap,Ai,Ax,Az,Symbolic,Control,Info)
integer,intent(in) :: m,n
integer,intent(in) :: Ap(*),Ai(*)
real(r8),intent(in) :: Ax(*),Az(*)
type(c_ptr) :: Symbolic
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
call s_umfpack_symbolic(m,n,Ap,Ai,Ax,Az,Symbolic,Control,Info)
end subroutine

subroutine umf4csym(m,n,Ap,Ai,Ax,Symbolic,Control,Info)
integer,intent(in) :: m,n
integer,intent(in) :: Ap(*),Ai(*)
complex(r8),intent(in) :: Ax(*)
type(c_ptr) :: Symbolic
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
call s_umfpack_symbolic(m,n,Ap,Ai,Ax,Symbolic,Control,Info)
end subroutine

subroutine umf4sym_ip(m,n,Ap,Ai,Ax,Symbolic,Control,Info)
integer,intent(in) :: m,n
integer,intent(in) :: Ap(*),Ai(*)
real(r8),intent(in) :: Ax(*)
integer(ip) :: Symbolic
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
type(c_ptr) :: c_Symbolic
c_Symbolic=c_null_ptr
call s_umfpack_symbolic(m,n,Ap,Ai,Ax,c_Symbolic,Control,Info)
Symbolic=transfer(c_Symbolic,Symbolic)
end subroutine

subroutine umf4zsym_ip(m,n,Ap,Ai,Ax,Az,Symbolic,Control,Info)
integer,intent(in) :: m,n
integer,intent(in) :: Ap(*),Ai(*)
real(r8),intent(in) :: Ax(*),Az(*)
integer(ip) :: Symbolic
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
type(c_ptr) :: c_Symbolic
c_Symbolic=c_null_ptr
call s_umfpack_symbolic(m,n,Ap,Ai,Ax,Az,c_Symbolic,Control,Info)
Symbolic=transfer(c_Symbolic,Symbolic)
end subroutine

subroutine umf4csym_ip(m,n,Ap,Ai,Ax,Symbolic,Control,Info)
integer,intent(in) :: m,n
integer,intent(in) :: Ap(*),Ai(*)
complex(r8),intent(in) :: Ax(*)
integer(ip) :: Symbolic
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
type(c_ptr) :: c_Symbolic
c_Symbolic=c_null_ptr
call s_umfpack_symbolic(m,n,Ap,Ai,Ax,c_Symbolic,Control,Info)
Symbolic=transfer(c_Symbolic,Symbolic)
end subroutine

subroutine umf4num(Ap,Ai,Ax,Symbolic,Numeric,Control,Info)
integer,intent(in) :: Ap(*),Ai(*)
real(r8),intent(in) :: Ax(*)
type(c_ptr) :: Symbolic,Numeric
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
call s_umfpack_numeric(Ap,Ai,Ax,Symbolic,Numeric,Control,Info)
end subroutine

subroutine umf4znum(Ap,Ai,Ax,Az,Symbolic,Numeric,Control,Info)
integer,intent(in) :: Ap(*),Ai(*)
real(r8),intent(in) :: Ax(*),Az(*)
type(c_ptr) :: Symbolic,Numeric
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
call s_umfpack_numeric(Ap,Ai,Ax,Az,Symbolic,Numeric,Control,Info)
end subroutine

subroutine umf4cnum(Ap,Ai,Ax,Symbolic,Numeric,Control,Info)
integer,intent(in) :: Ap(*),Ai(*)
complex(r8),intent(in) :: Ax(*)
type(c_ptr) :: Symbolic,Numeric
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
call s_umfpack_numeric(Ap,Ai,Ax,Symbolic,Numeric,Control,Info)
end subroutine

subroutine umf4num_ip(Ap,Ai,Ax,Symbolic,Numeric,Control,Info)
integer,intent(in) :: Ap(*),Ai(*)
real(r8),intent(in) :: Ax(*)
integer(ip) :: Symbolic,Numeric
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
type(c_ptr) :: c_Symbolic,c_Numeric
c_Symbolic=transfer(Symbolic,c_Symbolic)
c_Numeric=c_null_ptr
call s_umfpack_numeric(Ap,Ai,Ax,c_Symbolic,c_Numeric,Control,Info)
Numeric=transfer(c_Numeric,Numeric)
end subroutine

subroutine umf4znum_ip(Ap,Ai,Ax,Az,Symbolic,Numeric,Control,Info)
integer,intent(in) :: Ap(*),Ai(*)
real(r8),intent(in) :: Ax(*),Az(*)
integer(ip) :: Symbolic,Numeric
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
type(c_ptr) :: c_Symbolic,c_Numeric
c_Symbolic=transfer(Symbolic,c_Symbolic)
c_Numeric=c_null_ptr
call s_umfpack_numeric(Ap,Ai,Ax,Az,c_Symbolic,c_Numeric,Control,Info)
Numeric=transfer(c_Numeric,Numeric)
end subroutine

subroutine umf4cnum_ip(Ap,Ai,Ax,Symbolic,Numeric,Control,Info)
integer,intent(in) :: Ap(*),Ai(*)
complex(r8),intent(in) :: Ax(*)
integer(ip) :: Symbolic,Numeric
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
type(c_ptr) :: c_Symbolic,c_Numeric
c_Symbolic=transfer(Symbolic,c_Symbolic)
c_Numeric=c_null_ptr
call s_umfpack_numeric(Ap,Ai,Ax,c_Symbolic,c_Numeric,Control,Info)
Numeric=transfer(c_Numeric,Numeric)
end subroutine

subroutine umf4solr(sys,Ap,Ai,Ax,x,b,Numeric,Control,Info)
integer,intent(in) :: sys
integer,intent(in) :: Ap(*),Ai(*)
real(r8),intent(in) :: Ax(*)
real(r8),intent(out) :: x(*)
real(r8),intent(in) :: b(*)
type(c_ptr) :: Numeric
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
call s_umfpack_solve(sys,Ap,Ai,Ax,x,b,Numeric,Control,Info)
end subroutine

subroutine umf4zsolr(sys,Ap,Ai,Ax,Az,x,xz,b,bz,Numeric,Control,Info)
integer,intent(in) :: sys
integer,intent(in) :: Ap(*),Ai(*)
real(r8),intent(in) :: Ax(*),Az(*)
real(r8),intent(out) :: x(*),xz(*)
real(r8),intent(in) :: b(*),bz(*)
type(c_ptr) :: Numeric
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
call s_umfpack_solve(sys,Ap,Ai,Ax,Az,x,xz,b,bz,Numeric,Control,Info)
end subroutine

subroutine umf4csolr(sys,Ap,Ai,Ax,x,b,Numeric,Control,Info)
integer,intent(in) :: sys
integer,intent(in) :: Ap(*),Ai(*)
complex(r8),intent(in) :: Ax(*)
complex(r8),intent(out) :: x(*)
complex(r8),intent(in) :: b(*)
type(c_ptr) :: Numeric
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
call s_umfpack_solve(sys,Ap,Ai,Ax,x,b,Numeric,Control,Info)
end subroutine

subroutine umf4solr_ip(sys,Ap,Ai,Ax,x,b,Numeric,Control,Info)
integer,intent(in) :: sys
integer,intent(in) :: Ap(*),Ai(*)
real(r8),intent(in) :: Ax(*)
real(r8),intent(out) :: x(*)
real(r8),intent(in) :: b(*)
integer(ip) :: Numeric
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
type(c_ptr) :: c_Numeric
c_Numeric=transfer(Numeric,c_Numeric)
call s_umfpack_solve(sys,Ap,Ai,Ax,x,b,c_Numeric,Control,Info)
end subroutine

subroutine umf4zsolr_ip(sys,Ap,Ai,Ax,Az,x,xz,b,bz,Numeric,Control,Info)
integer,intent(in) :: sys
integer,intent(in) :: Ap(*),Ai(*)
real(r8),intent(in) :: Ax(*),Az(*)
real(r8),intent(out) :: x(*),xz(*)
real(r8),intent(in) :: b(*),bz(*)
integer(ip) :: Numeric
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
type(c_ptr) :: c_Numeric
c_Numeric=transfer(Numeric,c_Numeric)
call s_umfpack_solve(sys,Ap,Ai,Ax,Az,x,xz,b,bz,c_Numeric,Control,Info)
end subroutine

subroutine umf4csolr_ip(sys,Ap,Ai,Ax,x,b,Numeric,Control,Info)
integer,intent(in) :: sys
integer,intent(in) :: Ap(*),Ai(*)
complex(r8),intent(in) :: Ax(*)
complex(r8),intent(out) :: x(*)
complex(r8),intent(in) :: b(*)
integer(ip) :: Numeric
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
type(c_ptr) :: c_Numeric
c_Numeric=transfer(Numeric,c_Numeric)
call s_umfpack_solve(sys,Ap,Ai,Ax,x,b,c_Numeric,Control,Info)
end subroutine

subroutine umf4sol(sys,x,b,Numeric,Control,Info)
integer,intent(in) :: sys
real(r8),intent(out) :: x(*)
real(r8),intent(in) :: b(*)
type(c_ptr) :: Numeric
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
real(r8) :: c_Control(0:UMFPACK_CONTROL-1)
c_Control=Control
c_Control(UMFPACK_IRSTEP)=0.
! call umfpack_solve(sys,x,b,Numeric,c_Control,Info)
call s_umfpack_solve(sys,x=x,b=b,Numeric=Numeric,Control=c_Control,Info=Info)
end subroutine

subroutine umf4zsol(sys,x,xz,b,bz,Numeric,Control,Info)
integer,intent(in) :: sys
real(r8),intent(out) :: x(*),xz(*)
real(r8),intent(in) :: b(*),bz(*)
type(c_ptr) :: Numeric
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
real(r8) :: c_Control(0:UMFPACK_CONTROL-1)
c_Control=Control
c_Control(UMFPACK_IRSTEP)=0.
! call umfpack_solve(sys,x,xz,b,bz,Numeric,c_Control,Info)
call s_umfpack_solve(sys,xx=x,xz=xz,bx=b,bz=bz,Numeric=Numeric,Control=c_Control,Info=Info)
end subroutine

subroutine umf4csol(sys,x,b,Numeric,Control,Info)
integer,intent(in) :: sys
complex(r8),intent(out) :: x(*)
complex(r8),intent(in) :: b(*)
type(c_ptr) :: Numeric
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
real(r8) :: c_Control(0:UMFPACK_CONTROL-1)
c_Control=Control
c_Control(UMFPACK_IRSTEP)=0.
! call umfpack_solve(sys,x,b,Numeric,c_Control,Info)
call s_umfpack_solve(sys,xx=x,bx=b,Numeric=Numeric,Control=c_Control,Info=Info)
end subroutine

subroutine umf4sol_ip(sys,x,b,Numeric,Control,Info)
integer,intent(in) :: sys
real(r8),intent(out) :: x(*)
real(r8),intent(in) :: b(*)
integer(ip) :: Numeric
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
real(r8) :: c_Control(0:UMFPACK_CONTROL-1)
type(c_ptr) :: c_Numeric
c_Numeric=transfer(Numeric,c_Numeric)
c_Control=Control
c_Control(UMFPACK_IRSTEP)=0.
! call umfpack_solve(sys,x,b,c_Numeric,c_Control,Info)
call s_umfpack_solve(sys,x=x,b=b,Numeric=c_Numeric,Control=c_Control,Info=Info)
end subroutine

subroutine umf4zsol_ip(sys,x,xz,b,bz,Numeric,Control,Info)
integer,intent(in) :: sys
real(r8),intent(out) :: x(*),xz(*)
real(r8),intent(in) :: b(*),bz(*)
integer(ip) :: Numeric
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
real(r8) :: c_Control(0:UMFPACK_CONTROL-1)
type(c_ptr) :: c_Numeric
c_Numeric=transfer(Numeric,c_Numeric)
c_Control=Control
c_Control(UMFPACK_IRSTEP)=0.
! call umfpack_solve(sys,x,xz,b,bz,c_Numeric,c_Control,Info)
call s_umfpack_zi_solve(sys,xx=x,xz=xz,bx=b,bz=bz,Numeric=c_Numeric,Control=c_Control,Info=Info)
end subroutine

subroutine umf4csol_ip(sys,x,b,Numeric,Control,Info)
integer,intent(in) :: sys
complex(r8),intent(out) :: x(*)
complex(r8),intent(in) :: b(*)
integer(ip) :: Numeric
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
real(r8) :: c_Control(0:UMFPACK_CONTROL-1)
type(c_ptr) :: c_Numeric
c_Numeric=transfer(Numeric,c_Numeric)
c_Control=Control
c_Control(UMFPACK_IRSTEP)=0.
! call umfpack_solve(sys,x,b,c_Numeric,c_Control,Info)
call s_umfpack_ci_solve(sys,xx=x,bx=b,Numeric=c_Numeric,Control=c_Control,Info=Info)
end subroutine

subroutine umf4scal(x,b,Numeric,status)
real(r8),intent(out) :: x(*)
real(r8),intent(in) :: b(*)
type(c_ptr) :: Numeric
integer,intent(out) :: status
status=umfpack_scale_function(x,b,Numeric)
end subroutine

subroutine umf4zscal(x,xz,b,bz,Numeric,status)
real(r8),intent(out) :: x(*),xz(*)
real(r8),intent(in) :: b(*),bz(*)
type(c_ptr) :: Numeric
integer,intent(out) :: status
status=umfpack_scale_function(x,xz,b,bz,Numeric)
end subroutine

subroutine umf4cscal(x,b,Numeric,status)
complex(r8),intent(out) :: x(*)
complex(r8),intent(in) :: b(*)
type(c_ptr) :: Numeric
integer,intent(out) :: status
status=umfpack_scale_function(x,b,Numeric)
end subroutine

subroutine umf4scal_ip(x,b,Numeric,status)
real(r8),intent(out) :: x(*)
real(r8),intent(in) :: b(*)
integer(ip) :: Numeric
integer,intent(out) :: status
type(c_ptr) :: c_Numeric
c_Numeric=transfer(Numeric,c_Numeric)
status=umfpack_scale_function(x,b,c_Numeric)
end subroutine

subroutine umf4zscal_ip(x,xz,b,bz,Numeric,status)
real(r8),intent(out) :: x(*),xz(*)
real(r8),intent(in) :: b(*),bz(*)
integer(ip) :: Numeric
integer,intent(out) :: status
type(c_ptr) :: c_Numeric
c_Numeric=transfer(Numeric,c_Numeric)
status=umfpack_scale_function(x,xz,b,bz,c_Numeric)
end subroutine

subroutine umf4cscal_ip(x,b,Numeric,status)
complex(r8),intent(out) :: x(*)
complex(r8),intent(in) :: b(*)
integer(ip) :: Numeric
integer,intent(out) :: status
type(c_ptr) :: c_Numeric
c_Numeric=transfer(Numeric,c_Numeric)
status=umfpack_scale_function(x,b,c_Numeric)
end subroutine

subroutine umf4pinf(Control,Info)
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
call umfpack_report_info(Control,Info)
end subroutine

subroutine umf4zpinf(Control,Info)
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
call umfpack_report_info(Control,Info,"zi")
end subroutine

subroutine umf4cpinf(Control,Info)
real(r8) :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
call umfpack_report_info(Control,Info,"ci")
end subroutine

subroutine umf4fnum(Numeric)
type(c_ptr) :: Numeric
call umfpack_free_numeric(Numeric)
end subroutine

subroutine umf4zfnum(Numeric)
type(c_ptr) :: Numeric
call umfpack_free_numeric(Numeric,"zi")
end subroutine

subroutine umf4cfnum(Numeric)
type(c_ptr) :: Numeric
call umfpack_free_numeric(Numeric,"ci")
end subroutine

subroutine umf4fnum_ip(Numeric)
integer(ip) :: Numeric
type(c_ptr) :: c_Numeric
c_Numeric=transfer(Numeric,c_Numeric)
call umfpack_free_numeric(c_Numeric)
Numeric=transfer(c_Numeric,Numeric)
end subroutine

subroutine umf4zfnum_ip(Numeric)
integer(ip) :: Numeric
type(c_ptr) :: c_Numeric
c_Numeric=transfer(Numeric,c_Numeric)
call umfpack_free_numeric(c_Numeric,"zi")
Numeric=transfer(c_Numeric,Numeric)
end subroutine

subroutine umf4cfnum_ip(Numeric)
integer(ip) :: Numeric
type(c_ptr) :: c_Numeric
c_Numeric=transfer(Numeric,c_Numeric)
call umfpack_free_numeric(c_Numeric,"ci")
Numeric=transfer(c_Numeric,Numeric)
end subroutine

subroutine umf4fsym(Symbolic)
type(c_ptr) :: Symbolic
call umfpack_free_symbolic(Symbolic)
end subroutine

subroutine umf4zfsym(Symbolic)
type(c_ptr) :: Symbolic
call umfpack_free_symbolic(Symbolic,"zi")
end subroutine

subroutine umf4cfsym(Symbolic)
type(c_ptr) :: Symbolic
call umfpack_free_symbolic(Symbolic,"ci")
end subroutine

subroutine umf4fsym_ip(Symbolic)
integer(ip) :: Symbolic
type(c_ptr) :: c_Symbolic
c_Symbolic=transfer(Symbolic,c_Symbolic)
call umfpack_free_symbolic(c_Symbolic)
Symbolic=transfer(c_Symbolic,Symbolic)
end subroutine

subroutine umf4zfsym_ip(Symbolic)
integer(ip) :: Symbolic
type(c_ptr) :: c_Symbolic
c_Symbolic=transfer(Symbolic,c_Symbolic)
call umfpack_free_symbolic(c_Symbolic,"zi")
Symbolic=transfer(c_Symbolic,Symbolic)
end subroutine

subroutine umf4cfsym_ip(Symbolic)
integer(ip) :: Symbolic
type(c_ptr) :: c_Symbolic
c_Symbolic=transfer(Symbolic,c_Symbolic)
call umfpack_free_symbolic(c_Symbolic,"ci")
Symbolic=transfer(c_Symbolic,Symbolic)
end subroutine

subroutine umf4snum(Numeric,filenum,status)
type(c_ptr) :: Numeric
integer,intent(in) :: filenum
integer,intent(out) :: status
character(20) :: filename
write (filename,'(a,i0,a)') "n",filenum,".umf"
status=umfpack_save_numeric(Numeric,filename)
end subroutine

subroutine umf4zsnum(Numeric,filenum,status)
type(c_ptr) :: Numeric
integer,intent(in) :: filenum
integer,intent(out) :: status
character(20) :: filename
write (filename,'(a,i0,a)') "n",filenum,".umf"
status=umfpack_save_numeric(Numeric,filename,"zi")
end subroutine

subroutine umf4csnum(Numeric,filenum,status)
type(c_ptr) :: Numeric
integer,intent(in) :: filenum
integer,intent(out) :: status
character(20) :: filename
write (filename,'(a,i0,a)') "n",filenum,".umf"
status=umfpack_save_numeric(Numeric,filename,"ci")
end subroutine

subroutine umf4snum_ip(Numeric,filenum,status)
integer(ip) :: Numeric
integer,intent(in) :: filenum
integer,intent(out) :: status
type(c_ptr) :: c_Numeric
character(20) :: filename
write (filename,'(a,i0,a)') "n",filenum,".umf"
c_Numeric=transfer(Numeric,c_Numeric)
status=umfpack_save_numeric(c_Numeric,filename)
end subroutine

subroutine umf4zsnum_ip(Numeric,filenum,status)
integer(ip) :: Numeric
integer,intent(in) :: filenum
integer,intent(out) :: status
type(c_ptr) :: c_Numeric
character(20) :: filename
write (filename,'(a,i0,a)') "n",filenum,".umf"
c_Numeric=transfer(Numeric,c_Numeric)
status=umfpack_save_numeric(c_Numeric,filename,"zi")
end subroutine

subroutine umf4csnum_ip(Numeric,filenum,status)
integer(ip) :: Numeric
integer,intent(in) :: filenum
integer,intent(out) :: status
type(c_ptr) :: c_Numeric
character(20) :: filename
write (filename,'(a,i0,a)') "n",filenum,".umf"
c_Numeric=transfer(Numeric,c_Numeric)
status=umfpack_save_numeric(c_Numeric,filename,"ci")
end subroutine

subroutine umf4ssym(Symbolic,filenum,status)
type(c_ptr) :: Symbolic
integer,intent(in) :: filenum
integer,intent(out) :: status
character(20) :: filename
write (filename,'(a,i0,a)') "s",filenum,".umf"
status=umfpack_save_symbolic(Symbolic,filename)
end subroutine

subroutine umf4zssym(Symbolic,filenum,status)
type(c_ptr) :: Symbolic
integer,intent(in) :: filenum
integer,intent(out) :: status
character(20) :: filename
write (filename,'(a,i0,a)') "s",filenum,".umf"
status=umfpack_save_symbolic(Symbolic,filename,"zi")
end subroutine

subroutine umf4cssym(Symbolic,filenum,status)
type(c_ptr) :: Symbolic
integer,intent(in) :: filenum
integer,intent(out) :: status
character(20) :: filename
write (filename,'(a,i0,a)') "s",filenum,".umf"
status=umfpack_save_symbolic(Symbolic,filename,"ci")
end subroutine

subroutine umf4ssym_ip(Symbolic,filenum,status)
integer(ip) :: Symbolic
integer,intent(in) :: filenum
integer,intent(out) :: status
type(c_ptr) :: c_Symbolic
character(20) :: filename
write (filename,'(a,i0,a)') "s",filenum,".umf"
c_Symbolic=transfer(Symbolic,c_Symbolic)
status=umfpack_save_symbolic(c_Symbolic,filename)
end subroutine

subroutine umf4zssym_ip(Symbolic,filenum,status)
integer(ip) :: Symbolic
integer,intent(in) :: filenum
integer,intent(out) :: status
type(c_ptr) :: c_Symbolic
character(20) :: filename
write (filename,'(a,i0,a)') "s",filenum,".umf"
c_Symbolic=transfer(Symbolic,c_Symbolic)
status=umfpack_save_symbolic(c_Symbolic,filename,"zi")
end subroutine

subroutine umf4cssym_ip(Symbolic,filenum,status)
integer(ip) :: Symbolic
integer,intent(in) :: filenum
integer,intent(out) :: status
type(c_ptr) :: c_Symbolic
character(20) :: filename
write (filename,'(a,i0,a)') "s",filenum,".umf"
c_Symbolic=transfer(Symbolic,c_Symbolic)
status=umfpack_save_symbolic(c_Symbolic,filename,"ci")
end subroutine

subroutine umf4lnum(Numeric,filenum,status)
type(c_ptr) :: Numeric
integer,intent(in) :: filenum
integer,intent(out) :: status
character(20) :: filename
write (filename,'(a,i0,a)') "n",filenum,".umf"
status=umfpack_load_numeric(Numeric,filename)
end subroutine

subroutine umf4zlnum(Numeric,filenum,status)
type(c_ptr) :: Numeric
integer,intent(in) :: filenum
integer,intent(out) :: status
character(20) :: filename
write (filename,'(a,i0,a)') "n",filenum,".umf"
status=umfpack_load_numeric(Numeric,filename,"zi")
end subroutine

subroutine umf4clnum(Numeric,filenum,status)
type(c_ptr) :: Numeric
integer,intent(in) :: filenum
integer,intent(out) :: status
character(20) :: filename
write (filename,'(a,i0,a)') "n",filenum,".umf"
status=umfpack_load_numeric(Numeric,filename,"ci")
end subroutine

subroutine umf4lnum_ip(Numeric,filenum,status)
integer(ip) :: Numeric
integer,intent(in) :: filenum
integer,intent(out) :: status
type(c_ptr) :: c_Numeric
character(20) :: filename
write (filename,'(a,i0,a)') "n",filenum,".umf"
c_Numeric=c_null_ptr
status=umfpack_load_numeric(c_Numeric,filename)
Numeric=transfer(c_Numeric,Numeric)
end subroutine

subroutine umf4zlnum_ip(Numeric,filenum,status)
integer(ip) :: Numeric
integer,intent(in) :: filenum
integer,intent(out) :: status
type(c_ptr) :: c_Numeric
character(20) :: filename
write (filename,'(a,i0,a)') "n",filenum,".umf"
c_Numeric=c_null_ptr
status=umfpack_load_numeric(c_Numeric,filename,"zi")
Numeric=transfer(c_Numeric,Numeric)
end subroutine

subroutine umf4clnum_ip(Numeric,filenum,status)
integer(ip) :: Numeric
integer,intent(in) :: filenum
integer,intent(out) :: status
type(c_ptr) :: c_Numeric
character(20) :: filename
write (filename,'(a,i0,a)') "n",filenum,".umf"
c_Numeric=c_null_ptr
status=umfpack_load_numeric(c_Numeric,filename,"ci")
Numeric=transfer(c_Numeric,Numeric)
end subroutine

subroutine umf4lsym(Symbolic,filenum,status)
type(c_ptr) :: Symbolic
integer,intent(in) :: filenum
integer,intent(out) :: status
character(20) :: filename
write (filename,'(a,i0,a)') "s",filenum,".umf"
status=umfpack_load_symbolic(Symbolic,filename)
end subroutine

subroutine umf4zlsym(Symbolic,filenum,status)
type(c_ptr) :: Symbolic
integer,intent(in) :: filenum
integer,intent(out) :: status
character(20) :: filename
write (filename,'(a,i0,a)') "s",filenum,".umf"
status=umfpack_load_symbolic(Symbolic,filename,"zi")
end subroutine

subroutine umf4clsym(Symbolic,filenum,status)
type(c_ptr) :: Symbolic
integer,intent(in) :: filenum
integer,intent(out) :: status
character(20) :: filename
write (filename,'(a,i0,a)') "s",filenum,".umf"
status=umfpack_load_symbolic(Symbolic,filename,"ci")
end subroutine

subroutine umf4lsym_ip(Symbolic,filenum,status)
integer(ip) :: Symbolic
integer,intent(in) :: filenum
integer,intent(out) :: status
type(c_ptr) :: c_Symbolic
character(20) :: filename
write (filename,'(a,i0,a)') "s",filenum,".umf"
c_Symbolic=c_null_ptr
status=umfpack_load_symbolic(c_Symbolic,filename)
Symbolic=transfer(c_Symbolic,Symbolic)
end subroutine

subroutine umf4zlsym_ip(Symbolic,filenum,status)
integer(ip) :: Symbolic
integer,intent(in) :: filenum
integer,intent(out) :: status
type(c_ptr) :: c_Symbolic
character(20) :: filename
write (filename,'(a,i0,a)') "s",filenum,".umf"
c_Symbolic=c_null_ptr
status=umfpack_load_symbolic(c_Symbolic,filename,"zi")
Symbolic=transfer(c_Symbolic,Symbolic)
end subroutine

subroutine umf4clsym_ip(Symbolic,filenum,status)
integer(ip) :: Symbolic
integer,intent(in) :: filenum
integer,intent(out) :: status
type(c_ptr) :: c_Symbolic
character(20) :: filename
write (filename,'(a,i0,a)') "s",filenum,".umf"
c_Symbolic=c_null_ptr
status=umfpack_load_symbolic(c_Symbolic,filename,"ci")
Symbolic=transfer(c_Symbolic,Symbolic)
end subroutine

! ======================================================================
! Defined operator .umfpack.: associated procedures
! ======================================================================

function umfpack_di_operator_CSC(A,B) result (X)
type(tCSC_di),intent(in) :: A
real(c_double),intent(in) :: B(*)
real(c_double) :: X(size(A%Ap)-1)
type(c_ptr) :: Symbolic,Numeric
integer :: n
integer,parameter :: sys=UMFPACK_A
n=size(A%Ap)-1
Symbolic=c_null_ptr
Numeric=c_null_ptr
call s_umfpack_di_symbolic(n,n,A%Ap,A%Ai,A%Ax,Symbolic)
call s_umfpack_di_numeric(A%Ap,A%Ai,A%Ax,Symbolic,Numeric)
call umfpack_di_free_symbolic(Symbolic)
call s_umfpack_di_solve(sys,A%Ap,A%Ai,A%Ax,X,B,Numeric)
call umfpack_di_free_numeric(Numeric)
end function

function umfpack_zi_operator_CSC(A,B) result (X)
type(tCSC_zi),intent(in) :: A
type(tVec_zi),intent(in) :: B
type(tVec_zi) :: X
type(c_ptr) :: Symbolic,Numeric
integer :: n
integer,parameter :: sys=UMFPACK_A
n=size(A%Ap)-1
Symbolic=c_null_ptr
Numeric=c_null_ptr
allocate (X%x(n),X%z(n))
call s_umfpack_zi_symbolic(n,n,A%Ap,A%Ai,A%Ax,A%Az,Symbolic)
call s_umfpack_zi_numeric(A%Ap,A%Ai,A%Ax,A%Az,Symbolic,Numeric)
call umfpack_zi_free_symbolic(Symbolic)
call s_umfpack_zi_solve(sys,A%Ap,A%Ai,A%Ax,A%Az,X%x,X%z,B%x,B%z,Numeric)
call umfpack_zi_free_numeric(Numeric)
end function

function umfpack_ci_operator_CSC(A,B) result (X)
type(tCSC_ci),intent(in) :: A
complex(c_double),intent(in) :: B(*)
complex(c_double) :: X(size(A%Ap)-1)
type(c_ptr) :: Symbolic,Numeric
integer :: n
integer,parameter :: sys=UMFPACK_A
n=size(A%Ap)-1
Symbolic=c_null_ptr
Numeric=c_null_ptr
call s_umfpack_ci_symbolic(n,n,A%Ap,A%Ai,A%Ax,Symbolic)
call s_umfpack_ci_numeric(A%Ap,A%Ai,A%Ax,Symbolic,Numeric)
call umfpack_zi_free_symbolic(Symbolic)
call s_umfpack_ci_solve(sys,A%Ap,A%Ai,A%Ax,X,B,Numeric)
call umfpack_zi_free_numeric(Numeric)
end function

function umfpack_di_operator_CSR(A,B) result (X)
type(tCSR_di),intent(in) :: A
real(c_double),intent(in) :: B(*)
real(c_double) :: X(size(A%Ap)-1)
type(c_ptr) :: Symbolic,Numeric
integer :: n
integer,parameter :: sys=UMFPACK_Aat
n=size(A%Ap)-1
Symbolic=c_null_ptr
Numeric=c_null_ptr
call s_umfpack_di_symbolic(n,n,A%Ap,A%Ai,A%Ax,Symbolic)
call s_umfpack_di_numeric(A%Ap,A%Ai,A%Ax,Symbolic,Numeric)
call umfpack_di_free_symbolic(Symbolic)
call s_umfpack_di_solve(sys,A%Ap,A%Ai,A%Ax,X,B,Numeric)
call umfpack_di_free_numeric(Numeric)
end function

function umfpack_zi_operator_CSR(A,B) result (X)
type(tCSR_zi),intent(in) :: A
type(tVec_zi),intent(in) :: B
type(tVec_zi) :: X
type(c_ptr) :: Symbolic,Numeric
integer :: n
integer,parameter :: sys=UMFPACK_Aat
n=size(A%Ap)-1
Symbolic=c_null_ptr
Numeric=c_null_ptr
allocate (X%x(n),X%z(n))
call s_umfpack_zi_symbolic(n,n,A%Ap,A%Ai,A%Ax,A%Az,Symbolic)
call s_umfpack_zi_numeric(A%Ap,A%Ai,A%Ax,A%Az,Symbolic,Numeric)
call umfpack_zi_free_symbolic(Symbolic)
call s_umfpack_zi_solve(sys,A%Ap,A%Ai,A%Ax,A%Az,X%x,X%z,B%x,B%z,Numeric)
call umfpack_zi_free_numeric(Numeric)
end function

function umfpack_ci_operator_CSR(A,B) result (X)
type(tCSR_ci),intent(in) :: A
complex(c_double),intent(in) :: B(*)
complex(c_double) :: X(size(A%Ap)-1)
type(c_ptr) :: Symbolic,Numeric
integer :: n
integer,parameter :: sys=UMFPACK_Aat
n=size(A%Ap)-1
Symbolic=c_null_ptr
Numeric=c_null_ptr
call s_umfpack_ci_symbolic(n,n,A%Ap,A%Ai,A%Ax,Symbolic)
call s_umfpack_ci_numeric(A%Ap,A%Ai,A%Ax,Symbolic,Numeric)
call umfpack_zi_free_symbolic(Symbolic)
call s_umfpack_ci_solve(sys,A%Ap,A%Ai,A%Ax,X,B,Numeric)
call umfpack_zi_free_numeric(Numeric)
end function

function umfpack_di_operator_pCSC(A,B) result (X)
type(pCSC_di),intent(in) :: A
real(c_double),intent(in) :: B(:)
real(c_double) :: X(size(A%Ap)-1)
type(c_ptr) :: Symbolic,Numeric
integer :: n
integer,parameter :: sys=UMFPACK_A
n=size(A%Ap)-1
Symbolic=c_null_ptr
Numeric=c_null_ptr
call s_umfpack_di_symbolic(n,n,A%Ap,A%Ai,A%Ax,Symbolic)
call s_umfpack_di_numeric(A%Ap,A%Ai,A%Ax,Symbolic,Numeric)
call umfpack_di_free_symbolic(Symbolic)
call s_umfpack_di_solve(sys,A%Ap,A%Ai,A%Ax,X,B,Numeric)
call umfpack_di_free_numeric(Numeric)
end function

function umfpack_zi_operator_pCSC(A,B) result (X)
type(pCSC_zi),intent(in) :: A
type(pVec_zi),intent(in) :: B
type(tVec_zi) :: X
type(c_ptr) :: Symbolic,Numeric
integer :: n
integer,parameter :: sys=UMFPACK_A
n=size(A%Ap)-1
Symbolic=c_null_ptr
Numeric=c_null_ptr
allocate (X%x(n),X%z(n))
call s_umfpack_zi_symbolic(n,n,A%Ap,A%Ai,A%Ax,A%Az,Symbolic)
call s_umfpack_zi_numeric(A%Ap,A%Ai,A%Ax,A%Az,Symbolic,Numeric)
call umfpack_zi_free_symbolic(Symbolic)
call s_umfpack_zi_solve(sys,A%Ap,A%Ai,A%Ax,A%Az,X%x,X%z,B%x,B%z,Numeric)
call umfpack_zi_free_numeric(Numeric)
end function

function umfpack_ci_operator_pCSC(A,B) result (X)
type(pCSC_ci),intent(in) :: A
complex(c_double),intent(in) :: B(*)
complex(c_double) :: X(size(A%Ap)-1)
type(c_ptr) :: Symbolic,Numeric
integer :: n
integer,parameter :: sys=UMFPACK_A
n=size(A%Ap)-1
Symbolic=c_null_ptr
Numeric=c_null_ptr
call s_umfpack_ci_symbolic(n,n,A%Ap,A%Ai,A%Ax,Symbolic)
call s_umfpack_ci_numeric(A%Ap,A%Ai,A%Ax,Symbolic,Numeric)
call umfpack_zi_free_symbolic(Symbolic)
call s_umfpack_ci_solve(sys,A%Ap,A%Ai,A%Ax,X,B,Numeric)
call umfpack_zi_free_numeric(Numeric)
end function

function umfpack_di_operator_pCSR(A,B) result (X)
type(pCSR_di),intent(in) :: A
real(c_double),intent(in) :: B(:)
real(c_double) :: X(size(A%Ap)-1)
type(c_ptr) :: Symbolic,Numeric
integer :: n
integer,parameter :: sys=UMFPACK_Aat
n=size(A%Ap)-1
Symbolic=c_null_ptr
Numeric=c_null_ptr
call s_umfpack_di_symbolic(n,n,A%Ap,A%Ai,A%Ax,Symbolic)
call s_umfpack_di_numeric(A%Ap,A%Ai,A%Ax,Symbolic,Numeric)
call umfpack_di_free_symbolic(Symbolic)
call s_umfpack_di_solve(sys,A%Ap,A%Ai,A%Ax,X,B,Numeric)
call umfpack_di_free_numeric(Numeric)
end function

function umfpack_zi_operator_pCSR(A,B) result (X)
type(pCSR_zi),intent(in) :: A
type(pVec_zi),intent(in) :: B
type(tVec_zi) :: X
type(c_ptr) :: Symbolic,Numeric
integer :: n
integer,parameter :: sys=UMFPACK_Aat
n=size(A%Ap)-1
Symbolic=c_null_ptr
Numeric=c_null_ptr
allocate (X%x(n),X%z(n))
call s_umfpack_zi_symbolic(n,n,A%Ap,A%Ai,A%Ax,A%Az,Symbolic)
call s_umfpack_zi_numeric(A%Ap,A%Ai,A%Ax,A%Az,Symbolic,Numeric)
call umfpack_zi_free_symbolic(Symbolic)
call s_umfpack_zi_solve(sys,A%Ap,A%Ai,A%Ax,A%Az,X%x,X%z,B%x,B%z,Numeric)
call umfpack_zi_free_numeric(Numeric)
end function

function umfpack_ci_operator_pCSR(A,B) result (X)
type(pCSR_ci),intent(in) :: A
complex(c_double),intent(in) :: B(*)
complex(c_double) :: X(size(A%Ap)-1)
type(c_ptr) :: Symbolic,Numeric
integer :: n
integer,parameter :: sys=UMFPACK_Aat
n=size(A%Ap)-1
Symbolic=c_null_ptr
Numeric=c_null_ptr
call s_umfpack_ci_symbolic(n,n,A%Ap,A%Ai,A%Ax,Symbolic)
call s_umfpack_ci_numeric(A%Ap,A%Ai,A%Ax,Symbolic,Numeric)
call umfpack_zi_free_symbolic(Symbolic)
call s_umfpack_ci_solve(sys,A%Ap,A%Ai,A%Ax,X,B,Numeric)
call umfpack_zi_free_numeric(Numeric)
end function

type(pCSC_di) function make_CSC_di(Ap,Ai,Ax) result (result)
integer,target,intent(in) :: Ap(:),Ai(:)
real(r8),target,intent(in) :: Ax(:)
result%Ap=>Ap
result%Ai=>Ai
result%Ax=>Ax
end function

type(pCSC_zi) function make_CSC_zi(Ap,Ai,Ax,Az) result (result)
integer,target,intent(in) :: Ap(:),Ai(:)
real(r8),target,intent(in) :: Ax(:),Az(:)
result%Ap=>Ap
result%Ai=>Ai
result%Ax=>Ax
result%Az=>Az
end function

type(pCSC_ci) function make_CSC_ci(Ap,Ai,Ax) result (result)
integer,target,intent(in) :: Ap(:),Ai(:)
complex(r8),target,intent(in) :: Ax(:)
result%Ap=>Ap
result%Ai=>Ai
result%Ax=>Ax
end function

type(pCSR_di) function make_CSR_di(Ap,Ai,Ax) result (result)
integer,target,intent(in) :: Ap(:),Ai(:)
real(r8),target,intent(in) :: Ax(:)
result%Ap=>Ap
result%Ai=>Ai
result%Ax=>Ax
end function

type(pCSR_zi) function make_CSR_zi(Ap,Ai,Ax,Az) result (result)
integer,target,intent(in) :: Ap(:),Ai(:)
real(r8),target,intent(in) :: Ax(:),Az(:)
result%Ap=>Ap
result%Ai=>Ai
result%Ax=>Ax
result%Az=>Az
end function

type(pCSR_ci) function make_CSR_ci(Ap,Ai,Ax) result (result)
integer,target,intent(in) :: Ap(:),Ai(:)
complex(r8),target,intent(in) :: Ax(:)
result%Ap=>Ap
result%Ai=>Ai
result%Ax=>Ax
end function

type(pVec_zi) function make_Vec_zi(bx,bz) result (result)
real(r8),target,intent(in) :: bx(:),bz(:)
result%x=>bx
result%z=>bz
end function

! ======================================================================
end module
! ======================================================================
