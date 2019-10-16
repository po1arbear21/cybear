module arpack_m
  !! interface for arpack

  implicit none

  interface

    subroutine dsaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, &
                      ipntr, workd, workl, lworkl, info)
      !! Solve real symmetric sparse eigenvalue problem (EVP)

      integer,          intent(inout) :: ido
        !! Reverse communication flag (must be 0 on first call)
      character(len=1), intent(in)    :: bmat
        !! Type of B matrix: ('I': Standard EVP; 'G': Generalized EVP)
      integer,          intent(in)    :: n
        !! Dimension of the EVP
      character(len=2), intent(in)    :: which
        !! 'LA': nev largest algebraic eigenvalues
        !! 'SA': nev smallest algebraic eigenvalues
        !! 'LM': nev largest magnitude eigenvalues
        !! 'SM': nev smallest magnitude eigenvalues
        !! 'BE': nev eigenvalues from each end of spectrum
      integer,          intent(in)    :: nev
        !! Number of eigenvalues to compute
      real,             intent(in)    :: tol
        !! Relative tolerance
      real,             intent(inout) :: resid(n)
        !! Residual vector
      integer,          intent(in)    :: ncv
        !! Number of columns of the matrix v (<= n)
      integer,          intent(in)    :: ldv
        !! Leading dimension of v
      real,             intent(out)   :: v(ldv,ncv)
        !! Output Lanczos basis vectors
      integer,          intent(inout) :: iparam(11)
        !! Integer parameters
      integer,          intent(out)   :: ipntr(11)
        !! Integer pointers
      real,             intent(inout) :: workd(3*n)
        !! Workspace
      integer,          intent(in)    :: lworkl
        !! Size of additional workspace (at least ncv**2 + 8*ncv)
      real,             intent(out)   :: workl(lworkl)
        !! Additional workspace
      integer,          intent(inout) :: info
        !! Input: Set to 0 if initial residual vector is supplied
        !! Output: Success/Error flags
    end subroutine

    subroutine dnaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, &
                      ipntr, workd, workl, lworkl, info)
      !! Solve real unsymmetric sparse eigenvalue problem

      integer,          intent(inout) :: ido
        !! Reverse communication flag (must be 0 on first call)
      character(len=1), intent(in)    :: bmat
        !! Type of B matrix: ('I': Standard EVP; 'G': Generalized EVP)
      integer,          intent(in)    :: n
        !! Dimension of the EVP
      character(len=2), intent(in)    :: which
        !! 'LM': nev largest  magnitude eigenvalues
        !! 'SM': nev smallest magnitude eigenvalues
        !! 'LR': nev largest  real part eigenvalues
        !! 'SR': nev smallest real part eigenvalues
        !! 'LI': nev largest  imag part eigenvalues
        !! 'SI': nev smallest imag part eigenvalues
      integer,          intent(in)    :: nev
        !! Number of eigenvalues to compute
      real,             intent(in)    :: tol
        !! Relative tolerance
      real,             intent(inout) :: resid(n)
        !! Residual vector
      integer,          intent(in)    :: ncv
        !! Number of columns of the matrix v (<= n)
      integer,          intent(in)    :: ldv
        !! Leading dimension of v
      real,             intent(out)   :: v(ldv,ncv)
        !! Output Arnoldi basis vectors
      integer,          intent(inout) :: iparam(11)
        !! Integer parameters
      integer,          intent(out)   :: ipntr(14)
        !! Integer pointers
      real,             intent(inout) :: workd(3*n)
        !! Workspace
      integer,          intent(in)    :: lworkl
        !! Size of additional workspace (at least 3*ncv**2 + 6*ncv)
      real,             intent(out)   :: workl(lworkl)
        !! Additional workspace
      integer,          intent(inout) :: info
        !! Input: Set to 0 if initial residual vector is supplied
        !! Output: Success/Error flags
    end subroutine

    subroutine znaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, &
                      ipntr, workd, workl, lworkl, rwork, info)
      !! Solve complex unsymmetric sparse eigenvalue problem

      integer,          intent(inout) :: ido
        !! Reverse communication flag (must be 0 on first call)
      character(len=1), intent(in)    :: bmat
        !! Type of B matrix: ('I': Standard EVP; 'G': Generalized EVP)
      integer,          intent(in)    :: n
        !! Dimension of the EVP
      character(len=2), intent(in)    :: which
        !! 'LM': nev largest  magnitude eigenvalues
        !! 'SM': nev smallest magnitude eigenvalues
        !! 'LR': nev largest  real part eigenvalues
        !! 'SR': nev smallest real part eigenvalues
        !! 'LI': nev largest  imag part eigenvalues
        !! 'SI': nev smallest imag part eigenvalues
      integer,          intent(in)    :: nev
        !! Number of eigenvalues to compute
      real,             intent(in)    :: tol
        !! Relative tolerance
      complex,          intent(inout) :: resid(n)
        !! Residual vector
      integer,          intent(in)    :: ncv
        !! Number of columns of the matrix v (<= n)
      integer,          intent(in)    :: ldv
        !! Leading dimension of v
      complex,          intent(out)   :: v(ldv,ncv)
        !! Output Arnoldi basis vectors
      integer,          intent(inout) :: iparam(11)
        !! Integer parameters
      integer,          intent(out)   :: ipntr(14)
        !! Integer pointers
      complex,          intent(inout) :: workd(3*n)
        !! Workspace
      integer,          intent(in)    :: lworkl
        !! Size of additional workspace (at least 3*ncv**2 + 5*ncv)
      complex,          intent(out)   :: workl(lworkl)
        !! Additional workspace
      real,             intent(inout) :: rwork(ncv)
        !! Additional real workspace
      integer,          intent(inout) :: info
        !! Input: Set to 0 if initial residual vector is supplied
        !! Output: Success/Error flags
    end subroutine

  end interface

end module