!!----
!!----
!!---- The module CFML_EisPack is an adaption to CrysFML of the well-known
!!---- Linear Algebra collection of subroutines originally written in Fortran 77
!!---- and finally superseeded by LAPACK. Here we have adapted the translation to
!!---- Fortran 90 by John Burkardt, just transformed to a single module  and
!!---- cosmetic changes by J.Rodriguez-Carvajal (JRC) @ ILL
!!----
!!---- Modifications by JRC:
!!----    References and authors only at the beginning of the module
!!----    Special references have been kept in particular subroutines
!!----    Implicit none only at the beginning of the module
!!----    "end"'s followed by the name of the subroutine or function
!!----    Make explicit the intent of arguments
!!----    Use CFML_GlobalDeps
!!----    Introduce two global variables ERR_EisPack and ERR_EisPack_Mess for fatal error control
!!----    outside the subroutine by the calling program. The variable ierr is still in use with
!!----    different meanings depending on the subroutine.
!!----    Eliminate "go to"'s and numerical labels in HQR2
!!----
!!-
!!---- To do:  Modify to private and check which are the important subroutines to put as public
!!----         Avoid to pass the dimension in the arguments
!!----         Make interfaces to simplify the calls
!!----         Unify the style
!!----         Eliminate gotos and numerical labels
!!----
!!----  Licensing:
!!----
!!----    This code is distributed under the GNU LGPL license.
!!----
!!----  Authors:
!!----
!!----    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!!----    Klema, Moler.
!!----    FORTRAN90 version by John Burkardt
!!----
!!----  References:
!!----
!!----    James Wilkinson, Christian Reinsch,
!!----    Handbook for Automatic Computation,
!!----    Volume II, Linear Algebra, Part 2,
!!----    Springer, 1971,
!!----    ISBN: 0387054146,
!!----    LC: QA251.W67.
!!----
!!----    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!!----    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!!----    Matrix Eigensystem Routines, EISPACK Guide,
!!----    Lecture Notes in Computer Science, Volume 6,
!!----    Springer Verlag, 1976,
!!----    ISBN13: 978-3540075462,
!!----    LC: QA193.M37.
!!----
Module CFML_EisPack
  !Use CFML_GlobalDeps, only :: sp,dp

  Implicit None
  !Private
  Public

  logical,           public :: ERR_EisPack = .false.
  character(len=150),public :: ERR_EisPack_Mess = " "
  integer, parameter :: sp = selected_real_kind(6,30)
  integer, parameter :: dp = selected_real_kind(14,150)

  Contains


    !!---- subroutine bakvec ( n, t, e, m, z, ierr )
    !!----  integer,                      intent(in)     :: n,m
    !!----  real(kind=dp),dimension(n,3), intent(in)     :: t
    !!----  real(kind=dp),dimension(n),   intent(in out) :: e
    !!----  real(kind=dp),dimension(n,m), intent(in out) :: z
    !!----
    !!---- BAKVEC determines eigenvectors by reversing the FIGI transformation.
    !!----
    !!----  Discussion:
    !!----
    !!----    BAKVEC forms the eigenvectors of a nonsymmetric tridiagonal
    !!----    matrix by back transforming those of the corresponding symmetric
    !!----    matrix determined by FIGI.
    !!----
    !!----  Modified:
    !!----
    !!----    04 March 2018=
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input, real(kind=dp) T(N,3), contains the nonsymmetric matrix.  Its
    !!----    subdiagonal is stored in the positions 2:N of the first column,
    !!----    its diagonal in positions 1:N of the second column,
    !!----    and its superdiagonal in positions 1:N-1 of the third column.
    !!----    T(1,1) and T(N,3) are arbitrary.
    !!----
    !!----    Input/output, real(kind=dp) E(N).  On input, E(2:N) contains the
    !!----    subdiagonal elements of the symmetric matrix.  E(1) is arbitrary.
    !!----    On output, the contents of E have been destroyed.
    !!----
    !!----    Input, integer ::M, the number of eigenvectors to be back
    !!----    transformed.
    !!----
    !!----    Input/output, real(kind=dp) Z(N,M), contains the eigenvectors.
    !!----    On output, they have been transformed as requested.
    !!----
    !!----    Output, integer ::IERR, an error flag.
    !!----    0, for normal return,
    !!----    2*N+I, if E(I) is zero with T(I,1) or T(I-1,3) non-zero.
    !!----    In this case, the symmetric matrix is not similar
    !!----    to the original matrix, and the eigenvectors
    !!----    cannot be found by this program.
    !!----
    Subroutine bakvec ( n, t, e, m, z, ierr )
      integer,                      intent(in)     :: n,m
      real(kind=dp),dimension(n,3), intent(in)     :: t
      real(kind=dp),dimension(n),   intent(in out) :: e
      real(kind=dp),dimension(n,m), intent(in out) :: z
      !
      integer:: i,j,ierr

      ierr = 0

      if ( m == 0 ) return
      e(1) = 1.0_dp
      if ( n == 1 ) return

      do i = 2, n
        if ( e(i) == 0.0_dp ) then
          if ( t(i,1) /= 0.0_dp .or. t(i-1,3) /= 0.0_dp ) then
            ierr = 2 * n + i
            return
          end if
          e(i) = 1.0_dp
        else
          e(i) = e(i-1) * e(i) / t(i-1,3)
        end if
      end do
      do j = 1, m
        z(2:n,j) = z(2:n,j) * e(2:n)
      end do

    End Subroutine bakvec

    !!---- Subroutine balanc ( n, a, low, igh, scal )
    !!----  integer,                      intent(in)     :: n
    !!----  real(kind=dp), dimension(n,n),intent(in out) :: a
    !!----  integer,                      intent(out)    :: low, igh
    !!----  real(kind=dp), dimension(n),  intent(out)    :: scal
    !!----
    !!----  BALANC balances a real matrix before eigenvalue calculations.
    !!----
    !!----  Discussion:
    !!----
    !!----    BALANC balances a real matrix and isolates eigenvalues.
    !!----
    !!----    Suppose that the principal submatrix in rows LOW through IGH
    !!----    has been balanced, that P(J) denotes the index interchanged
    !!----    with J during the permutation step, and that the elements
    !!----    of the diagonal matrix used are denoted by D(I,J).  Then
    !!----
    !!----      scal(J) = P(J),    J = 1,...,LOW-1,
    !!----               = D(J,J),  J = LOW,...,IGH,
    !!----               = P(J)     J = IGH+1,...,N.
    !!----
    !!----    The order in which the interchanges are made is N to IGH+1,
    !!----    then 1 to LOW-1.
    !!----
    !!----    Note that 1 is returned for LOW if IGH is zero formally.
    !!----
    !!----  Modified:
    !!----
    !!----    13 February 2018, 30 November 2019
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input/output, real(kind=dp) A(N,N), the N by N matrix.  On output,
    !!----    the matrix has been balanced.
    !!----
    !!----    Output, integer ::LOW, IGH, indicate that A(I,J) is equal to
    !!----    zero if
    !!----    (1) I is greater than J and
    !!----    (2) J=1,...,LOW-1 or I=IGH+1,...,N.
    !!----
    !!----    Output, real(kind=dp) scal(N), contains information determining the
    !!----    permutations and scaling factors used.
    !!----
    Subroutine balanc ( n, a, low, igh, scal )
      integer,                      intent(in)     :: n
      real(kind=dp), dimension(n,n),intent(in out) :: a
      integer,                      intent(out)    :: low, igh
      real(kind=dp), dimension(n),  intent(out)    :: scal
      !
      real(kind=dp) :: b2,c,f,g,r,s,t,radixx
      logical :: done, noconv, swap
      integer :: i,j,k,l,m

      radixx = 16.0_dp
      b2 = radixx * radixx
      j = 0
      m = 0
      k = 1
      l = n
    !
    !  Search for rows isolating an eigenvalue and push them down.
    !
      done = .false.

      do while ( .not. done )

        do j = l, 1, -1

          swap = .true.

          do i = 1, l
            if ( i /= j ) then
              if ( a(j,i) /= 0.0_dp ) then
                swap = .false.
                exit
              end if
            end if
          end do

          if ( swap ) then

            m = l
            scal(m) = j

            if ( j /= m ) then

              do i = 1, l
                t      = a(i,j)
                a(i,j) = a(i,m)
                a(i,m) = t
              end do

              do i = k, n
                t      = a(j,i)
                a(j,i) = a(m,i)
                a(m,i) = t
              end do

            end if

            if ( l == 1 ) then
              low = k
              igh = l
              return
            end if

            l = l - 1
            if ( l < 1 ) then
              done = .true.
            end if
            exit

          else if ( j == 1 ) then
            done = .true.
            exit
          end if

        end do

      end do
      !
      !  Search for columns isolating an eigenvalue and push them left.
      !
      done = .false.

      do while ( .not. done )

        do j = k, l

          swap = .true.

          do i = k, l
            if ( i /= j ) then
              if ( a(i,j) /= 0.0_dp ) then
                swap = .false.
                exit
              end if
            end if
          end do

          if ( swap ) then

            m = k
            scal(m) = j

            if ( j /= m ) then

              do i = 1, l
                t      = a(i,j)
                a(i,j) = a(i,m)
                a(i,m) = t
              end do

              do i = k, n
                t      = a(j,i)
                a(j,i) = a(m,i)
                a(m,i) = t
              end do

            end if

            k = k + 1
            if ( l < k ) then
              done = .true.
            end if
            exit

          else

            if ( j == l ) then
              done = .true.
              exit
            end if

          end if

        end do

      end do
      !
      !  Balance the submatrix in rows K to L.
      !
      scal(k:l) = 1.0_dp
      !
      !  Iterative loop for norm reduction.
      !
      noconv = .true.

      do while ( noconv )

        noconv = .false.

        do i = k, l

          c = 0.0_dp
          r = 0.0_dp

          do j = k, l
            if ( j /= i ) then
              c = c + abs ( a(j,i) )
              r = r + abs ( a(i,j) )
            end if
          end do
          !
          !  Guard against zero C or R due to underflow.
          !
          if ( c /= 0.0_dp .and. r /= 0.0_dp ) then

            g = r / radixx
            f = 1.0_dp
            s = c + r

            do while ( c < g )
              f = f * radixx
              c = c * b2
            end do

            g = r * radixx

            do while ( g <= c )
              f = f / radixx
              c = c / b2
            end do
            !
            !  Balance.
            !
            if ( ( c + r ) / f < 0.95_dp * s ) then

              g = 1.0_dp / f
              scal(i) = scal(i) * f
              noconv = .true.

              a(i,k:n) = a(i,k:n) * g
              a(1:l,i) = a(1:l,i) * f

            end if

          end if

        end do

      end do

      low = k
      igh = l

    End Subroutine balanc

    !!---- subroutine balbak ( n, low, igh, scal, m, z )
    !!---- BALBAK determines eigenvectors by undoing the BALANC transformation.
    !!----
    !!----  Discussion:
    !!----
    !!----    BALBAK forms the eigenvectors of a real general matrix by
    !!----    back transforming those of the corresponding balanced matrix
    !!----    determined by BALANC.
    !!----
    !!----  Modified:
    !!----
    !!----    04 March 2018
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input, integer ::LOW, IGH, column indices determined by BALANC.
    !!----
    !!----    Input, real(kind=dp) scal(N), contains information determining
    !!----    the permutations and scaling factors used by BALANC.
    !!----
    !!----    Input, integer ::M, the number of columns of Z to be
    !!----    back-transformed.
    !!----
    !!----    Input/output, real(kind=dp) Z(N,M), contains the real and imaginary
    !!----    parts of the eigenvectors, which, on return, have been back-transformed.
    !!----
    Subroutine balbak ( n, low, igh, scal, m, z )

      integer,                       intent(in) :: n, low, igh,m
      real(kind=dp), dimension(n),   intent(in) :: scal
      real(kind=dp), dimension(n,m), intent(out):: z

      integer       :: i,j,k,ii
      real(kind=dp) :: t

      if ( m <= 0 ) return

      if ( igh /= low ) then
        do i = low, igh
          z(i,1:m) = z(i,1:m) * scal(i)
        end do
      end if

      do ii = 1, n
        i = ii
        if ( i < low .or. igh < i ) then
          if ( i < low ) then
            i = low - ii
          end if
          k = int ( scal(i) )
          if ( k /= i ) then
            do j = 1, m
              t      = z(i,j)
              z(i,j) = z(k,j)
              z(k,j) = t
            end do
          end if
        end if
      end do

    End Subroutine balbak

    !!---- subroutine bandr ( n, mb, a, d, e, e2, matz, z )
    !!----  integer,                        intent(in)     :: n,mb
    !!----  real(kind=dp), dimension(n,mb), intent(in out) :: a
    !!----  real(kind=dp), dimension(n),    intent(out)    :: d, e, e2
    !!----  logical.                        intent(in)     :: matz
    !!----  real(kind=dp), dimension(n,n),  intent(out)    :: z
    !!----
    !!---- BANDR reduces a symmetric band matrix to symmetric tridiagonal form.
    !!----
    !!---- Discussion:
    !!----
    !!----   BANDR reduces a real symmetric band matrix
    !!----   to a symmetric tridiagonal matrix using and optionally
    !!----   accumulating orthogonal similarity transformations.
    !!----
    !!----
    !!---- Modified:
    !!----
    !!----   04 March 2018
    !!----
    !!---- Arguments:
    !!----
    !!----   Input, integer ::N, the order of the matrix.
    !!----
    !!----   Input, integer ::MB, is the (half) band width of the matrix,
    !!----   defined as the number of adjacent diagonals, including the principal
    !!----   diagonal, required to specify the non-zero portion of the
    !!----   lower triangle of the matrix.
    !!----
    !!----   Input/output, real(kind=dp) A(N,MB).  On input, contains the lower
    !!----   triangle of the symmetric band input matrix stored as an N by MB array.
    !!----   Its lowest subdiagonal is stored in the last N+1-MB positions of the first
    !!----   column, its next subdiagonal in the last N+2-MB positions of the second
    !!----   column, further subdiagonals similarly, and finally its principal diagonal
    !!----   in the N positions of the last column.  Contents of storages not part of
    !!----   the matrix are arbitrary.  On output, A has been destroyed, except for
    !!----   its last two columns which contain a copy of the tridiagonal matrix.
    !!----
    !!----   Output, real(kind=dp) D(N), the diagonal elements of the tridiagonal
    !!----   matrix.
    !!----
    !!----   Output, real(kind=dp) E(N), the subdiagonal elements of the tridiagonal
    !!----   matrix in E(2:N).  E(1) is set to zero.
    !!----
    !!----   Output, real(kind=dp) E2(N), contains the squares of the corresponding
    !!----   elements of E.  E2 may coincide with E if the squares are not needed.
    !!----
    !!----   Input, logical MATZ, should be true if the transformation matrix is
    !!----   to be accumulated, and false otherwise.
    !!----
    !!----   Output, real(kind=dp) Z(N,N), the orthogonal transformation matrix
    !!----   produced in the reduction if MATZ is true.  Otherwise, Z is
    !!----   not referenced.
    !!----
    Subroutine bandr ( n, mb, a, d, e, e2, matz, z )
      integer,                        intent(in)     :: n,mb
      real(kind=dp), dimension(n,mb), intent(in out) :: a
      real(kind=dp), dimension(n),    intent(out)    :: d, e, e2
      logical,                        intent(in)     :: matz
      real(kind=dp), dimension(n,n),  intent(out)    :: z
      !
      real(kind=dp) ::  b1,b2,c2,f1,f2,dmin,dminrt,g,s2,u
      integer :: i1,i2,j,k,j1,j2,kr,l,m1,maxl,maxr,mr,r,ugl

      dmin = epsilon ( dmin )
      dminrt = sqrt ( dmin )
      !
      !  Initialize the diagonal scaling matrix.
      !
      d(1:n) = 1.0_dp

      if ( matz ) then
        call r8mat_identity ( n, z )
      end if

      m1 = mb - 1

      if ( m1 < 1 ) then
        d(1:n) = a(1:n,mb)
        e(1:n) = 0.0_dp
        e2(1:n) = 0.0_dp
        return
      end if

      if ( m1 /= 1 ) then

        do k = 1, n - 2

          maxr = min ( m1, n - k )

          do r = maxr, 2, -1

            kr = k + r
            mr = mb - r
            g = a(kr,mr)
            a(kr-1,1) = a(kr-1,mr+1)
            ugl = k

            do j = kr, n, m1

              j1 = j - 1
              j2 = j1 - 1

              if ( g == 0.0_dp ) then
                exit
              end if

              b1 = a(j1,1) / g
              b2 = b1 * d(j1) / d(j)
              s2 = 1.0_dp / ( 1.0_dp + b1 * b2 )

              if ( s2 < 0.5_dp ) then

                b1 = g / a(j1,1)
                b2 = b1 * d(j) / d(j1)
                c2 = 1.0_dp - s2
                d(j1) = c2 * d(j1)
                d(j) = c2 * d(j)
                f1 = 2.0_dp * a(j,m1)
                f2 = b1 * a(j1,mb)
                a(j,m1) = - b2 * ( b1 * a(j,m1) - a(j,mb) ) - f2 + a(j,m1)
                a(j1,mb) = b2 * ( b2 * a(j,mb) + f1 ) + a(j1,mb)
                a(j,mb) = b1 * ( f2 - f1 ) + a(j,mb)

                do l = ugl, j2
                  i2 = mb - j + l
                  u = a(j1,i2+1) + b2 * a(j,i2)
                  a(j,i2) = - b1 * a(j1,i2+1) + a(j,i2)
                  a(j1,i2+1) = u
                end do

                ugl = j
                a(j1,1) = a(j1,1) + b2 * g

                if ( j /= n ) then

                  maxl = min ( m1, n - j1 )

                  do l = 2, maxl
                    i1 = j1 + l
                    i2 = mb - l
                    u = a(i1,i2) + b2 * a(i1,i2+1)
                    a(i1,i2+1) = - b1 * a(i1,i2) + a(i1,i2+1)
                    a(i1,i2) = u
                  end do

                  i1 = j + m1

                  if ( i1 <= n ) then
                    g = b2 * a(i1,1)
                  end if

                end if

                if ( matz ) then

                  do l = 1, n
                    u = z(l,j1) + b2 * z(l,j)
                    z(l,j) = - b1 * z(l,j1) + z(l,j)
                    z(l,j1) = u
                  end do

                end if

              else

                u = d(j1)
                d(j1) = s2 * d(j)
                d(j) = s2 * u
                f1 = 2.0_dp * a(j,m1)
                f2 = b1 * a(j,mb)
                u = b1 * ( f2 - f1 ) + a(j1,mb)
                a(j,m1) = b2 * ( b1 * a(j,m1) - a(j1,mb) ) + f2 - a(j,m1)
                a(j1,mb) = b2 * ( b2 * a(j1,mb) + f1 ) + a(j,mb)
                a(j,mb) = u

                do l = ugl, j2
                  i2 = mb - j + l
                  u = b2 * a(j1,i2+1) + a(j,i2)
                  a(j,i2) = - a(j1,i2+1) + b1 * a(j,i2)
                  a(j1,i2+1) = u
                end do

                ugl = j
                a(j1,1) = b2 * a(j1,1) + g

                if ( j /= n ) then

                  maxl = min ( m1, n - j1 )

                  do l = 2, maxl
                    i1 = j1 + l
                    i2 = mb - l
                    u = b2 * a(i1,i2) + a(i1,i2+1)
                    a(i1,i2+1) = - a(i1,i2) + b1 * a(i1,i2+1)
                    a(i1,i2) = u
                  end do

                  i1 = j + m1

                  if ( i1 <= n ) then
                    g = a(i1,1)
                    a(i1,1) = b1 * a(i1,1)
                  end if

                end if

                if ( matz ) then

                  do l = 1, n
                    u = b2 * z(l,j1) + z(l,j)
                    z(l,j) = - z(l,j1) + b1 * z(l,j)
                    z(l,j1) = u
                  end do

                end if

              end if

            end do

          end do
         !
         !  Rescale to avoid underflow or overflow.
         !
          if ( mod ( k, 64 ) == 0 ) then

            do j = k, n

              if ( d(j) < dmin ) then

                maxl = max ( 1, mb + 1 - j )

                a(j,maxl:m1) = dminrt * a(j,maxl:m1)

                if ( j /= n ) then

                  maxl = min ( m1, n - j )

                  do l = 1, maxl
                    i1 = j + l
                    i2 = mb - l
                    a(i1,i2) = dminrt * a(i1,i2)
                  end do

                end if

                if ( matz ) then
                  z(1:n,j) = dminrt * z(1:n,j)
                end if

                a(j,mb) = dmin * a(j,mb)
                d(j) = d(j) / dmin

              end if

            end do

          end if

        end do

      end if
      !
      !  Form square root of scaling matrix.
      !
      e(2:n) = sqrt ( d(2:n) )

      if ( matz ) then

        do k = 2, n
          z(1:n,k) = z(1:n,k) * e(k)
        end do

      end if

      u = 1.0_dp

      do j = 2, n
        a(j,m1) = u * e(j) * a(j,m1)
        u = e(j)
        e2(j) = a(j,m1) * a(j,m1)
        a(j,mb) = d(j) * a(j,mb)
        d(j) = a(j,mb)
        e(j) = a(j,m1)
      end do

      d(1) = a(1,mb)
      e(1) = 0.0_dp
      e2(1) = 0.0_dp

    End Subroutine bandr

    !!----  Subroutine bandv ( n, mbw, a, e21, m, w, z, ierr )
    !!----  BANDV finds eigenvectors from eigenvalues, for a real symmetric band matrix.
    !!----
    !!----  Discussion:
    !!----
    !!----    BANDV finds those eigenvectors of a real symmetric
    !!----    band matrix corresponding to specified eigenvalues, using inverse
    !!----    iteration.
    !!----
    !!----    The routine may also be used to solve systems of linear equations with a
    !!----    symmetric or non-symmetric band coefficient matrix.
    !!----
    !!----
    !!----  Modified:
    !!----
    !!----    27 February 2018
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input, integer ::MBW, the number of columns of the array A used
    !!----    to store the band matrix.  If the matrix is symmetric, MBW is its (half)
    !!----    band width, denoted MB and defined as the number of adjacent
    !!----    diagonals, including the principal diagonal, required to
    !!----    specify the non-zero portion of the lower triangle of the
    !!----    matrix.  If the routine is being used to solve systems
    !!----    of linear equations and the coefficient matrix is not
    !!----    symmetric, it must however have the same number of adjacent
    !!----    diagonals above the main diagonal as below, and in this
    !!----    case, MBW=2*MB-1.
    !!----
    !!----    Input, real(kind=dp) A(N,MBW), the lower triangle of the symmetric
    !!----    band input matrix stored as an N by MB array.  Its lowest subdiagonal is
    !!----    stored in the last N+1-MB positions of the first column, its next
    !!----    subdiagonal in the last N+2-MB positions of the second column, further
    !!----    subdiagonals similarly, and finally its principal diagonal in the N
    !!----    positions of column MB.  If the routine is being used to solve systems
    !!----    of linear equations, and the coefficient matrix is not symmetric, A is
    !!----    N by 2*MB-1 instead, with lower triangle as above and with its first
    !!----    superdiagonal stored in the first N-1 positions of column MB+1, its
    !!----    second superdiagonal in the first N-2 positions of column MB+2, further
    !!----    superdiagonals similarly, and finally its highest superdiagonal in
    !!----    the first N+1-MB positions of the last column.  Contents of storages
    !!----    not part of the matrix are arbitrary.
    !!----
    !!----    Input, real(kind=dp) E21, specifies the ordering of the eigenvalues
    !!----    and contains 0.0 if the eigenvalues are in ascending order, or 2.0 if
    !!----    the eigenvalues are in descending order.  If the routine is being used
    !!----    to solve systems of linear equations, E21 should be set to 1.0
    !!----    if the coefficient matrix is symmetric and to -1.0 if not.
    !!----
    !!----    Input, integer ::M, the number of specified eigenvalues or the
    !!----    number of systems of linear equations.
    !!----
    !!----    Input, real(kind=dp) W(M), contains the M eigenvalues in ascending or
    !!----    descending order.  If the routine is being used to solve systems of
    !!----    linear equations (A-W(1:M)*I) * X(1:M) = B(1:M), where I is the identity
    !!----    matrix, W should be set accordingly.
    !!----
    !!----    Input/output, real(kind=dp) Z(N,M).  On input, the constant matrix
    !!----    columns B(1:M), if the routine is used to solve systems of linear
    !!----    equations.  On output, the associated set of orthogonal eigenvectors.
    !!----    Any vector which fails to converge is set to zero.  If the
    !!----    routine is used to solve systems of linear equations,
    !!----    Z contains the solution matrix columns X(1:M).
    !!----
    !!----    Output, integer ::IERR, error flag.
    !!----    0, for normal return,
    !!----    -R, if the eigenvector corresponding to the R-th eigenvalue fails to
    !!----    converge, or if the R-th system of linear equations is nearly singular.
    !!----
    Subroutine bandv ( n, mbw, a, e21, m, w, z, ierr )
      integer,                        intent(in)     :: n, mbw, m
      real(kind=dp), dimension(n,mbw),intent(in out) :: a
      real(kind=dp),                  intent(in)     :: e21
      real(kind=dp), dimension(m),    intent(in)     :: w
      real(kind=dp), dimension(n,m),  intent(in out) :: z
      integer,                        intent(out)    :: ierr
      !
      real(kind=dp) :: eps2, eps3, eps4, norm, order
      real(kind=dp) :: t,u,uk,v,x0,x1,xu
      integer       :: i,j,k,m1,m21,mb, maxj, maxk, group
      integer       :: ii,ij,ij1,its,jj,kj,kj1,r

      real(kind=dp), dimension(n*(2*mbw-1)) :: rv
      real(kind=dp), dimension(n)           :: rv6

      ierr = 0

      if ( m == 0 ) return

      rv6(1:n) = 0.0_dp

      x0 = 0.0_dp

      if ( e21 < 0.0_dp ) then
        mb = ( mbw + 1 ) / 2
      else
        mb = mbw
      end if

      m1 = mb - 1
      m21 = m1 + mb
      order = 1.0_dp - abs ( e21 )
      !
      !  Find vectors by inverse iteration.
      !
      do r = 1, m

        its = 1
        x1 = w(r)
        !
        !  Compute norm of matrix.
        !
        if ( r == 1 ) then

          norm = 0.0_dp

          do j = 1, mb

            jj = mb + 1 - j
            kj = jj + m1
            ij = 1

            v = 0.0_dp
            do i = mb + 1 - j, n
              v = v + abs ( a(i,j) )
              if ( e21 < 0.0_dp ) then
                v = v + abs ( a(ij,kj) )
                ij = ij + 1
              end if
            end do

            norm = max ( norm, v )

          end do

          if ( e21 < 0.0_dp ) then
            norm = 0.5_dp * norm
          end if
          !
          !  EPS2 is the criterion for grouping,
          !  EPS3 replaces zero pivots and equal roots are modified by eps3,
          !  EPS4 is taken very small to avoid overflow.
          !
          if ( norm == 0.0_dp ) then
            norm = 1.0_dp
          end if

          eps2 = 0.001_dp * norm * abs ( order)
          eps3 = abs ( norm ) * epsilon ( norm )
          uk = n
          uk = sqrt ( uk )
          eps4 = uk * eps3
          group = 0
         !
         !  Look for close or coincident roots.
         !
        else

          if ( eps2 <= abs ( x1 - x0 ) ) then
            group = 0
          else
            group = group + 1

            if ( order * ( x1 - x0 ) <= 0.0_dp ) then
              x1 = x0 + order * eps3
            end if
          end if

        end if
        !
        !  Expand matrix, subtract eigenvalue, and initialize vector.
        !
        do i = 1, n

          ij = i + min ( 0, i - m1 ) * n
          kj = ij + mb * n
          ij1 = kj + m1 * n

          do j = 1, m1

            if ( ij <= m1 ) then
              if ( ij <= 0 ) then
                rv(ij1) = 0.0_dp
                ij1 = ij1 + n
              end if
            else
              rv(ij) = a(i,j)
            end if

            ij = ij + n
            ii = i + j

            if ( ii <= n ) then

              jj = mb - j

              if ( e21 < 0.0_dp ) then
                ii = i
                jj = mb + j
              end if

              rv(kj) = a(ii,jj)
              kj = kj + n

            end if

          end do

          rv(ij) = a(i,mb) - x1
          rv6(i) = eps4
          if ( order == 0.0_dp ) then
            rv6(i) = z(i,r)
          end if

        end do

        if ( m1 /= 0 ) then
    !
    !  Elimination with interchanges.
    !
          do i = 1, n

            ii = i + 1
            maxk = min ( i + m1 - 1, n )
            maxj = min ( n - i, m21 - 2 ) * n

            do k = i, maxk

              kj1 = k
              j = kj1 + n
              jj = j + maxj

              do kj = j, jj, n
                rv(kj1) = rv(kj)
                kj1 = kj
              end do

              rv(kj1) = 0.0_dp

            end do

            if ( i < n ) then

              u = 0.0_dp
              maxk = min ( i + m1, n )
              maxj = min ( n - ii, m21 - 2 ) * n

              do j = i, maxk
                if ( abs ( u ) <= abs ( rv(j) ) ) then
                  u = rv(j)
                  k = j
                end if
              end do

              j = i + n
              jj = j + maxj

              if ( k /= i ) then

                kj = k

                do ij = i, jj, n
                  t      = rv(ij)
                  rv(ij) = rv(kj)
                  rv(kj) = t
                  kj = kj + n
                end do

                if ( order == 0.0_dp ) then
                  t      = rv6(i)
                  rv6(i) = rv6(k)
                  rv6(k) = t
                end if

              end if

              if ( u /= 0.0_dp ) then

                do k = ii, maxk

                  v = rv(k) / u
                  kj = k

                  do ij = j, jj, n
                    kj = kj + n
                    rv(kj) = rv(kj) - v * rv(ij)
                  end do

                  if ( order == 0.0_dp ) then
                    rv6(k) = rv6(k) - v * rv6(i)
                  end if

                end do

              end if

            end if

          end do

        end if
        !
        !  Back substitution.
        !
        do

          do i = n, 1, -1

            maxj = min ( n + 1 - i, m21 )

            if ( maxj /= 1 ) then

              ij1 = i
              j = ij1 + n
              jj = j + ( maxj - 2 ) * n

              do ij = j, jj, n
                ij1 = ij1 + 1
                rv6(i) = rv6(i) - rv(ij) * rv6(ij1)
              end do

            end if

            v = rv(i)
            !
            !  Error: nearly singular linear system.
            !
            if ( abs ( v ) < eps3 ) then
              if ( order == 0.0_dp ) then
                ierr = - r
              end if
              v = sign ( eps3, v )
            end if

            rv6(i) = rv6(i) / v

          end do

          xu = 1.0_dp

          if ( order == 0.0_dp ) then
            z(1:n,r) = rv6(1:n) * xu
            x0 = x1
            exit
          end if
          !
          !  Orthogonalize with respect to previous members of group.
          !
          do j = r - group, r - 1
            xu = dot_product ( rv6(1:n), z(1:n,j) )
            rv6(1:n) = rv6(1:n) - xu * z(1:n,j)
          end do

          norm = sum ( abs ( rv6(1:n) ) )
          !
          !  Choose a new starting vector.
          !
          if ( 0.1_dp <= norm ) then
            u = 0.0_dp
            do i = 1, n
              u = pythag ( u, rv6(i) )
            end do
            xu = 1.0_dp / u
            z(1:n,r) = rv6(1:n) * xu
            x0 = x1
            exit
          else if ( n <= its ) then
            ierr = - r
            xu = 0.0_dp
            z(1:n,r) = rv6(1:n) * xu
            x0 = x1
            exit
          else
            its = its + 1
            xu = eps4 / ( uk + 1.0_dp )
            rv6(1) = eps4
            rv6(2:n) = xu
            rv6(its) = rv6(its) - eps4 * uk
          end if

        end do

      end do

    End Subroutine bandv




    !!---- subroutine bisect ( n, eps1, d, e, e2, t1, t2, mm, m, w, ind, ierr )
    !!---- BISECT computes some eigenvalues of a real symmetric tridiagonal matrix.
    !!----
    !!---- Discussion:
    !!----
    !!----   BISECT finds those eigenvalues of a real symmetric tridiagonal matrix
    !!----   which lie in a specified interval, using bisection.
    !!----
    !!----   In the original code, the lower and upper seach bounds were
    !!----   copied, then modified, and then restored at the end.  But they
    !!----   should really be input only!
    !!----
    !!---- Modified:
    !!----
    !!----   14 March 2018
    !!----
    !!----
    !!---- Arguments:
    !!----
    !!----   Input, integer ::N, the order of the matrix.
    !!----
    !!----   Input/output, real(kind=dp) EPS1, is an absolute error tolerance for
    !!----   the computed eigenvalues.  If the input EPS1 is non-positive, it is reset
    !!----   for each submatrix to a default value, namely, minus the product of the
    !!----   relative machine precision and the 1-norm of the submatrix.
    !!----
    !!----   Input, real(kind=dp) D(N), the diagonal elements of the input matrix.
    !!----
    !!----   Input, real(kind=dp) E(N), contains in E(2:N) the subdiagonal elements
    !!----   of the matrix.  E(1) is arbitrary.
    !!----
    !!----   Input/output, real(kind=dp) E2(N).  On input, the squares of the
    !!----   corresponding elements of E.  E2(1) is arbitrary.  On output, elements of
    !!----   E2, corresponding to elements of E regarded as negligible, have been
    !!----   replaced by zero, causing the matrix to split into a direct sum of
    !!----   submatrices.  E2(1) is also set to zero.
    !!----
    !!----   Input, real(kind=dp) T1, T2, define the interval to be searched for
    !!----   eigenvalues.  If T1 is not less than T2, no eigenvalues will be found.
    !!----
    !!----   Input, integer ::MM, an upper bound for the number of
    !!----   eigenvalues in the interval.  Warning: if more than MM eigenvalues are
    !!----   determined to lie in the interval, an error return is made with no
    !!----   eigenvalues found.
    !!----
    !!----   Output, integer ::M, the number of eigenvalues determined to lie
    !!----   in (LB,UB).
    !!----
    !!----   Output, real(kind=dp) W(M), the eigenvalues in ascending order.
    !!----
    !!----   Output, integer ::IND(MM), contains in its first M positions
    !!----   the submatrix indices associated with the corresponding eigenvalues in W:
    !!----   1 for eigenvalues belonging to the first submatrix from the top, 2 for
    !!----   those belonging to the second submatrix, and so on.
    !!----
    !!----   Output, integer ::IERR, error flag.
    !!----   0, for normal return,
    !!----   3*N+1, if M exceeds MM.
    !!----
    Subroutine bisect ( n, eps1, d, e, e2, t1, t2, mm, m, w, ind, ierr )
      integer,                     intent(in)     :: n
      real(kind=dp),               intent(in out) :: eps1
      real(kind=dp), dimension(n), intent(in out) :: d, e, e2
      real(kind=dp),               intent(in)     :: t1,t2
      integer,                     intent(in)     :: mm
      integer,                     intent (out)   :: m
      real(kind=dp), dimension(mm),intent(out)    :: w
      integer,       dimension(mm),intent(out)    :: ind
      integer,                     intent(out)    :: ierr
      !
      integer                    ::i,j,k,l,m1,m2,p,q,r,s,tag
      real(kind=dp)              :: tst1,tst2,u,ub,v,x0,x1,xu,lb
      real(kind=dp),dimension(n) :: rv4, rv5


      ierr = 0
      s = 0
      tag = 0
      lb = t1
      ub = t2
      !
      !  Look for small sub-diagonal entries.
      !
      e2(1) = 0.0_dp

      do i = 2, n

        tst1 = abs ( d(i) ) + abs ( d(i-1) )
        tst2 = tst1 + abs ( e(i) )

        if ( tst2 <= tst1 ) then
          e2(i) = 0.0_dp
        end if

      end do
      !
      !  Determine the number of eigenvalues in the interval.
      !
      p = 1
      q = n

      x1 = ub
      s = sturm_sequence ( d, e, e2, n, p, q, x1 )
      m = s

      x1 = lb
      s = sturm_sequence ( d, e, e2, n, p, q, x1 )
      m = m - s

      if ( mm < m ) then
        ierr = 3 * n + 1
        return
      end if

      q = 0
      r = 0
      !
      !  Establish and process next submatrix, refining
      !  interval by the Gerschgorin bounds.
      !
      do

        if ( r == m ) return

        tag = tag + 1
        p = q + 1
        xu = d(p)
        x0 = d(p)
        u = 0.0_dp

        do q = p, n

          x1 = u
          u = 0.0_dp
          v = 0.0_dp

          if ( q < n ) then
            u = abs ( e(q+1) )
            v = e2(q+1)
          end if

          xu = min ( d(q) - ( x1 + u ), xu )
          x0 = max ( d(q) + ( x1 + u ), x0 )

          if ( v == 0.0_dp ) exit

        end do

        x1 = max ( abs ( xu ), abs ( x0 ) ) * epsilon ( x1 )
        if ( eps1 <= 0.0_dp ) eps1 = - x1
        !
        !  Check for an isolated root within interval.
        !
        if ( p == q ) then

          if ( d(p) < t1 .or. t2 <= d(p) ) then
            if ( q < n ) then
              cycle
            else
              return
            end if
          end if

          m1 = p
          m2 = p
          rv5(p) = d(p)

        else

          x1 = x1 * ( q - p + 1 )
          lb = max ( t1, xu - x1 )
          ub = min ( t2, x0 + x1 )
          x1 = lb
          s = sturm_sequence ( d, e, e2, n, p, q, x1 )
          m1 = s + 1
          x1 = ub
          s = sturm_sequence ( d, e, e2, n, p, q, x1 )
          m2 = s

          if ( m2 < m1 ) then
            if ( q < n ) then
              cycle
            else
              return
            end if
          end if
          !
          !  Find roots by bisection.
          !
          x0 = ub
          rv5(m1:m2) = ub
          rv4(m1:m2) = lb
          !
          !  Loop for the K-th eigenvalue.
          !
          k = m2

          do

            xu = lb
            do i = k, m1, -1
              if ( xu < rv4(i) ) then
                xu = rv4(i)
                exit
              end if
            end do

            x0 = min ( x0, rv5(k) )
            !
            !  Next bisection step.
            !
            do

              x1 = ( xu + x0 ) * 0.5_dp

              if ( ( x0 - xu ) <= abs ( eps1 ) ) then
                exit
              end if

              tst1 = 2.0_dp * ( abs ( xu ) + abs ( x0 ) )
              tst2 = tst1 + ( x0 - xu )
              if ( tst2 == tst1 ) then
                exit
              end if

              s = sturm_sequence ( d, e, e2, n, p, q, x1 )
              !
              !  Refine intervals.
              !
              if ( k <= s ) then
                x0 = x1
                cycle
              end if

              xu = x1

              if ( s < m1 ) then
                rv4(m1) = x1
                cycle
              end if

              rv4(s+1) = x1

              if ( x1 < rv5(s) ) then
                rv5(s) = x1
              end if

            end do
            !
            !  K-th eigenvalue found.
            !
            rv5(k) = x1
            k = k - 1

            if ( k < m1 ) then
              exit
            end if

          end do

        end if
        !
        !  Order eigenvalues tagged with their submatrix associations.
        !
        s = r
        r = r + m2 - m1 + 1
        j = 1
        k = m1

        do l = 1, r

          if ( j <= s ) then

            if ( m2 < k ) then
              exit
            end if

            if ( w(l) <= rv5(k) ) then
              j = j + 1
              cycle
            end if

            do i = l + s - j, l, -1
              w(i+1) = w(i)
              ind(i+1) = ind(i)
            end do

          end if

          w(l) = rv5(k)
          ind(l) = tag
          k = k + 1

        end do

        if ( n <= q ) exit

      end do

    End Subroutine bisect



    !!---- Subroutine bqr ( nm, n, mb, a, t, r, ierr )
    !!----
    !!----  BQR finds the smallest eigenvalue of a real symmetric band matrix.
    !!----
    !!----  Discussion:
    !!----
    !!----    BQR finds the eigenvalue of smallest magnitude of a real
    !!----    symmetric band matrix using the QR algorithm with shifts of origin.
    !!----    Consecutive calls can be made to find further eigenvalues.
    !!----
    !!----    Note that for a subsequent call, N should be replaced by N-1, but
    !!----    MB should not be altered even when it exceeds the current N.
    !!----
    !!----
    !!----  Modified:
    !!----
    !!----    20 March 2018
    !!----
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::NM, the leading dimension of the array.
    !!----    This should be the order of the full matrix.
    !!----
    !!----    Input, integer ::N, the current order of the reduced matrix.
    !!----    On first call, N is the order of the full matrix.  On each subsequent
    !!----    call, N should be reduced by 1.
    !!----
    !!----    Input, integer ::MB, the (half) band width of the matrix,
    !!----    defined as the number of adjacent diagonals, including the principal
    !!----    diagonal, required to specify the non-zero portion of the
    !!----    lower triangle of the matrix.
    !!----
    !!----    Input/output, real(kind=dp) A(N,MB).  On input, A contains the lower
    !!----    triangle of the symmetric band input matrix stored as an N by MB array.
    !!----    Its lowest subdiagonal is stored in the last N+1-MB positions of the first
    !!----    column, its next subdiagonal in the last N+2-MB positions of the
    !!----    second column, further subdiagonals similarly, and finally its principal
    !!----    diagonal in the N positions of the last column.  Contents of storages
    !!----    not part of the matrix are arbitrary.  On a subsequent call, its output
    !!----    contents from the previous call should be passed.  On output, A contains
    !!----    the transformed band matrix.  The matrix A+T*I derived from the output
    !!----    parameters is similar to the input A+T*I to within rounding errors.
    !!----    Its last row and column are null as long as IERR is zero.
    !!----
    !!----    Input/output, real(kind=dp) T.  On input, T specifies the shift (of
    !!----    eigenvalues) applied to the diagonal of A in forming the input matrix.
    !!----    What is actually determined is the eigenvalue nearest to T of A+T*I, where
    !!----    I is the identity matrix.  On a subsequent call, the output value of T
    !!----    from the previous call should be passed if the next nearest eigenvalue
    !!----    is sought.  On output, T contains the computed eigenvalue of A+T*I,
    !!----    as long as IERR is zero.
    !!----
    !!----    Input/output, real(kind=dp) R.  On input for the first call, R should
    !!----    be specified as zero, and as its output value from the previous call
    !!----    on a subsequent call.  It is used to determine when the last row and
    !!----    column of the transformed band matrix can be regarded as negligible.
    !!----    On output, R contains the maximum of its input value and the norm of the
    !!----    last column of the input matrix A.
    !!----
    !!----    Output, integer ::IERR, error flag.
    !!----    0, normal return.
    !!----    N, if the eigenvalue has not been determined after 30 iterations.
    !!----

    Subroutine bqr ( nm, n, mb, a, t, r, ierr )
      integer,                         intent(in)     :: nm, n, mb
      real(kind=dp), dimension(nm,mb), intent(in out) :: a
      real(kind=dp),                   intent(in out) :: t, r
      integer, intent(out) :: ierr


      real(kind=dp) :: f,g,q,s, scal, tst1, tst2
      integer       :: i ,ii,ik,imult,its,j,jk,jm,k,kj,kj1,kk,km,l,ll,m,m1,m2,&
                       m21,m3,m31,m4,mk,mn,mz
      real(kind=dp), dimension(2*mb*mb+4*mb-3) :: rv

      ierr = 0
      m1 = min ( mb, n )
      m = m1 - 1
      m2 = m + m
      m21 = m2 + 1
      m3 = m21 + m
      m31 = m3 + 1
      m4 = m31 + m2
      mn = m + n
      mz = mb - m1
      its = 0
    !
    !  Test for convergence.
    !
      do

        g = a(n,mb)

        if ( m == 0 ) then
          exit
        end if

        f = 0.0_dp
        do k = 1, m
          mk = k + mz
          f = f + abs ( a(n,mk) )
        end do

        if ( its == 0 .and. r < f ) then
          r = f
        end if

        tst1 = r
        tst2 = tst1 + f

        if ( tst2 <= tst1 ) then
          exit
        end if

        if ( 30 <= its ) then
          ierr = n
          return
        end if

        its = its + 1
        !
        !  Form shift from bottom 2 by 2 minor.
        !
        if ( f <= 0.25_dp * r .or. 5 <= its ) then

          f = a(n,mb-1)

          if ( f /= 0.0_dp ) then
            q = ( a(n-1,mb) - g ) / ( 2.0_dp * f )
            s = pythag ( q, 1.0_dp )
            g = g - f / ( q + sign ( s, q ) )
          end if

          t = t + g

          a(1:n,mb) = a(1:n,mb) - g

        end if

        rv(m31:m4) = 0.0_dp

        do_ii: do ii = 1, mn

          i = ii - m

          if ( n < ii .and. i <= 0 ) then

            l = max ( 1, m1 + 1 - i )

          else

            if ( n < ii ) then

              l = max ( 1, m1 + 1 - i )
            !
            !  Perform additional steps.
            !
              rv(1:m21) = 0.0_dp
              ll = min ( m1, n - ii + m1 )
              !
              !  Get row of triangular factor R.
              !
              do kk = 1, ll
                k = kk - 1
                km = k + m1
                ik = i + k
                mk = mb - k
                rv(km) = a(ik,mk)
              end do
              !
              !  Post-multiply with Householder reflections.
              !
              ll = m1
              imult = 1

            else
              !
              !  Form column of shifted matrix A-G*I.
              !
              l = max ( 1, 2 - i )

              rv(1:m3) = 0.0_dp

              do k = l, m1
                km = k + m
                mk = k + mz
                rv(km) = a(ii,mk)
              end do

              ll = min ( m, n - ii )

              do k = 1, ll
                km = k + m21
                ik = ii + k
                mk = mb - k
                rv(km) = a(ik,mk)
              end do
              !
              !  Pre-multiply with Householder reflections.
              !
              ll = m2
              imult = 0

            end if
           !
           !  Multiplication procedure.
           !

            do_140: do
                kj = m4 - m1

                do j = 1, ll

                  kj = kj + m1
                  jm = j + m3

                  if ( rv(jm) /= 0.0_dp ) then

                    f = 0.0_dp
                    do k = 1, m1
                      kj = kj + 1
                      jk = j + k - 1
                      f = f + rv(kj) * rv(jk)
                    end do
                    f = f / rv(jm)

                    kj = kj - m1

                    do k = 1, m1
                      kj = kj + 1
                      jk = j + k - 1
                      rv(jk) = rv(jk) - rv(kj) * f
                    end do

                    kj = kj - m1

                  end if

                end do

                if ( imult /= 0 ) then

                  do k = l, m1
                    mk = k + mz
                    a(i,mk) = rv(k)
                  end do

                  if ( 1 < l ) then
                    l = l - 1
                  end if

                  kj1 = m4 + l * m1

                  do j = l, m2

                    jm = j + m3
                    rv(jm) = rv(jm+1)

                    do k = 1, m1
                      kj1 = kj1 + 1
                      kj = kj1 - m1
                      rv(kj) = rv(kj1)
                    end do

                  end do

                  cycle do_ii

                end if
                !
                !  Householder reflection.
                !
                f = rv(m21)
                s = 0.0_dp
                rv(m4) = 0.0_dp
                scal = sum ( abs ( rv(m21:m3) ) )

                if ( scal /= 0.0_dp ) then

                  do k = m21, m3
                    s = s + ( rv(k) / scal )**2
                  end do

                  s = scal * scal * s
                  g = - sign ( sqrt ( s ), f )
                  rv(m21) = g
                  rv(m4) = s - f * g
                  kj = m4 + m2 * m1 + 1
                  rv(kj) = f - g

                  do k = 2, m1
                    kj = kj + 1
                    km = k + m2
                    rv(kj) = rv(km)
                  end do

                end if
                !
                !  Save column of triangular factor R.
                !
                do k = l, m1
                  km = k + m
                  mk = k + mz
                  a(ii,mk) = rv(km)
                end do

                l = max ( 1, m1 + 1 - i )

                if ( 0 < i ) then
                  !
                  !  Perform additional steps.
                  !
                  rv(1:m21) = 0.0_dp
                  ll = min ( m1, n - ii + m1 )
                  !
                  !  Get row of triangular factor R.
                  !
                  do kk = 1, ll
                    k = kk - 1
                    km = k + m1
                    ik = i + k
                    mk = mb - k
                    rv(km) = a(ik,mk)
                  end do
                  !
                  !  Post-multiply with Householder reflections.
                  !
                  ll = m1
                  imult = 1
                  cycle do_140
                end if

                exit

            end do do_140

          end if
          !
          !  Update Householder reflections.
          !
          if ( 1 < l ) then
            l = l - 1
          end if

          kj1 = m4 + l * m1

          do j = l, m2

            jm = j + m3
            rv(jm) = rv(jm+1)

            do k = 1, m1
              kj1 = kj1 + 1
              kj = kj1 - m1
              rv(kj) = rv(kj1)
            end do

          end do

        end do do_ii

      end do

      t = t + g
      a(1:n,mb) = a(1:n,mb) - g
      do k = 1, m1
        mk = k + mz
        a(n,mk) = 0.0_dp
      end do

    End Subroutine bqr


    !!---- Subroutine cbabk2 ( n, low, igh, scal, m, zr, zi )
    !!----
    !!----  CBABK2 finds eigenvectors by undoing the CBAL transformation.
    !!----
    !!----  Discussion:
    !!----
    !!----    CBABK2 forms the eigenvectors of a complex general matrix by
    !!----    back transforming those of the corresponding balanced matrix determined
    !!----    by CBAL.
    !!----
    !!----  Modified:
    !!----
    !!----    11 February 2018
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input, integer ::LOW, IGH, values determined by CBAL.
    !!----
    !!----    Input, real(kind=dp) scal(N), information determining the permutations
    !!----    and scaling factors used by CBAL.
    !!----
    !!----    Input, integer ::M, the number of eigenvectors to be back
    !!----    transformed.
    !!----
    !!----    Input/output, real(kind=dp) ZR(N,M), ZI(N,M).  On input, the real
    !!----    and imaginary parts, respectively, of the eigenvectors to be back
    !!----    transformed in their first M columns.  On output, the transformed
    !!----    eigenvectors.
    !!----
    Subroutine cbabk2 ( n, low, igh, scal, m, zr, zi )
      integer,                     intent(in)     :: n, low, igh
      real(kind=dp), dimension(n), intent(in)     :: scal
      integer,                     intent(in)     :: m
      real(kind=dp),dimension(n,m),intent(in out) :: zi
      real(kind=dp),dimension(n,m),intent(in out) :: zr

      integer ::i,ii,j,k
      real(kind=dp) :: s, t

      if ( m == 0 ) return

      if ( igh /= low ) then

        do i = low, igh

          s = scal(i)

          zr(i,1:m) = zr(i,1:m) * s
          zi(i,1:m) = zi(i,1:m) * s

        end do

      end if

      do ii = 1, n

        i = ii

        if ( i < low .or. igh < i ) then

          if ( i < low ) then
            i = low - ii
          end if

          k = int ( scal(i) )

          if ( k /= i ) then

            do j = 1, m
              t       = zr(i,j)
              zr(i,j) = zr(k,j)
              zr(k,j) = t
              t       = zi(i,j)
              zi(i,j) = zi(k,j)
              zi(k,j) = t
            end do

          end if

        end if

      end do

    End Subroutine cbabk2



    !!----  Subroutine cbal ( n, ar, ai, low, igh, scal )
    !!----
    !!----  CBAL balances a complex matrix before eigenvalue calculations.
    !!----
    !!----  Discussion:
    !!----
    !!----    CBAL balances a complex matrix and isolates eigenvalues whenever possible.
    !!----
    !!----    Suppose that the principal submatrix in rows low through igh
    !!----    has been balanced, that P(J) denotes the index interchanged
    !!----    with J during the permutation step, and that the elements
    !!----    of the diagonal matrix used are denoted by D(I,J).  Then
    !!----      scal(J) = P(J),    for J = 1,...,LOW-1
    !!----               = D(J,J)       J = LOW,...,IGH
    !!----               = P(J)         J = IGH+1,...,N.
    !!----    The order in which the interchanges are made is N to IGH+1,
    !!----    then 1 to LOW-1.
    !!----
    !!----    Note that 1 is returned for IGH if IGH is zero formally.
    !!----
    !!----  Modified:
    !!----
    !!----    24 March 2018
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input/output, real(kind=dp) AR(N,N), AI(N,N).  On input, the real and
    !!----    imaginary parts of the complex matrix to be balanced.  On output,
    !!----    the real and imaginary parts of the balanced matrix.
    !!----
    !!----    Output, integer ::LOW, IGH, are values such that AR(I,J)
    !!----    and AI(I,J) are zero if I is greater than J and either J=1,...,LOW-1 or
    !!----    I=IGH+1,...,N.
    !!----
    !!----    Output, real(kind=dp) scal(N), information determining the
    !!----    permutations and scaling factors used.
    !!----
    Subroutine cbal ( n, ar, ai, low, igh, scal )
      integer,                       intent(in)     :: n
      real(kind=dp), dimension(n,n), intent(in out) :: ar,ai
      integer,                       intent(   out) :: low, igh
      real(kind=dp), dimension(n),   intent(   out) :: scal
      !
      real(kind=dp) :: b2,c,f,g,r,radixx,s,t
      integer :: i,j,k,l,m
      logical  :: jump,jump2,noconv

      radixx = 16.0_dp
      j = 0
      m = 0

      b2 = radixx * radixx
      k = 1
      l = n

      do

        jump2 = .true.

        do j = l, 1, -1

          jump = .true.
          do i = 1, l
            if ( i /= j ) then
              if ( ar(j,i) /= 0.0_dp .or. ai(j,i) /= 0.0_dp ) then
                jump = .false.
                exit
              end if
            end if
          end do

          if ( jump ) then

            m = l

            scal(m) = j

            if ( j /= m ) then

              do i = 1, l
                t       = ar(i,j)
                ar(i,j) = ar(i,m)
                ar(i,m) = t
                t       = ai(i,j)
                ai(i,j) = ai(i,m)
                ai(i,m) = t
              end do

              do i = k, n
                t       = ar(j,i)
                ar(j,i) = ar(m,i)
                ar(m,i) = t
                t       = ai(j,i)
                ai(j,i) = ai(m,i)
                ai(m,i) = t
              end do

            end if

            if ( l == 1 ) then
              low = k
              igh = l
              return
            end if

            l = l - 1

            jump2 = .false.
            exit

          end if

        end do

        if ( jump2 ) exit

      end do
      !
      !  Search for columns isolating an eigenvalue and push them left.
      !
      do

        jump2 = .true.

        do j = k, l

          jump = .true.
          do i = k, l
            if ( i /= j ) then
              if ( ar(i,j) /= 0.0_dp .or. ai(i,j) /= 0.0_dp ) then
                jump = .false.
                exit
              end if
            end if
          end do

          if ( jump ) then

            m = k

            scal(m) = j

            if ( j /= m ) then

              do i = 1, l
                t       = ar(i,j)
                ar(i,j) = ar(i,m)
                ar(i,m) = t
                t       = ai(i,j)
                ai(i,j) = ai(i,m)
                ai(i,m) = t
              end do

              do i = k, n
                t       = ar(j,i)
                ar(j,i) = ar(m,i)
                ar(m,i) = t
                t       = ai(j,i)
                ai(j,i) = ai(m,i)
                ai(m,i) = t
              end do

            end if

            k = k + 1
            jump2 = .false.
            exit

          end if

        end do

        if ( jump2 ) exit

      end do
      !
      !  Now balance the submatrix in rows k to l.
      !
      scal(k:l) = 1.0_dp
      !
      !  Iterative loop for norm reduction.
      !
      do

        noconv = .false.

        do i = k, l

          c = 0.0_dp
          r = 0.0_dp

          do j = k, l
            if ( j /= i ) then
              c = c + abs ( ar(j,i) ) + abs ( ai(j,i) )
              r = r + abs ( ar(i,j) ) + abs ( ai(i,j) )
            end if
          end do
          !
          !  Guard against zero C or R due to underflow.
          !
          if ( c /= 0.0_dp .and. r /= 0.0_dp ) then

            g = r / radixx
            f = 1.0_dp
            s = c + r

            do while ( c < g )
              f = f * radixx
              c = c * b2
            end do

            g = r * radixx

            do while  ( g <= c )
              f = f / radixx
              c = c / b2
            end do
            !
            !  Now balance.
            !
            if ( ( c + r ) / f < 0.95_dp * s ) then

              g = 1.0_dp / f
              scal(i) = scal(i) * f
              noconv = .true.

              ar(i,k:n) = ar(i,k:n) * g
              ai(i,k:n) = ai(i,k:n) * g

              ar(1:l,i) = ar(1:l,i) * f
              ai(1:l,i) = ai(1:l,i) * f

            end if

          end if

        end do

        if ( .not. noconv ) exit

      end do

      low = k
      igh = l

    End Subroutine cbal

    !!----  Subroutine cdiv ( ar, ai, br, bi, cr, ci )
    !!----
    !!----  CDIV emulates complex division, using real arithmetic.
    !!----
    !!----  Discussion:
    !!----
    !!----    CDIV performs complex division:
    !!----
    !!----      (CR,CI) = (AR,AI) / (BR,BI)
    !!----
    !!----  Modified:
    !!----
    !!----    18 October 2009
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, real(kind=dp) AR, AI, the real and imaginary parts of
    !!----    the numerator.
    !!----
    !!----    Input, real(kind=dp) BR, BI, the real and imaginary parts of
    !!----    the denominator.
    !!----
    !!----    Output, real(kind=dp) CR, CI, the real and imaginary parts of
    !!----    the result.
    !!----
    Subroutine cdiv ( ar, ai, br, bi, cr, ci )
      real(kind=dp), intent(in)  :: ar, ai, br, bi
      real(kind=dp), intent(out) :: cr, ci

      real(kind=dp) :: ais,ars,bis,brs,s

      s = abs ( br ) + abs ( bi )

      ars = ar / s
      ais = ai / s
      brs = br / s
      bis = bi / s

      s = brs * brs + bis * bis
      cr = ( ars * brs + ais * bis ) / s
      ci = ( ais * brs - ars * bis ) / s

    End Subroutine cdiv



    !!---- Subroutine cg_lr ( n, ar, ai, wr, wi, matz, zr, zi, ierr )
    !!----
    !!----  CG_LR gets eigenvalues and eigenvectors of a complex general matrix.
    !!----
    !!----  Discussion:
    !!----
    !!----    CG_LR calls EISPACK routines to find the eigenvalues and eigenvectors
    !!----    of a complex general matrix, using elementary transformations.
    !!----
    !!----  Modified:
    !!----
    !!----    04 March 2018
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input/output, real(kind=dp) AR(N,N), AI(N,N).  On input, the real and
    !!----    imaginary parts of the complex matrix.  On output, AR and AI
    !!----    have been overwritten by other information.
    !!----
    !!----    Output, real(kind=dp) WR(N), WI(N), the real and imaginary parts
    !!----    of the eigenvalues.
    !!----
    !!----    Input, logical MATZ, is false if only eigenvalues are desired,
    !!----    and true if both eigenvalues and eigenvectors are desired.
    !!----
    !!----    Output, real(kind=dp) ZR(N,N), ZI(N,N), the real and imaginary parts,
    !!----    respectively, of the eigenvectors, if MATZ is not zero.
    !!----
    !!----    Output, integer ::IERR, an error completion code described in
    !!----    the documentation for COMLR and COMLR2.  The normal completion code
    !!----    is zero.
    !!----
    Subroutine cg_lr ( n, ar, ai, wr, wi, matz, zr, zi, ierr )
      integer,                       intent(in)     :: n
      real(kind=dp), dimension(n,n), intent(in out) :: ar, ai
      real(kind=dp), dimension(n,n), intent(   out) :: wr, wi
      logical,                       intent(in)     :: matz
      real(kind=dp), dimension(n,n), intent(   out) :: zr, zi
      integer,                       intent(   out) :: ierr
      !
      real(kind=dp), dimension(n) :: fv1
      integer,       dimension(n) :: inter
      integer ::is1,is2

      call cbal ( n, ar, ai, is1, is2, fv1 )

      inter(1:n) = 0
      call comhes ( n, is1, is2, ar, ai, inter )

      if ( .not. matz ) then

        call comlr ( n, is1, is2, ar, ai, wr, wi, ierr )

        if ( ierr /= 0 ) then
            ERR_EisPack=.true.
            ERR_EisPack_Mess=" CG_LR - Fatal error!: Nonzero error return from COMLR."
            return
        end if

      else

        call comlr2 ( n, is1, is2, inter, ar, ai, wr, wi, zr, zi, ierr )

        if ( ierr /= 0 ) then
            ERR_EisPack=.true.
            ERR_EisPack_Mess=" CG_LR - Fatal error!: Nonzero error return from COMLR2."
            return
        end if

        call cbabk2 ( n, is1, is2, fv1, n, zr, zi )

      end if

    End Subroutine cg_lr



    !!----  Subroutine cg_qr ( n, ar, ai, wr, wi, matz, zr, zi, ierr )
    !!----
    !!----  CG_QR gets eigenvalues and eigenvectors of a complex general matrix.
    !!----
    !!----  Discussion:
    !!----
    !!----    CG_QR calls EISPACK routines to find the eigenvalues and eigenvectors
    !!----    of a complex general matrix, using unitary transformations.
    !!----
    !!----  Modified:
    !!----
    !!----    09 February 2018
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input/output, real(kind=dp) AR(N,N), AI(N,N).  On input, the real and
    !!----    imaginary parts of the complex matrix.  On output, AR and AI
    !!----    have been overwritten by other information.
    !!----
    !!----    Output, real(kind=dp) WR(N), WI(N), the real and imaginary parts
    !!----    of the eigenvalues.
    !!----
    !!----    Input, logical MATZ, is false if only eigenvalues are desired,
    !!----    and true if both eigenvalues and eigenvectors are desired.
    !!----
    !!----    Output, real(kind=dp) ZR(N,N), ZI(N,N), the real and imaginary parts,
    !!----    respectively, of the eigenvectors, if MATZ is not zero.
    !!----
    !!----    Output, integer ::IERR, an error completion code described in
    !!----    the documentation for COMQR and COMQR2.  The normal completion code
    !!----    is zero.
    !!----
    Subroutine cg_qr ( n, ar, ai, wr, wi, matz, zr, zi, ierr )
      integer,                       intent(in)     :: n
      real(kind=dp), dimension(n,n), intent(in out) :: ar, ai
      real(kind=dp), dimension(n,n), intent(   out) :: wr, wi
      logical,                       intent(in)     :: matz
      real(kind=dp), dimension(n,n), intent(   out) :: zr, zi
      integer,                       intent(   out) :: ierr

      real(kind=dp), dimension(n) :: fv1,fv2,fv3
      integer :: is1, is2

      call cbal ( n, ar, ai, is1, is2, fv1 )

      call corth ( n, is1, is2, ar, ai, fv2, fv3 )

      if ( .not. matz ) then

        call comqr ( n, is1, is2, ar, ai, wr, wi, ierr )

        if ( ierr /= 0 ) then
            ERR_EisPack=.true.
            ERR_EisPack_Mess=" CG_QR - Fatal error!: Nonzero error return from COMQR."
            return
        end if

      else

        call comqr2 ( n, is1, is2, fv2, fv3, ar, ai, wr, wi, zr, zi, ierr )

        if ( ierr /= 0 ) then
            ERR_EisPack=.true.
            ERR_EisPack_Mess=" CG_QR - Fatal error!: Nonzero error return from COMQR2."
            return
        end if

        call cbabk2 ( n, is1, is2, fv1, n, zr, zi )

      end if

    End Subroutine cg_qr



    !!----  Subroutine ch ( n, ar, ai, w, matz, zr, zi, ierr )
    !!----
    !!----  CH gets eigenvalues and eigenvectors of a complex Hermitian matrix.
    !!----
    !!----  Discussion:
    !!----
    !!----    CH finds the eigenvalues and eigenvectors of a complex hermitian matrix.
    !!----
    !!----
    !!----  Modified:
    !!----
    !!----    09 February 2018
    !!----
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input/output, real(kind=dp) AR(N,N), AI(N,N).  On input, the real and
    !!----    imaginary parts of the complex matrix.  On output, AR and AI
    !!----    have been overwritten by other information.
    !!----
    !!----    Output, real(kind=dp) W(N), the eigenvalues in ascending order.
    !!----
    !!----    Input, logical MATZ, is false if only eigenvalues are desired,
    !!----    and true if both eigenvalues and eigenvectors are desired.
    !!----
    !!----    Output, real(kind=dp) ZR(N,N), ZI(N,N), the real and imaginary parts,
    !!----    respectively, of the eigenvectors, if MATZ is true.
    !!----
    !!----    Output, integer ::IERR, an error completion code described in
    !!----    the documentation for TQLRAT and TQL2.  The normal completion code is zero.
    !!----
    Subroutine ch ( n, ar, ai, w, matz, zr, zi, ierr )
      integer,                       intent(in)     :: n
      real(kind=dp), dimension(n,n), intent(in out) :: ar, ai
      real(kind=dp), dimension(n,n), intent(   out) :: w
      logical,                       intent(in)     :: matz
      real(kind=dp), dimension(n,n), intent(   out) :: zr, zi
      integer,                       intent(   out) :: ierr

      real(kind=dp), dimension(2,n) :: fm1
      real(kind=dp), dimension(n)   :: fv1, fv2

      call htridi ( n, ar, ai, w, fv1, fv2, fm1 )

      if ( .not. matz ) then

        call tqlrat ( n, w, fv2, ierr )

        if ( ierr /= 0 ) then
            ERR_EisPack=.true.
            ERR_EisPack_Mess=" CH - Fatal error!: Nonzero error return from TLQRAT"
            return
        end if

      else

        call r8mat_identity ( n, zr )

        call tql2 ( n, w, fv1, zr, ierr )

        if ( ierr /= 0 ) then
            ERR_EisPack=.true.
            ERR_EisPack_Mess=" CH - Fatal error!: Nonzero error return from TLQ2"
            return
        end if

        call htribk ( n, ar, ai, fm1, n, zr, zi )

      end if

    end Subroutine ch



    !!----  Subroutine ch3 ( n, a, d, matz, zr, zi, ierr )
    !!----
    !!----  CH3 gets eigenvalues and eigenvectors of a complex Hermitian matrix.
    !!----
    !!----  Discussion:
    !!----
    !!----    CH3 finds the eigenvalues and eigenvectors of a complex hermitian matrix.
    !!----
    !!----  Modified:
    !!----
    !!----    05 February 2018
    !!----
    !!----  Author:
    !!----
    !!----    John Burkardt
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input/output, real(kind=dp) A(N,N).  On input, the lower triangle of
    !!----    the complex hermitian input matrix.  The real parts of the matrix elements
    !!----    are stored in the full lower triangle of A, and the imaginary parts are
    !!----    stored in the transposed positions of the strict upper triangle of A.  No
    !!----    storage is required for the zero imaginary parts of the diagonal elements.
    !!----    On output, A contains information about the unitary transformations
    !!----    used in the reduction.
    !!----
    !!----    Output, real(kind=dp) D(N), the eigenvalues in ascending order.
    !!----
    !!----    Input, logical MATZ, is false if only eigenvalues are desired,
    !!----    and true if both eigenvalues and eigenvectors are desired.
    !!----
    !!----    Output, real(kind=dp) ZR(N,N), ZI(N,N), the real and imaginary parts,
    !!----    respectively, of the eigenvectors, if MATZ is true.
    !!----
    !!----    Output, integer ::IERR, an error completion code described in
    !!----    the documentation for TQLRAT and TQL2.  The normal completion code is zero.
    !!----
    Subroutine ch3 ( n, a, d, matz, zr, zi, ierr )
      integer,                       intent(in)     :: n
      real(kind=dp), dimension(n,n), intent(in out) :: a
      real(kind=dp), dimension(n),   intent(   out) :: d
      logical,                       intent(in)     :: matz
      real(kind=dp), dimension(n,n), intent(   out) :: zr, zi
      integer,                       intent(   out) :: ierr

      real(kind=dp), dimension(n)   :: e, e2
      real(kind=dp), dimension(2,n) :: tau

      call htrid3 ( n, a, d, e, e2, tau )

      if ( .not. matz ) then

        call tqlrat ( n, d, e2, ierr )

        if ( ierr /= 0 ) then
            ERR_EisPack=.true.
            ERR_EisPack_Mess=" CH3 - Fatal error!: Nonzero error return from TQLRAT."
            return
        end if

      else

        call r8mat_identity ( n, zr )

        call tql2 ( n, d, e, zr, ierr )

        if ( ierr /= 0 ) then
            ERR_EisPack=.true.
            ERR_EisPack_Mess=" CH3 - Fatal error!: Nonzero error return from TQL2."
            return
        end if

        call htrib3 ( n, a, tau, n, zr, zi )

      end if

    End Subroutine ch3


    !!----  Subroutine cinvit ( n, ar, ai, wr, wi, select, mm, m, zr, zi, ierr )
    !!----
    !!----  CINVIT gets eigenvectors from eigenvalues for a complex Hessenberg matrix.
    !!----
    !!----  Discussion:
    !!----
    !!----    CINVIT finds those eigenvectors of a complex upper Hessenberg matrix
    !!----    corresponding to specified eigenvalues, using inverse iteration.
    !!----
    !!----  Modified:
    !!----
    !!----    11 February 2018
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input, real(kind=dp) AR(N,N), AI(N,N), the real and imaginary parts of
    !!----    the complex Hessenberg matrix.
    !!----
    !!----    Input/output, real(kind=dp) WR(N), WI(N).  On input, the real and
    !!----    imaginary parts of the eigenvalues of the matrix.  The eigenvalues must
    !!----    be stored in a manner identical to that of COMLR, which
    !!----    recognizes possible splitting of the matrix.  On output, WR may have been
    !!----    altered since close eigenvalues are perturbed slightly in searching for
    !!----    independent eigenvectors.
    !!----
    !!----    Input, logical SELECT(N), specifies the eigenvectors to be found.  The
    !!----    eigenvector corresponding to the J-th eigenvalue is specified by
    !!----    setting SELECT(J) to TRUE.
    !!----
    !!----    Input, integer ::MM, an upper bound for the number of
    !!----    eigenvectors to be found.
    !!----
    !!----    Output, integer ::M, the number of eigenvectors actually found.
    !!----
    !!----    Output, real(kind=dp) ZR(N,MM), ZI(N,MM), the real and imaginary parts
    !!----    of the eigenvectors.  The eigenvectors are normalized so that the
    !!----    component of largest magnitude is 1.
    !!----    Any vector which fails the acceptance test is set to zero.
    !!----
    !!----    Output, integer ::IERR, error flag.
    !!----    0, for normal return,
    !!----    -(2*N+1), if more than MM eigenvectors have been specified,
    !!----    -K, if the iteration corresponding to the K-th value fails,
    !!----    -(N+K), if both error situations occur.
    !!----
    Subroutine cinvit ( n, ar, ai, wr, wi, select, mm, m, zr, zi, ierr )
      integer,                       intent(in)     :: n
      real(kind=dp), dimension(n,n), intent(in)     :: ar,ai
      real(kind=dp), dimension(n),   intent(in out) :: wr,wi
      logical,       dimension(n),   intent(in)     :: select
      integer,                       intent(in)     :: mm
      integer,                       intent(   out) :: m
      real(kind=dp), dimension(n,mm),intent(   out) :: zr, zi
      integer,                       intent(   out) :: ierr
      !
      real(kind=dp) :: eps3, growto,ilambd,norm,normv,rlambd,t,ukroot,x,y
      integer       :: s,uk,i,its,j,k,mp
      logical :: repeat
      real(kind=dp), dimension(n,n) :: rm1,rm2
      real(kind=dp), dimension(n)   :: rv1,rv2

      ierr = 0
      uk = 0
      s = 1

      do k = 1, n

        if ( .not. select(k) ) then
          cycle
        end if

        if ( mm < s ) then
          if ( ierr /= 0 ) then
            ierr = ierr - n
          end if
          if ( ierr == 0 ) then
            ierr = - ( 2 * n + 1 )
          end if
          m = s - 1
          return
        end if

        if ( uk < k ) then
          !
          !  Check for possible splitting.
          !
          do uk = k, n
            if ( uk == n ) then
              exit
            end if
            if ( ar(uk+1,uk) == 0.0_dp .and. ai(uk+1,uk) == 0.0_dp ) then
              exit
            end if
          end do
          !
          !  Compute infinity norm of leading UK by UK (Hessenberg) matrix.
          !
          norm = 0.0_dp
          do i = 1, uk
            x = 0.0_dp
            do j = max ( 1, i - 1 ), uk
              x = x + pythag ( ar(i,j), ai(i,j) )
            end do
            norm = max ( norm, x )
          end do
          !
          !  EPS3 replaces zero pivot in decomposition
          !  and close roots are modified by EPS3.
          !
          if ( norm == 0.0_dp ) then
            norm = 1.0_dp
          end if

          eps3 = abs ( norm ) * epsilon ( eps3 )
          !
          !  GROWTO is the criterion for growth.
          !
          ukroot = real ( uk, kind = 8 )
          ukroot = sqrt ( ukroot )
          growto = 0.1_dp / ukroot

        end if

        rlambd = wr(k)
        ilambd = wi(k)
        !
        !  Perturb eigenvalue if it is close to any previous eigenvalue.
        !
        if ( 1 < k ) then

          do

            repeat = .false.

            do i = k - 1, 1, -1
              if ( select(i) .and. &
                   abs ( wr(i) - rlambd ) < eps3 .and. &
                   abs ( wi(i) - ilambd ) < eps3 ) then
                rlambd = rlambd + eps3
                repeat = .true.
                exit
              end if
            end do

            if ( .not. repeat ) then
              wr(k) = rlambd
              exit
            end if

          end do

        end if

        do i = 1, uk
          do j = 1, i - 2
            rm1(i,j) = 0.0_dp
            rm2(i,j) = 0.0_dp
          end do
          do j = max ( i - 1, 1 ), uk
            rm1(i,j) = ar(i,j)
            rm2(i,j) = ai(i,j)
          end do
          rm1(i,i) = rm1(i,i) - rlambd
          rm2(i,i) = rm2(i,i) - ilambd
          rv1(i) = eps3
        end do
        !
        !  Triangular decomposition with interchanges, replacing zero pivots by eps3.
        !
        do i = 2, uk

          mp = i - 1

          if ( pythag ( rm1(mp,mp), rm2(mp,mp) ) < &
               pythag ( rm1(i,mp), rm2(i,mp) ) ) then

            do j = i - 1, uk
              t         = rm1(i,j)
              rm1(i,j)  = rm1(mp,j)
              rm1(mp,j) = t
              t         = rm2(i,j)
              rm2(i,j)  = rm2(mp,j)
              rm2(mp,j) = t
            end do

          end if

          if ( rm1(mp,mp) == 0.0_dp .and. rm2(mp,mp) == 0.0_dp ) then
            rm1(mp,mp) = eps3
          end if

          call cdiv ( rm1(i,mp), rm2(i,mp), rm1(mp,mp), rm2(mp,mp), x, y )

          if ( x /= 0.0_dp .or. y /= 0.0_dp ) then

            do j = i, uk
              rm1(i,j) = rm1(i,j) - x * rm1(mp,j) + y * rm2(mp,j)
              rm2(i,j) = rm2(i,j) - x * rm2(mp,j) - y * rm1(mp,j)
            end do

          end if

        end do

        if ( rm1(uk,uk) == 0.0_dp .and. rm2(uk,uk) == 0.0_dp ) then
          rm1(uk,uk) = eps3
        end if

        its = 0
        !
        !  Back substitution.
        !
        do

          do i = uk, 1, -1

            x = rv1(i)
            y = 0.0_dp
            do j = i + 1, uk
              x = x - rm1(i,j) * rv1(j) + rm2(i,j) * rv2(j)
              y = y - rm1(i,j) * rv2(j) - rm2(i,j) * rv1(j)
            end do

            call cdiv ( x, y, rm1(i,i), rm2(i,i), rv1(i), rv2(i) )

          end do
          !
          !  Acceptance test for eigenvector and normalization.
          !
          its = its + 1
          norm = 0.0_dp
          normv = 0.0_dp

          do i = 1, uk
            x = pythag ( rv1(i), rv2(i) )
            if ( normv < x ) then
              normv = x
              j = i
            end if
            norm = norm + x
          end do
          !
          !  Accept vector.
          !
          if ( growto <= norm ) then

            x = rv1(j)
            y = rv2(j)

            do i = 1, uk
              call cdiv ( rv1(i), rv2(i), x, y, zr(i,s), zi(i,s) )
            end do

            if ( uk /= n ) then
              j = uk + 1
              zr(j:n,s) = 0.0_dp
              zi(j:n,s) = 0.0_dp
            end if
            s = s + 1
            exit
           !
           !  Choose a new starting vector.
           !
          else if ( its < uk ) then

            x = ukroot
            y = eps3 / ( x + 1.0_dp )

            rv1(1) = eps3
            rv1(2:uk) = y

            j = uk - its + 1
            rv1(j) = rv1(j) - eps3 * x
           !
           !  Error: unaccepted eigenvector.
           !
          else

            j = 1
            ierr = - k
           !
           !  Set remaining vector components to zero.
           !
            zr(j:n,s) = 0.0_dp
            zi(j:n,s) = 0.0_dp
            s = s + 1
            exit
          end if

        end do

      end do

      m = s - 1

    End Subroutine cinvit



    !!----  Subroutine combak ( n, low, igh, ar, ai, inter, m, zr, zi )
    !!----
    !!----  COMBAK determines eigenvectors by undoing the COMHES transformation.
    !!----
    !!----  Discussion:
    !!----
    !!----    COMBAK forms the eigenvectors of a complex general matrix by
    !!----    back transforming those of the corresponding upper Hessenberg matrix
    !!----    determined by COMHES.
    !!----
    !!----  Modified:
    !!----
    !!----    04 March 2018
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input, integer ::LOW, IGH, are determined by the balancing
    !!----    routine CBAL.  If CBAL is not used, set LOW = 1 and IGH = to the order
    !!----    of the matrix.
    !!----
    !!----    Input, real(kind=dp) AR(N,IGH), AI(N,IGH), the multipliers which
    !!----    were used in the reduction by COMHES in their lower triangles below
    !!----    the subdiagonal.
    !!----
    !!----    Input, integer ::INTER(IGH), information on the rows and
    !!----    columns interchanged in the reduction by COMHES.
    !!----
    !!----    Input, integer ::M, the number of eigenvectors to be back
    !!----    transformed.
    !!----
    !!----    Input/output, real(kind=dp) ZR(N,M), ZI(N,M).  On input, the real
    !!----    and imaginary parts of the eigenvectors to be back transformed.  On
    !!----    output, the real and imaginary parts of the transformed eigenvectors.
    !!----
    Subroutine combak ( n, low, igh, ar, ai, inter, m, zr, zi )
      integer,                         intent(in)     :: n,low, igh
      real(kind=dp), dimension(n,igh), intent(in)     :: ar,ai
      integer,       dimension(igh),   intent(in)     :: inter
      integer,                         intent(in)     :: m
      real(kind=dp), dimension(n,m),   intent(in out) :: zr, zi

      integer       :: i,j,mp
      real(kind=dp) :: t,xi,xr

      if ( m == 0 ) return
      if ( igh - 1 < low + 1 ) return

      do mp = igh - 1, low + 1, -1

        do i = mp + 1, igh

          xr = ar(i,mp-1)
          xi = ai(i,mp-1)

          if ( xr /= 0.0_dp .or. xi /= 0.0_dp ) then
            zr(i,1:m) = zr(i,1:m) + xr * zr(mp,1:m) - xi * zi(mp,1:m)
            zi(i,1:m) = zi(i,1:m) + xr * zi(mp,1:m) + xi * zr(mp,1:m)
          end if

        end do

        i = inter(mp)

        if ( i /= mp ) then

          do j = 1, m
            t        = zr(i,j)
            zr(i,j)  = zr(mp,j)
            zr(mp,j) = t
            t        = zi(i,j)
            zi(i,j)  = zi(mp,j)
            zi(mp,j) = t
          end do

        end if

      end do

    End Subroutine combak


    !!----  Subroutine comhes ( n, low, igh, ar, ai, inter )
    !!----
    !!----  COMHES transforms a complex general matrix to upper Hessenberg form.
    !!----
    !!----  Discussion:
    !!----
    !!----    COMHES is given a complex general matrix and reduces a submatrix in rows
    !!----    and columns LOW through IGH to upper Hessenberg form by
    !!----    stabilized elementary similarity transformations.
    !!----
    !!----  Modified:
    !!----
    !!----    04 March 2018
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input, integer ::LOW, IGH, are determined by the balancing
    !!----    routine CBAL.  If CBAL is not used, set LOW = 1 and IGH = N.
    !!----
    !!----    Input/output, real(kind=dp) AR(N,N), AI(N,N).  On input, the real and
    !!----    imaginary parts of the complex input matrix.  On output, the real and
    !!----    imaginary parts of the Hessenberg matrix.  The multipliers which were
    !!----    used in the reduction are stored in the remaining triangles under the
    !!----    Hessenberg matrix.
    !!----
    !!----    Output, integer ::INTER(IGH), information on the rows and
    !!----    columns interchanged in the reduction.
    !!----
    Subroutine comhes ( n, low, igh, ar, ai, inter )
      integer,                         intent(in)     :: n,low, igh
      real(kind=dp), dimension(n,n),   intent(in out) :: ar,ai
      integer,       dimension(igh),   intent(out)    :: inter

      integer :: i,j,m
      real(kind=dp) :: t,xi,xr,yi,yr

      do m = low + 1, igh - 1
    !
    !  Choose the pivot I.
    !
        xr = 0.0_dp
        xi = 0.0_dp
        i = m

        do j = m, igh

          if ( abs ( xr )        + abs ( xi ) < &
               abs ( ar(j,m-1) ) + abs ( ai(j,m-1) ) ) then
            xr = ar(j,m-1)
            xi = ai(j,m-1)
            i = j
          end if

        end do

        inter(m) = i
    !
    !  Interchange rows and columns of AR and AI.
    !
        if ( i /= m ) then

          do j = m - 1, n
            t       = ar(i,j)
            ar(i,j) = ar(m,j)
            ar(m,j) = t
            t       = ai(i,j)
            ai(i,j) = ai(m,j)
            ai(m,j) = t
          end do

          do j = 1, igh
            t       = ar(j,i)
            ar(j,i) = ar(j,m)
            ar(j,m) = t
            t       = ai(j,i)
            ai(j,i) = ai(j,m)
            ai(j,m) = t
          end do

        end if
        !
        !  Carry out the transformation.
        !
        if ( xr /= 0.0_dp .or. xi /= 0.0_dp ) then

          do i = m + 1, igh

            yr = ar(i,m-1)
            yi = ai(i,m-1)

            if ( yr /= 0.0_dp .or. yi /= 0.0_dp ) then

              call cdiv ( yr, yi, xr, xi, yr, yi )
              ar(i,m-1) = yr
              ai(i,m-1) = yi

              do j = m, n
                ar(i,j) = ar(i,j) - yr * ar(m,j) + yi * ai(m,j)
                ai(i,j) = ai(i,j) - yr * ai(m,j) - yi * ar(m,j)
              end do

              ar(1:igh,m) = ar(1:igh,m) + yr * ar(1:igh,i) - yi * ai(1:igh,i)
              ai(1:igh,m) = ai(1:igh,m) + yr * ai(1:igh,i) + yi * ar(1:igh,i)

            end if

          end do

        end if

      end do

    End Subroutine comhes



    !!----  Subroutine comlr ( n, low, igh, hr, hi, wr, wi, ierr )
    !!----
    !!----  COMLR gets all eigenvalues of a complex upper Hessenberg matrix.
    !!----
    !!----  Discussion:
    !!----
    !!----    COMLR finds the eigenvalues of a complex upper Hessenberg matrix by the
    !!----    modified LR method.
    !!----
    !!----  Modified:
    !!----
    !!----    03 February 2018
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input, integer ::LOW, IGH, are determined by the balancing
    !!----    routine CBAL.  If CBAL is not used, set LOW = 1 and IGH = N.
    !!----
    !!----    Input/output, real(kind=dp) HR(N,N), HI(N,N).  On input, the real and
    !!----    imaginary parts of the complex upper Hessenberg matrix.  Their lower
    !!----    triangles below the subdiagonal contain the multipliers which were used
    !!----    in the reduction by COMHES if performed.  On output, the upper Hessenberg
    !!----    portions of HR and HI have been destroyed.  Therefore, they must be
    !!----    saved before calling COMLR if subsequent calculation of eigenvectors
    !!----    is to be performed.
    !!----
    !!----    Output, real(kind=dp) WR(N), WI(N), the real and imaginary parts of the
    !!----    eigenvalues.  If an error exit is made, the eigenvalues should be correct
    !!----    for indices IERR+1,...,N.
    !!----
    !!----    Output, integer ::IERR, error flag.
    !!----    0, for normal return,
    !!----    J, if the limit of 30*N iterations is exhausted while the J-th
    !!----      eigenvalue is being sought.
    !!----
    Subroutine comlr ( n, low, igh, hr, hi, wr, wi, ierr )
      integer,                         intent(in)     :: n,low, igh
      real(kind=dp), dimension(n,n),   intent(in out) :: hr,hi
      real(kind=dp), dimension(n),     intent(out)    :: wr,wi
      integer,                         intent(out)    :: ierr
      !
      integer :: en,i,itn,its,j,l,m
      real(kind=dp) :: ai,ar,si,sr,t,ti,tr,tst1,tst2,xi,xr,yi,yr,zzi,zzr

      ierr = 0
      !
      !  Store roots isolated by CBAL.
      !
      do i = 1, n
        if ( i < low .or. igh < i ) then
          wr(i) = hr(i,i)
          wi(i) = hi(i,i)
        end if
      end do

      en = igh
      tr = 0.0_dp
      ti = 0.0_dp
      itn = 30 * n
      !
      !  Search for next eigenvalue.
      !
      if ( en < low ) then
        return
      end if

      its = 0
      !
      !  Look for single small sub-diagonal element.
      !
      do

        do l = en, low, - 1

          if ( l == low ) then
            exit
          end if

          tst1 = abs ( hr(l-1,l-1) ) + abs ( hi(l-1,l-1) ) + abs ( hr(l,l) ) &
            + abs ( hi(l,l) )
          tst2 = tst1 + abs ( hr(l,l-1) ) + abs ( hi(l,l-1) )

          if ( tst2 == tst1 ) then
            exit
          end if

        end do
        !
        !  A root found.
        !
        if ( l == en ) then
          wr(en) = hr(en,en) + tr
          wi(en) = hi(en,en) + ti
          en = en - 1
          if ( en < low ) then
            return
          end if
          its = 0
          cycle
        end if

        if ( itn == 0 ) then
          ierr = en
          return
        end if

        if ( its == 10 .or. its == 20 ) then
          sr = abs ( hr(en,en-1) ) + abs ( hr(en-1,en-2) )
          si = abs ( hi(en,en-1) ) + abs ( hi(en-1,en-2) )
        else

          sr = hr(en,en)
          si = hi(en,en)
          xr = hr(en-1,en) * hr(en,en-1) - hi(en-1,en) * hi(en,en-1)
          xi = hr(en-1,en) * hi(en,en-1) + hi(en-1,en) * hr(en,en-1)

          if ( xr /= 0.0_dp .or. xi /= 0.0_dp ) then

            yr = ( hr(en-1,en-1) - sr) / 2.0_dp
            yi = ( hi(en-1,en-1) - si) / 2.0_dp
            ar = yr * yr - yi * yi + xr
            ai = 2.0_dp * yr * yi + xi
            call csroot ( ar, ai, zzr, zzi )

            if ( yr * zzr + yi * zzi < 0.0_dp ) then
              zzr = - zzr
              zzi = - zzi
            end if

            ar = yr + zzr
            ai = yi + zzi
            call cdiv ( xr, xi, ar, ai, xr, xi )
            sr = sr - xr
            si = si - xi

          end if

        end if

        do i = low, en
          hr(i,i) = hr(i,i) - sr
          hi(i,i) = hi(i,i) - si
        end do

        tr = tr + sr
        ti = ti + si
        its = its + 1
        itn = itn - 1
        !
        !  Look for two consecutive small sub-diagonal elements.
        !
        xr = abs ( hr(en-1,en-1) ) + abs ( hi(en-1,en-1) )
        yr = abs ( hr(en,en-1) ) + abs ( hi(en,en-1) )
        zzr = abs ( hr(en,en) ) + abs ( hi(en,en) )

        do m = en - 1, l, -1

          if ( m == l ) exit

          yi = yr
          yr = abs ( hr(m,m-1) ) + abs ( hi(m,m-1) )
          xi = zzr
          zzr = xr
          xr = abs ( hr(m-1,m-1) ) + abs ( hi(m-1,m-1) )
          tst1 = zzr / yi * ( zzr + xr + xi )
          tst2 = tst1 + yr
          if ( tst2 == tst1 ) exit
        end do
        !
        !  Triangular decomposition H=L*R.
        !
        do i = m + 1, en

          xr = hr(i-1,i-1)
          xi = hi(i-1,i-1)
          yr = hr(i,i-1)
          yi = hi(i,i-1)

          if (  abs ( yr ) + abs ( yi ) <= abs ( xr ) + abs ( xi ) ) then

            call cdiv ( yr, yi, xr, xi, zzr, zzi )
            wr(i) = - 1.0_dp

          else
            !
            !  Interchange rows of HR and HI.
            !
            do j = i - 1, en
              t         = hr(i-1,j)
              hr(i-1,j) = hr(i,j)
              hr(i,j)   = t
              t         = hi(i-1,j)
              hi(i-1,j) = hi(i,j)
              hi(i,j)   = t
            end do

            call cdiv ( xr, xi, yr, yi, zzr, zzi )
            wr(i) = 1.0_dp

          end if

          hr(i,i-1) = zzr
          hi(i,i-1) = zzi

          do j = i, en
            hr(i,j) = hr(i,j) - zzr * hr(i-1,j) + zzi * hi(i-1,j)
            hi(i,j) = hi(i,j) - zzr * hi(i-1,j) - zzi * hr(i-1,j)
          end do

        end do
        !
        !  Composition R*L=H.
        !
        do j = m + 1, en

          xr = hr(j,j-1)
          xi = hi(j,j-1)
          hr(j,j-1) = 0.0_dp
          hi(j,j-1) = 0.0_dp
          !
          !  Interchange columns of HR and HI, if necessary.
          !
          if ( 0.0_dp < wr(j) ) then

            do i = l, j
              t         = hr(i,j-1)
              hr(i,j-1) = hr(i,j)
              hr(i,j)   = t
              t         = hi(i,j-1)
              hi(i,j-1) = hi(i,j)
              hi(i,j)   = t
            end do

          end if

          do i = l, j
            hr(i,j-1) = hr(i,j-1) + xr * hr(i,j) - xi * hi(i,j)
            hi(i,j-1) = hi(i,j-1) + xr * hi(i,j) + xi * hr(i,j)
          end do

        end do

      end do

    End Subroutine comlr



    !!----  Subroutine comlr2 ( n, low, igh, inter, hr, hi, wr, wi, zr, zi, ierr )
    !!----
    !!----  COMLR2 gets eigenvalues/vectors of a complex upper Hessenberg matrix.
    !!----
    !!----  Discussion:
    !!----
    !!----    COMLR2 finds the eigenvalues and eigenvectors of a complex
    !!----    upper Hessenberg matrix by the modified LR method.
    !!----
    !!----    The eigenvectors of a complex general matrix can also be found if
    !!----    COMHES has been used to reduce this general matrix to Hessenberg form.
    !!----
    !!----  Modified:
    !!----
    !!----    10 February 2018
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input, integer ::LOW, IGH, are determined by the balancing
    !!----    routine CBAL.  If CBAL is not used, set LOW = 1 and IGH = N.
    !!----
    !!----    Input, integer ::INTER(IGH), information on the rows and columns
    !!----    interchanged in the reduction by COMHES, if performed.  If the
    !!----    eigenvectors of the Hessenberg matrix are desired, set INTER(J)=J for these
    !!----    elements.
    !!----
    !!----    Input/output, real(kind=dp) HR(N,N), HI(N,N).  On input, the real
    !!----    and imaginary parts of the complex upper Hessenberg matrix.  Their lower
    !!----    triangles below the subdiagonal contain the multipliers which were used in
    !!----    the reduction by COMHES, if performed.  If the eigenvectors of the
    !!----    Hessenberg matrix are desired, these elements must be set to zero.  On
    !!----    output, the upper Hessenberg portions of HR and HI have been destroyed,
    !!----    but the location HR(1,1) contains the norm of the triangularized matrix.
    !!----
    !!----    Output, real(kind=dp) WR(N), WI(N), the real and imaginary parts of the
    !!----    eigenvalues.  If an error exit is made, the eigenvalues should be
    !!----    correct for indices IERR+1,...,N.
    !!----
    !!----    Output, real(kind=dp) ZR(N,N), ZI(N,N), the real and imaginary parts
    !!----    of the eigenvectors.  The eigenvectors are unnormalized.  If an error exit
    !!----    is made, none of the eigenvectors has been found.
    !!----
    !!----    Output, integer ::IERR, error flag.
    !!----    0, for normal return,
    !!----    J, if the limit of 30*N iterations is exhausted while the J-th
    !!----      eigenvalue is being sought.
    !!----
    Subroutine comlr2 ( n, low, igh, inter, hr, hi, wr, wi, zr, zi, ierr )
      integer,                         intent(in)     :: n,low, igh
      integer,       dimension(igh),   intent(in)     :: inter
      real(kind=dp), dimension(n,n),   intent(in out) :: hr,hi
      real(kind=dp), dimension(n),     intent(out)    :: wr,wi
      real(kind=dp), dimension(n,n),   intent(out)    :: zr,zi
      integer,                         intent(out)    :: ierr

      integer       :: en,enm1,i,itn,its,j,k,l,m
      real(kind=dp) :: norm,si,sr,t,ti,tr,tst1,tst2,xi,xr,yi,yr,zzi,zzr

      ierr = 0
      !
      !  Initialize the eigenvector matrix.
      !
      call r8mat_identity ( n, zr )

      zi(1:n,1:n) = 0.0_dp
      !
      !  Form the matrix of accumulated transformations from the information left
      !  by COMHES.
      !
      do i = igh - 1, low + 1, -1

        do k = i + 1, igh
          zr(k,i) = hr(k,i-1)
          zi(k,i) = hi(k,i-1)
        end do

        j = inter(i)

        if ( i /= j ) then

          do k = i, igh
            zr(i,k) = zr(j,k)
            zi(i,k) = zi(j,k)
            zr(j,k) = 0.0_dp
            zi(j,k) = 0.0_dp
          end do

          zr(j,i) = 1.0_dp

        end if

      end do
      !
      !  Store roots isolated by CBAL.
      !
      do i = 1, n
        if ( i < low .or. igh < i ) then
          wr(i) = hr(i,i)
          wi(i) = hi(i,i)
        end if
      end do

      en = igh
      tr = 0.0_dp
      ti = 0.0_dp
      itn = 30 * n
       !
       !  Search for next eigenvalue.
       !
      its = 0
      enm1 = en - 1
      !
      !  Look for single small sub-diagonal element.
      !
      do

        do l = en, low, -1

          if ( l == low ) then
            exit
          end if

          tst1 = abs ( hr(l-1,l-1) ) + abs ( hi(l-1,l-1) ) + abs ( hr(l,l) ) &
            + abs ( hi(l,l) )
          tst2 = tst1 + abs ( hr(l,l-1) ) + abs ( hi(l,l-1) )

          if ( tst2 == tst1 ) then
            exit
          end if

        end do
        !
        !  A root found.
        !
        if ( l == en ) then
          hr(en,en) = hr(en,en) + tr
          wr(en) = hr(en,en)
          hi(en,en) = hi(en,en) + ti
          wi(en) = hi(en,en)
          en = enm1
          if ( en < low ) then
            exit
          end if
          its = 0
          enm1 = en - 1
          cycle
        end if

        if ( itn == 0 ) then
          ierr = en
          return
        end if
        !
        !  Form shift.
        !
        if ( its == 10 .or. its == 20 ) then

          sr = abs ( hr(en,enm1) ) + abs ( hr(enm1,en-2) )
          si = abs ( hi(en,enm1) ) + abs ( hi(enm1,en-2) )

        else

          sr = hr(en,en)
          si = hi(en,en)
          xr = hr(enm1,en) * hr(en,enm1) - hi(enm1,en) * hi(en,enm1)
          xi = hr(enm1,en) * hi(en,enm1) + hi(enm1,en) * hr(en,enm1)

          if ( xr /= 0.0_dp .or. xi /= 0.0_dp ) then

            yr = ( hr(enm1,enm1) - sr ) / 2.0_dp
            yi = ( hi(enm1,enm1) - si ) / 2.0_dp
            call csroot ( yr**2 - yi**2 + xr, 2.0_dp * yr * yi + xi, zzr, zzi )

            if ( yr * zzr + yi * zzi < 0.0_dp ) then
              zzr = - zzr
              zzi = - zzi
            end if

            call cdiv ( xr, xi, yr + zzr, yi + zzi, xr, xi )
            sr = sr - xr
            si = si - xi

          end if

        end if

        do i = low, en
          hr(i,i) = hr(i,i) - sr
          hi(i,i) = hi(i,i) - si
        end do

        tr = tr + sr
        ti = ti + si
        its = its + 1
        itn = itn - 1
        !
        !  Look for two consecutive small sub-diagonal elements.
        !
        xr = abs ( hr(enm1,enm1) ) + abs ( hi(enm1,enm1) )
        yr = abs ( hr(en,enm1) ) + abs ( hi(en,enm1) )
        zzr = abs ( hr(en,en) ) + abs ( hi(en,en) )

        do m = enm1, l, -1
          if ( m == l ) then
            exit
          end if
          yi = yr
          yr = abs ( hr(m,m-1) ) + abs ( hi(m,m-1) )
          xi = zzr
          zzr = xr
          xr = abs ( hr(m-1,m-1) ) + abs ( hi(m-1,m-1) )
          tst1 = zzr / yi * ( zzr + xr + xi )
          tst2 = tst1 + yr
          if ( tst2 == tst1 ) then
            exit
          end if
        end do
        !
        !  Triangular decomposition H=L*R.
        !
        do i = m + 1, en

          xr = hr(i-1,i-1)
          xi = hi(i-1,i-1)
          yr = hr(i,i-1)
          yi = hi(i,i-1)
          !
          !  Interchange rows of HR and HI.
          !
          if ( abs ( xr ) + abs ( xi) < abs ( yr ) + abs ( yi ) ) then

            do j = i - 1, n
              t         = hr(i-1,j)
              hr(i-1,j) = hr(i,j)
              hr(i,j)   = t
              t         = hi(i-1,j)
              hi(i-1,j) = hi(i,j)
              hi(i,j)   = t
            end do

            call cdiv ( xr, xi, yr, yi, zzr, zzi )
            wr(i) = 1.0_dp

          else

            call cdiv ( yr, yi, xr, xi, zzr, zzi )
            wr(i) = - 1.0_dp

          end if

          hr(i,i-1) = zzr
          hi(i,i-1) = zzi

          do j = i, n
            hr(i,j) = hr(i,j) - zzr * hr(i-1,j) + zzi * hi(i-1,j)
            hi(i,j) = hi(i,j) - zzr * hi(i-1,j) - zzi * hr(i-1,j)
          end do

        end do
        !
        !  Composition R*L=H.
        !
        do j = m + 1, en

          xr = hr(j,j-1)
          xi = hi(j,j-1)
          hr(j,j-1) = 0.0_dp
          hi(j,j-1) = 0.0_dp
          !
          !  Interchange columns of HR, HI, ZR, and ZI.
          !
          if ( 0.0_dp < wr(j) ) then

            do i = 1, j
              t         = hr(i,j-1)
              hr(i,j-1) = hr(i,j)
              hr(i,j)   = t
              t         = hi(i,j-1)
              hi(i,j-1) = hi(i,j)
              hi(i,j)   = t
            end do

            do i = low, igh
              t         = zr(i,j-1)
              zr(i,j-1) = zr(i,j)
              zr(i,j)   = t
              t         = zi(i,j-1)
              zi(i,j-1) = zi(i,j)
              zi(i,j)   = t
            end do

          end if

          do i = 1, j
            hr(i,j-1) = hr(i,j-1) + xr * hr(i,j) - xi * hi(i,j)
            hi(i,j-1) = hi(i,j-1) + xr * hi(i,j) + xi * hr(i,j)
          end do
          !
          !  Accumulate transformations.
          !
          do i = low, igh
            zr(i,j-1) = zr(i,j-1) + xr * zr(i,j) - xi * zi(i,j)
            zi(i,j-1) = zi(i,j-1) + xr * zi(i,j) + xi * zr(i,j)
          end do

        end do

      end do
      !
      !  All roots found.
      !  Backsubstitute to find vectors of upper triangular form.
      !
      norm = 0.0_dp
      do i = 1, n
        do j = i, n
          tr = abs ( hr(i,j) ) + abs ( hi(i,j) )
          norm = max ( norm, tr )
        end do
      end do

      hr(1,1) = norm
      if ( n == 1 ) then
        return
      end if

      if ( norm == 0.0_dp ) then
        return
      end if

      do en = n, 2, -1

        xr = wr(en)
        xi = wi(en)
        hr(en,en) = 1.0_dp
        hi(en,en) = 0.0_dp
        enm1 = en - 1

        do i = en - 1, 1, -1

          zzr = 0.0_dp
          zzi = 0.0_dp

          do j = i + 1, en
            zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
            zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
          end do

          yr = xr - wr(i)
          yi = xi - wi(i)

          if ( yr == 0.0_dp .and. yi == 0.0_dp ) then

            tst1 = norm
            yr = tst1

            do
              yr = 0.01_dp * yr
              tst2 = norm + yr
              if ( tst2 <=  tst1 ) then
                exit
              end if
            end do

          end if

          call cdiv ( zzr, zzi, yr, yi, hr(i,en), hi(i,en) )
          !
          !  Overflow control.
          !
          tr = abs ( hr(i,en) ) + abs ( hi(i,en) )

          if ( tr /= 0.0_dp ) then

            tst1 = tr
            tst2 = tst1 + 1.0_dp / tst1

            if ( tst2 <= tst1 ) then
              hr(i:en,en) = hr(i:en,en) / tr
              hi(i:en,en) = hi(i:en,en) / tr
            end if

          end if

        end do

      end do
      !
      !  End backsubstitution.
      !
      enm1 = n - 1
      !
      !  Vectors of isolated roots.
      !
      do i = 1, n - 1

        if ( i < low .or. igh < i ) then

          zr(i,i+1:n) = hr(i,i+1:n)
          zi(i,i+1:n) = hi(i,i+1:n)

        end if

      end do
      !
      !  Multiply by transformation matrix to give vectors of original full matrix.
      !
      do j = n, low + 1, -1

        m = min ( j, igh )

        do i = low, igh
          zzr = 0.0_dp
          zzi = 0.0_dp
          do k = low, m
            zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
            zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
          end do
          zr(i,j) = zzr
          zi(i,j) = zzi
        end do

      end do

    End Subroutine comlr2



    !!----  Subroutine comqr ( n, low, igh, hr, hi, wr, wi, ierr )
    !!----
    !!----  COMQR gets eigenvalues of a complex upper Hessenberg matrix.
    !!----
    !!----  Discussion:
    !!----
    !!----    COMQR finds the eigenvalues of a complex upper Hessenberg matrix by
    !!----    the QR method.
    !!----
    !!----  Modified:
    !!----
    !!----    10 February 2018
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input, integer ::LOW, IGH, are determined by the balancing
    !!----    routine CBAL.  If CBAL is not used, set LOW = 1 and IGH = N.
    !!----
    !!----    Input/output, real(kind=dp) HR(N,N), HI(N,N).  On input, the real
    !!----    and imaginary parts of the complex upper Hessenberg matrix.  Their lower
    !!----    triangles below the subdiagonal contain information about the unitary
    !!----    transformations used in the reduction by CORTH, if performed.  On output,
    !!----    the upper Hessenberg portions of HR and HI have been destroyed.
    !!----    Therefore, they must be saved before calling COMQR if subsequent
    !!----    calculation of eigenvectors is to be performed.
    !!----
    !!----    Output, real(kind=dp) WR(N), WI(N), the real and imaginary parts of the
    !!----    eigenvalues.  If an error exit is made, the eigenvalues should be
    !!----    correct for indices IERR+1,...,N.
    !!----
    !!----    Output, integer ::IERR, error flag.
    !!----    0, for normal return,
    !!----    J, if the limit of 30*N iterations is exhausted while the J-th
    !!----       eigenvalue is being sought.
    !!----
    Subroutine comqr ( n, low, igh, hr, hi, wr, wi, ierr )
      integer,                         intent(in)     :: n,low, igh
      real(kind=dp), dimension(n,n),   intent(in out) :: hr,hi
      real(kind=dp), dimension(n),     intent(out)    :: wr,wi
      integer,                         intent(out)    :: ierr

      integer :: en,enm1,i,itn,its,j,l,ll
      real(kind=dp) :: norm,si,sr,ti,tr,tst1,tst2,xi,xr,yi,yr,zzi,zzr

      ierr = 0
    !
    !  Create real subdiagonal elements.
    !
      l = low + 1

      do i = low + 1, igh

        ll = min ( i + 1, igh )

        if ( hi(i,i-1) /= 0.0_dp ) then

          norm = pythag ( hr(i,i-1), hi(i,i-1) )
          yr = hr(i,i-1) / norm
          yi = hi(i,i-1) / norm
          hr(i,i-1) = norm
          hi(i,i-1) = 0.0_dp

          do j = i, igh
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
          end do

          do j = low, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
          end do

        end if

      end do
      !
      !  Store roots isolated by CBAL.
      !
      do i = 1, n
        if ( i < low .or. igh < i ) then
          wr(i) = hr(i,i)
          wi(i) = hi(i,i)
        end if
      end do

      en = igh
      tr = 0.0_dp
      ti = 0.0_dp
      itn = 30 * n
      !
      !  Search for next eigenvalue.
      !
      if ( en < low ) return

      its = 0
      enm1 = en - 1
      !
      !  Look for single small sub-diagonal element.
      !
      do

        do l = en, low, -1
          if ( l == low ) then
            exit
          end if
          tst1 = abs ( hr(l-1,l-1) ) + abs ( hi(l-1,l-1) ) + abs ( hr(l,l) ) &
            + abs ( hi(l,l) )
          tst2 = tst1 + abs ( hr(l,l-1) )
          if ( tst2 == tst1 ) then
            exit
          end if
        end do
        !
        !  A root found.
        !
        if ( l == en ) then
          wr(en) = hr(en,en) + tr
          wi(en) = hi(en,en) + ti
          en = enm1
          if ( en < low ) then
            exit
          end if
          its = 0
          enm1 = en - 1
          cycle
        end if

        if ( itn == 0 ) then
          ierr = en
          exit
        end if

        if ( its == 10 .or. its == 20 ) then

          sr = abs ( hr(en,enm1) ) + abs ( hr(enm1,en-2) )
          si = 0.0_dp

        else

          sr = hr(en,en)
          si = hi(en,en)
          xr = hr(enm1,en) * hr(en,enm1)
          xi = hi(enm1,en) * hr(en,enm1)

          if ( xr /= 0.0_dp .or. xi /= 0.0_dp ) then

            yr = ( hr(enm1,enm1) - sr ) / 2.0_dp
            yi = ( hi(enm1,enm1) - si ) / 2.0_dp

            call csroot ( yr**2-yi**2+xr, 2.0_dp*yr*yi+xi, zzr, zzi )

            if ( yr * zzr + yi * zzi < 0.0_dp ) then
              zzr = - zzr
              zzi = - zzi
            end if

            call cdiv ( xr, xi, yr+zzr, yi+zzi, xr, xi )
            sr = sr - xr
            si = si - xi

          end if

        end if

        do i = low, en
          hr(i,i) = hr(i,i) - sr
          hi(i,i) = hi(i,i) - si
        end do

        tr = tr + sr
        ti = ti + si
        its = its + 1
        itn = itn - 1
        !
        !  Reduce to triangle (rows).
        !
        do i = l + 1, en

          sr = hr(i,i-1)
          hr(i,i-1) = 0.0_dp
          norm = pythag ( pythag ( hr(i-1,i-1), hi(i-1,i-1) ), sr )
          xr = hr(i-1,i-1) / norm
          wr(i-1) = xr
          xi = hi(i-1,i-1) / norm
          wi(i-1) = xi
          hr(i-1,i-1) = norm
          hi(i-1,i-1) = 0.0_dp
          hi(i,i-1) = sr / norm

          do j = i, en
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
          end do

        end do

        si = hi(en,en)

        if ( si /= 0.0_dp ) then
          norm = pythag ( hr(en,en), si )
          sr = hr(en,en) / norm
          si = si / norm
          hr(en,en) = norm
          hi(en,en) = 0.0_dp
        end if
        !
        !  Inverse operation (columns).
        !
        do j = l + 1, en

          xr = wr(j-1)
          xi = wi(j-1)

          do i = l, j

            yr = hr(i,j-1)
            yi = 0.0_dp
            zzr = hr(i,j)
            zzi = hi(i,j)
            if ( i /= j ) then
              yi = hi(i,j-1)
              hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
            end if
            hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi

          end do

        end do

        if ( si /= 0.0_dp ) then

          do i = l, en
            yr = hr(i,en)
            yi = hi(i,en)
            hr(i,en) = sr * yr - si * yi
            hi(i,en) = sr * yi + si * yr
          end do

        end if

      end do

    End Subroutine comqr

    !!----  Subroutine comqr2 ( n, low, igh, ortr, orti, hr, hi, wr, wi, zr, zi, ierr )
    !!----
    !!----  COMQR2 gets eigenvalues/vectors of a complex upper Hessenberg matrix.
    !!----
    !!----  Discussion:
    !!----
    !!----    COMQR2 finds the eigenvalues and eigenvectors
    !!----    of a complex upper Hessenberg matrix by the QR
    !!----    method.  The eigenvectors of a complex general matrix
    !!----    can also be found if CORTH has been used to reduce
    !!----    this general matrix to Hessenberg form.
    !!----
    !!----  Modified:
    !!----
    !!----    04 March 2018
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input, integer ::LOW, IGH, are determined by the balancing
    !!----    routine CBAL.  If CBAL is not used, set LOW = 1 and IGH = N.
    !!----
    !!----    Input/output, real(kind=dp) ORTR(N), ORTI(N).  On input, information
    !!----    about the unitary transformations used in the reduction by CORTH, if
    !!----    performed.  If the eigenvectors of the Hessenberg matrix are desired, set
    !!----    ORTR(J) and ORTI(J) to 0.0_dp for these elements.  On output, these arrays
    !!----    have been overwritten.
    !!----
    !!----    Input/output, real(kind=dp) HR(N,N), HI(N,N).  On input, the real and
    !!----    imaginary parts of the complex upper Hessenberg matrix.  Their lower
    !!----    triangles below the subdiagonal contain further information about the
    !!----    transformations which were used in the reduction by CORTH, if performed.
    !!----    If the eigenvectors of the Hessenberg matrix are desired, these elements
    !!----    may be arbitrary.
    !!----
    !!----    Output, real(kind=dp) WR(N), WI(N), the real and imaginary parts of the
    !!----    eigenvalues.  If an error exit is made, the eigenvalues should be
    !!----    correct for indices IERR+1,...,N.
    !!----
    !!----    Output, real(kind=dp) ZR(N,N), ZI(N,N), the real and imaginary parts of
    !!----    the eigenvectors.  The eigenvectors are unnormalized.  If an error exit
    !!----    is made, none of the eigenvectors has been found.
    !!----
    !!----    Output, integer ::IERR, error flag.
    !!----    0, for normal return,
    !!----    J, if the limit of 30*N iterations is exhausted while the J-th
    !!----      eigenvalue is being sought.
    !!----
    Subroutine comqr2 ( n, low, igh, ortr, orti, hr, hi, wr, wi, zr, zi, ierr )
      integer,                         intent(in)     :: n,low, igh
      real(kind=dp), dimension(igh),   intent(in out) :: ortr, orti
      real(kind=dp), dimension(n,n),   intent(in out) :: hr,hi
      real(kind=dp), dimension(n),     intent(out)    :: wr,wi
      real(kind=dp), dimension(n,n),   intent(out)    :: zr,zi
      integer,                         intent(out)    :: ierr

      integer :: en,enm1,i,itn,its,j,k,l,ll,m
      real(kind=dp) :: norm,si,sr,ti,tr,tst1,tst2,xi,xr,yi,yr,zzi,zzr

      ierr = 0
      !
      !  Initialize eigenvector matrix.
      !
      call r8mat_identity ( n, zr )

      zi(1:n,1:n) = 0.0_dp
      !
      !  Form the matrix of accumulated transformations from the information
      !  left by CORTH.
      !
      if ( 0 < igh - low - 1 ) then

        do i = igh - 1, low + 1, -1

          if ( ortr(i) == 0.0_dp .and. orti(i) == 0.0_dp ) cycle

          if ( hr(i,i-1) == 0.0_dp .and. hi(i,i-1) == 0.0_dp ) cycle
          !
          !  Norm below is negative of H formed in CORTH.
          !
          norm = hr(i,i-1) * ortr(i) + hi(i,i-1) * orti(i)

          do k = i + 1, igh
            ortr(k) = hr(k,i-1)
            orti(k) = hi(k,i-1)
          end do

          do j = i, igh

            sr = 0.0_dp
            si = 0.0_dp
            do k = i, igh
              sr = sr + ortr(k) * zr(k,j) + orti(k) * zi(k,j)
              si = si + ortr(k) * zi(k,j) - orti(k) * zr(k,j)
            end do
            sr = sr / norm
            si = si / norm

            do k = i, igh
              zr(k,j) = zr(k,j) + sr * ortr(k) - si * orti(k)
              zi(k,j) = zi(k,j) + sr * orti(k) + si * ortr(k)
            end do

          end do

        end do

      end if
      !
      !  Create real subdiagonal elements.
      !
      if ( 0 <= igh - low - 1 ) then

        l = low + 1

        do i = low + 1, igh

          ll = min ( i + 1, igh )

          if ( hi(i,i-1) == 0.0_dp ) cycle

          norm = pythag ( hr(i,i-1), hi(i,i-1) )
          yr = hr(i,i-1) / norm
          yi = hi(i,i-1) / norm
          hr(i,i-1) = norm
          hi(i,i-1) = 0.0_dp

          do j = i, n
            si =      yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
          end do

          do j = 1, ll
            si =      yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
          end do

          do j = low, igh
            si =      yr * zi(j,i) + yi * zr(j,i)
            zr(j,i) = yr * zr(j,i) - yi * zi(j,i)
            zi(j,i) = si
          end do

        end do

      end if
      !
      !  Store roots isolated by CBAL.
      !
      do i = 1, n
        if ( i < low .or. igh < i ) then
          wr(i) = hr(i,i)
          wi(i) = hi(i,i)
        end if
      end do

      en = igh
      tr = 0.0_dp
      ti = 0.0_dp
      itn = 30 * n
      !
      !  Search for next eigenvalue.
      !
      its = 0
      enm1 = en - 1
      !
      !  Look for single small sub-diagonal element.
      !
      do

        do l = en, low + 1, -1
          tst1 = abs ( hr(l-1,l-1) ) + abs ( hi(l-1,l-1) ) + abs ( hr(l,l) ) &
            + abs ( hi(l,l) )
          tst2 = tst1 + abs ( hr(l,l-1) )
          if ( tst2 == tst1 ) then
            exit
          end if
        end do
        !
        !  A root found.
        !
        if ( l == en ) then
          hr(en,en) = hr(en,en) + tr
          wr(en) = hr(en,en)
          hi(en,en) = hi(en,en) + ti
          wi(en) = hi(en,en)
          en = enm1
          if ( en < low ) exit
          its = 0
          enm1 = en - 1
          cycle
        end if

        if ( itn == 0 ) then
          ierr = en
          return
        end if

        if ( its == 10 .or. its == 20 ) then

          sr = abs ( hr(en,enm1) ) + abs ( hr(enm1,en-2) )
          si = 0.0_dp

        else

          sr = hr(en,en)
          si = hi(en,en)
          xr = hr(enm1,en) * hr(en,enm1)
          xi = hi(enm1,en) * hr(en,enm1)

          if ( xr /= 0.0_dp .or. xi /= 0.0_dp ) then

            yr = ( hr(enm1,enm1) - sr ) / 2.0_dp
            yi = ( hi(enm1,enm1) - si ) / 2.0_dp

            call csroot ( yr**2 - yi**2 + xr, 2.0_dp * yr * yi + xi, zzr, zzi )

            if ( yr * zzr + yi * zzi < 0.0_dp ) then
              zzr = - zzr
              zzi = - zzi
            end if

            call cdiv ( xr, xi, yr + zzr, yi + zzi, xr, xi )
            sr = sr - xr
            si = si - xi

          end if

        end if

        do i = low, en
          hr(i,i) = hr(i,i) - sr
          hi(i,i) = hi(i,i) - si
        end do

        tr = tr + sr
        ti = ti + si
        its = its + 1
        itn = itn - 1
        !
        !  Reduce to triangle (rows).
        !
        do i = l + 1, en

          sr = hr(i,i-1)
          hr(i,i-1) = 0.0_dp
          norm = pythag ( pythag ( hr(i-1,i-1), hi(i-1,i-1) ), sr )
          xr = hr(i-1,i-1) / norm
          wr(i-1) = xr
          xi = hi(i-1,i-1) / norm
          wi(i-1) = xi
          hr(i-1,i-1) = norm
          hi(i-1,i-1) = 0.0_dp
          hi(i,i-1) = sr / norm

          do j = i, n
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
          end do

        end do

        si = hi(en,en)

        if ( si /= 0.0_dp ) then

          norm = pythag ( hr(en,en), si )
          sr = hr(en,en) / norm
          si = si / norm
          hr(en,en) = norm
          hi(en,en) = 0.0_dp

          do j = en + 1, n
            yr = hr(en,j)
            yi = hi(en,j)
            hr(en,j) = sr * yr + si * yi
            hi(en,j) = sr * yi - si * yr
          end do

        end if
        !
        !  Inverse operation (columns).
        !
        do j = l + 1, en

          xr = wr(j-1)
          xi = wi(j-1)

          do i = 1, j

            yr = hr(i,j-1)
            yi = 0.0_dp
            zzr = hr(i,j)
            zzi = hi(i,j)

            if ( i /= j ) then
              yi = hi(i,j-1)
              hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
            end if

            hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi

          end do

          do i = low, igh
            yr = zr(i,j-1)
            yi = zi(i,j-1)
            zzr = zr(i,j)
            zzi = zi(i,j)
            zr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            zi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
            zr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            zi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
          end do

        end do

        if ( si /= 0.0_dp ) then

          do i = 1, en
            yr = hr(i,en)
            yi = hi(i,en)
            hr(i,en) = sr * yr - si * yi
            hi(i,en) = sr * yi + si * yr
          end do

          do i = low, igh
            yr = zr(i,en)
            yi = zi(i,en)
            zr(i,en) = sr * yr - si * yi
            zi(i,en) = sr * yi + si * yr
          end do

        end if

      end do
      !
      !  All roots found.
      !  Backsubstitute to find vectors of upper triangular form.
      !
      norm = 0.0_dp

      do i = 1, n
        do j = i, n
          tr = abs ( hr(i,j) ) + abs ( hi(i,j) )
          norm = max ( norm, tr )
        end do
      end do

      if ( n == 1 ) return

      if ( norm == 0.0_dp ) return

      do en = n, 2, -1

        xr = wr(en)
        xi = wi(en)
        hr(en,en) = 1.0_dp
        hi(en,en) = 0.0_dp
        enm1 = en - 1

        do i = en - 1, 1, -1

          zzr = 0.0_dp
          zzi = 0.0_dp

          do j = i + 1, en
            zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
            zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
          end do

          yr = xr - wr(i)
          yi = xi - wi(i)

          if ( yr == 0.0_dp .and. yi == 0.0_dp ) then

            tst1 = norm
            yr = tst1

            do
              yr = 0.01_dp * yr
              tst2 = norm + yr
              if ( tst2 <= tst1 ) exit
            end do

          end if

          call cdiv ( zzr, zzi, yr, yi, hr(i,en), hi(i,en) )
          !
          !  Overflow control.
          !
          tr = abs ( hr(i,en) ) + abs ( hi(i,en) )

          if ( tr /= 0.0_dp ) then

            tst1 = tr
            tst2 = tst1 + 1.0_dp / tst1

            if ( tst2 <= tst1 ) then

              do j = i, en
                hr(j,en) = hr(j,en) / tr
                hi(j,en) = hi(j,en) / tr
              end do

            end if

          end if

        end do

      end do
      !
      !  End backsubstitution.
      !
      enm1 = n - 1
      !
      !  Vectors of isolated roots.
      !
      do i = 1, n - 1

        if ( i < low .or. igh < i ) then

          do j = i + 1, n
            zr(i,j) = hr(i,j)
            zi(i,j) = hi(i,j)
          end do

        end if

      end do
      !
      !  Multiply by transformation matrix to give vectors of original full matrix.
      !
      do j = n, low + 1, -1

        m = min ( j, igh )

        do i = low, igh

          zzr = 0.0_dp
          zzi = 0.0_dp
          do k = low, m
            zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
            zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
          end do

          zr(i,j) = zzr
          zi(i,j) = zzi

        end do

      end do

    End Subroutine comqr2



    !!----  Subroutine cortb ( n, low, igh, ar, ai, ortr, orti, m, zr, zi )
    !!----
    !!----  CORTB determines eigenvectors by undoing the CORTH transformation.
    !!----
    !!----  Discussion:
    !!----
    !!----    CORTB forms the eigenvectors of a complex general matrix by back
    !!----    transforming those of the corresponding upper Hessenberg matrix
    !!----    determined by CORTH.
    !!----
    !!----  Modified:
    !!----
    !!----    24 February 2018
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input, integer ::LOW, IGH, are determined by the balancing
    !!----    routine CBAL.  If CBAL is not used, set LOW = 1 and IGH to the order
    !!----    of the matrix.
    !!----
    !!----    Input, real(kind=dp) AR(N,IGH), AI(N,IGH), information about the
    !!----    unitary transformations used in the reduction by CORTH in their strict
    !!----    lower triangles.
    !!----
    !!----    Input/output, real(kind=dp) ORTR(IGH), ORTI(IGH).  On input, further
    !!----    information about the transformations used in the reduction by CORTH.  On
    !!----    output, ORTR and ORTI have been further altered.
    !!----
    !!----    Input, integer ::M, the number of columns of ZR and ZI to be
    !!----    back transformed.
    !!----
    !!----    Input/output, real(kind=dp) ZR(N,M), ZI(N,M).  On input, the real and
    !!----    imaginary parts of the eigenvectors to be back transformed.  On output,
    !!----    the real and imaginary parts of the transformed eigenvectors.
    !!----
    Subroutine cortb ( n, low, igh, ar, ai, ortr, orti, m, zr, zi )
      integer,                         intent(in)     :: n,low, igh
      real(kind=dp), dimension(n,n),   intent(in)     :: ar,ai
      real(kind=dp), dimension(igh),   intent(in out) :: ortr, orti
      integer,                         intent(in)     :: m
      real(kind=dp), dimension(n,m),   intent(in out) :: zr,zi

      integer :: i,j,mp
      real(kind=dp) :: gi, gr, h

      if ( m == 0 ) return

      if ( igh - 1 < low + 1 ) return

      do mp = igh - 1, low + 1, -1

        if ( ar(mp,mp-1) /= 0.0_dp .or. ai(mp,mp-1) /= 0.0_dp ) then

          h = ar(mp,mp-1) * ortr(mp) + ai(mp,mp-1) * orti(mp)

          ortr(mp+1:igh) = ar(mp+1:igh,mp-1)
          orti(mp+1:igh) = ai(mp+1:igh,mp-1)

          do j = 1, m

            gr = ( dot_product ( ortr(mp:igh), zr(mp:igh,j) ) &
                 + dot_product ( orti(mp:igh), zi(mp:igh,j) ) ) / h

            gi = ( dot_product ( ortr(mp:igh), zi(mp:igh,j) ) &
                 - dot_product ( orti(mp:igh), zr(mp:igh,j) ) ) / h

            do i = mp, igh
              zr(i,j) = zr(i,j) + gr * ortr(i) - gi * orti(i)
              zi(i,j) = zi(i,j) + gr * orti(i) + gi * ortr(i)
            end do

          end do

        end if

      end do

    End Subroutine cortb



    !!----  Subroutine corth ( n, low, igh, ar, ai, ortr, orti )
    !!----
    !!----  CORTH transforms a complex general matrix to upper Hessenberg form.
    !!----
    !!----  Discussion:
    !!----
    !!----    CORTH is given a complex general matrix and reduces a submatrix situated
    !!----    in rows and columns LOW through IGH to upper Hessenberg form by
    !!----    unitary similarity transformations.
    !!----
    !!----  Modified:
    !!----
    !!----    04 March 2018
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input, integer ::LOW, IGH, are determined by the balancing
    !!----    routine CBAL.  If CBAL is not used, set LOW = 1 and IGH = N.
    !!----
    !!----    Input/output, real(kind=dp) AR(N,N), AI(N,N).  On input, the real and
    !!----    imaginary parts of the complex input matrix.  On output, the real and
    !!----    imaginary parts of the Hessenberg matrix.  Information about the unitary
    !!----    transformations used in the reduction is stored in the remaining
    !!----    triangles under the Hessenberg matrix.
    !!----
    !!----    Output, real(kind=dp) ORTR(IGH), ORTI(IGH), further information about
    !!----    the transformations.
    !!----
    Subroutine corth ( n, low, igh, ar, ai, ortr, orti )
      integer,                         intent(in)     :: n,low, igh
      real(kind=dp), dimension(n,n),   intent(in out) :: ar,ai
      real(kind=dp), dimension(igh),   intent(out)    :: ortr, orti

      real(kind=dp) :: f,fi,fr,g,h,scal
      integer :: i,j,m

      if ( igh - 1 < low + 1 ) return

      do m = low + 1, igh - 1

        h = 0.0_dp
        ortr(m) = 0.0_dp
        orti(m) = 0.0_dp
        scal = 0.0_dp
        !
        !  scal column.
        !
        do i = m, igh
          scal = scal + abs ( ar(i,m-1) ) + abs ( ai(i,m-1) )
        end do

        if ( scal == 0.0_dp ) cycle

        do i = igh, m, -1
          ortr(i) = ar(i,m-1) / scal
          orti(i) = ai(i,m-1) / scal
          h = h + ortr(i) * ortr(i) + orti(i) * orti(i)
        end do

        g = sqrt ( h )
        f = pythag ( ortr(m), orti(m) )

        if ( f /= 0.0_dp ) then
          h = h + f * g
          g = g / f
          ortr(m) = ( 1.0_dp + g ) * ortr(m)
          orti(m) = ( 1.0_dp + g ) * orti(m)
        else
          ortr(m) = g
          ar(m,m-1) = scal
        end if
        !
        !  Form (I-(U*Ut)/h) * A.
        !
        do j = m, n

          fr = 0.0_dp
          fi = 0.0_dp
          do i = igh, m, -1
            fr = fr + ortr(i) * ar(i,j) + orti(i) * ai(i,j)
            fi = fi + ortr(i) * ai(i,j) - orti(i) * ar(i,j)
          end do
          fr = fr / h
          fi = fi / h

          ar(m:igh,j) = ar(m:igh,j) - fr * ortr(m:igh) + fi * orti(m:igh)
          ai(m:igh,j) = ai(m:igh,j) - fr * orti(m:igh) - fi * ortr(m:igh)

        end do
        !
        !  Form (I-(U*Ut)/h) * A * (I-(U*Ut)/h)
        !
        do i = 1, igh

          fr = 0.0_dp
          fi = 0.0_dp
          do j = igh, m, - 1
            fr = fr + ortr(j) * ar(i,j) - orti(j) * ai(i,j)
            fi = fi + ortr(j) * ai(i,j) + orti(j) * ar(i,j)
          end do
          fr = fr / h
          fi = fi / h

          ar(i,m:igh) = ar(i,m:igh) - fr * ortr(m:igh) - fi * orti(m:igh)
          ai(i,m:igh) = ai(i,m:igh) + fr * orti(m:igh) - fi * ortr(m:igh)

        end do

        ortr(m) = scal * ortr(m)
        orti(m) = scal * orti(m)
        ar(m,m-1) = - g * ar(m,m-1)
        ai(m,m-1) = - g * ai(m,m-1)

      end do

    End Subroutine corth



    !!----  Subroutine csroot ( xr, xi, yr, yi )
    !!----
    !!----  CSROOT computes the complex square root of a complex quantity.
    !!----
    !!----  Discussion:
    !!----
    !!----    CSROOT chooses the branch of the square function so that
    !!----      0.0 <= YR
    !!----    and
    !!----      sign ( YI ) == sign ( XI )
    !!----
    !!----  Modified:
    !!----
    !!----    18 October 2009
    !!----
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, real(kind=dp) XR, XI, the real and imaginary parts of the
    !!----    quantity whose square root is desired.
    !!----
    !!----    Output, real(kind=dp) YR, YI, the real and imaginary parts of the
    !!----    square root.
    !!----
    Subroutine csroot ( xr, xi, yr, yi )
      real(kind=dp), intent(in)  :: xr, xi
      real(kind=dp), intent(out) :: yr, yi

      real(kind=dp) :: s,ti,tr

      tr = xr
      ti = xi
      s = sqrt ( 0.5_dp * ( pythag ( tr, ti ) + abs ( tr ) ) )

      if ( 0.0_dp <= tr ) then
        yr = s
      end if

      if ( ti < 0.0_dp ) then
        s = -s
      end if

      if ( tr <= 0.0_dp ) then
        yi = s
      end if

      if ( tr < 0.0_dp ) then
        yr = 0.5_dp * ( ti / yi )
      else if ( 0.0_dp < tr ) then
        yi = 0.5_dp * ( ti / yr )
      end if

    End Subroutine csroot



    !!----  Subroutine elmbak ( n, low, igh, a, ind, m, z )
    !!----
    !!----  ELMBAK determines eigenvectors by undoing the ELMHES transformation.
    !!----
    !!----  Discussion:
    !!----
    !!----    ELMBAK forms the eigenvectors of a real general matrix by back
    !!----    transforming those of the corresponding upper Hessenberg matrix
    !!----    determined by ELMHES.
    !!----
    !!----  Modified:
    !!----
    !!----    04 March 2018
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input, integer ::LOW, IGH, integers determined by the balancing
    !!----    routine BALANC.  If BALANC has not been used, set LOW = 1 and
    !!----    IGH equal to the order of the matrix.
    !!----
    !!----    Input, real(kind=dp) A(N,IGH), the multipliers which were used in the
    !!----    reduction by ELMHES in its lower triangle below the subdiagonal.
    !!----
    !!----    Input, integer ::IND(IGH), information on the rows and columns
    !!----    interchanged in the reduction by ELMHES.
    !!----
    !!----    Input, integer ::M, the number of columns of Z to be back
    !!----    transformed.
    !!----
    !!----    Input/output, real(kind=dp) Z(N,M).  On input, the real and imaginary
    !!----    parts of the eigenvectors to be back transformed.  On output, the real and
    !!----    imaginary parts of the transformed eigenvectors.
    !!----
    Subroutine elmbak ( n, low, igh, a, ind, m, z )
      integer,                         intent(in)     :: n,low, igh
      real(kind=dp), dimension(n,igh), intent(in)     :: a
      integer,       dimension(igh),   intent(in)     :: ind
      integer,                         intent(in)     :: m
      real(kind=dp), dimension(n,m),   intent(in out) :: z

      integer :: i,j,mm,mp
      real(kind=dp) :: t,x


      if ( m == 0 ) return

      if ( igh - 1 < low + 1 ) return

      do mm = low + 1, igh - 1

        mp = low + igh - mm

        do i = mp + 1, igh

          x = a(i,mp-1)
          if ( x /= 0.0_dp ) then
            z(i,1:m) = z(i,1:m) + x * z(mp,1:m)
          end if

        end do

        i = ind(mp)

        if ( i /= mp ) then

          do j = 1, m
            t       = z(i,j)
            z(i,j)  = z(mp,j)
            z(mp,j) = t
          end do

        end if

      end do

    End Subroutine elmbak


    !!----  Subroutine elmhes ( n, low, igh, a, ind )
    !!----
    !!----  ELMHES transforms a real general matrix to upper Hessenberg form.
    !!----
    !!----  Discussion:
    !!----
    !!----    ELMHES is given a real general matrix and reduces a submatrix
    !!----    situated in rows and columns LOW through IGH to upper Hessenberg
    !!----    form by stabilized elementary similarity transformations.
    !!----
    !!----
    !!----  Modified:
    !!----
    !!----    11 March 2018
    !!----
    !!----
    !!----  Reference:
    !!----
    !!----    Martin, James Wilkinson,
    !!----    ELMHES,
    !!----    Numerische Mathematik,
    !!----    Volume 12, pages 349-368, 1968.
    !!----
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input, integer ::LOW, IGH, are determined by the balancing
    !!----    routine BALANC.  If BALANC has not been used, set LOW = 1, IGH = N.
    !!----
    !!----    Input/output, real(kind=dp) A(N,N).  On input, the matrix to be
    !!----    reduced.  On output, the Hessenberg matrix.  The multipliers
    !!----    which were used in the reduction are stored in the
    !!----    remaining triangle under the Hessenberg matrix.
    !!----
    !!----    Output, integer ::IND(N), contains information on the rows
    !!----    and columns interchanged in the reduction.  Only elements LOW through
    !!----    IGH are used.
    !!----
    Subroutine elmhes ( n, low, igh, a, ind )
      integer,                         intent(in)     :: n,low, igh
      real(kind=dp), dimension(n,n),   intent(in out) :: a
      integer,       dimension(igh),   intent(out)    :: ind

      integer :: i,j,m
      real(kind=dp) :: t,x,y

      ind(1:igh) = 0

      do m = low + 1, igh - 1

        x = 0.0_dp
        i = m

        do j = m, igh
          if ( abs ( x ) < abs ( a(j,m-1) ) ) then
            x = a(j,m-1)
            i = j
          end if
        end do

        ind(m) = i
        !
        !  Interchange rows and columns of the matrix.
        !
        if ( i /= m ) then

          do j = m - 1, n
            t      = a(i,j)
            a(i,j) = a(m,j)
            a(m,j) = t
          end do

          do j = 1, igh
            t      = a(j,i)
            a(j,i) = a(j,m)
            a(j,m) = t

          end do

        end if

        if ( x /= 0.0_dp ) then

          do i = m + 1, igh

            y = a(i,m-1)

            if ( y /= 0.0_dp ) then

              y = y / x
              a(i,m-1) = y

              a(i,m:n) = a(i,m:n) - y * a(m,m:n)

              a(1:igh,m) = a(1:igh,m) + y * a(1:igh,i)

            end if

          end do

        end if

      end do

    End Subroutine elmhes



    !!----  Subroutine eltran ( n, low, igh, a, ind, z )
    !!----
    !!----  ELTRAN accumulates similarity transformations used by ELMHES.
    !!----
    !!----  Discussion:
    !!----
    !!----    ELTRAN accumulates the stabilized elementary similarity transformations
    !!----    used in the reduction of a real general matrix to upper Hessenberg form
    !!----    by ELMHES.
    !!----
    !!----  Modified:
    !!----
    !!----    04 March 2018
    !!----
    !!----
    !!----  Reference:
    !!----
    !!----    Peters, James Wilkinson,
    !!----    ELMTRANS,
    !!----    Numerische Mathematik,
    !!----    Volume 16, pages 181-204, 1970.
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input, integer ::LOW, IGH, are determined by the balancing
    !!----    routine BALANC.  If BALANC has not been used, set LOW = 1, IGH = N.
    !!----
    !!----    Input, real(kind=dp) A(N,IGH), the multipliers which were used in the
    !!----    reduction by ELMHES in its lower triangle below the subdiagonal.
    !!----
    !!----    Input, integer ::IND(IGH), information on the rows and columns
    !!----    interchanged in the reduction by ELMHES.
    !!----
    !!----    Output, real(kind=dp) Z(N,N), the transformation matrix produced in the
    !!----    reduction by ELMHES.
    !!----
    Subroutine eltran ( n, low, igh, a, ind, z )
      integer,                         intent(in)     :: n,low, igh
      real(kind=dp), dimension(n,igh), intent(in)     :: a
      integer,       dimension(igh),   intent(in)     :: ind
      real(kind=dp), dimension(n,n),   intent(out)    :: z


      integer :: i,kl,mm,mp
      !
      !  Initialize Z to the identity matrix.
      !
      call r8mat_identity ( n, z )

      kl = igh - low - 1

      if ( igh - 1 < low + 1 )return

      do mm = 1, igh - low - 1

         mp = igh - mm

         z(mp+1:igh,mp) = a(mp+1:igh,mp-1)

         i = ind(mp)

         if ( i /= mp ) then

           z(mp,mp:igh) = z(i,mp:igh)

           z(i,mp) = 1.0_dp
           z(i,mp+1:igh) = 0.0_dp

         end if

      end do

    End Subroutine eltran



    !!----  Subroutine figi ( n, t, d, e, e2, ierr )
    !!----
    !!----  FIGI transforms a real nonsymmetric tridiagonal matrix to symmetric form.
    !!----
    !!----  Discussion:
    !!----
    !!----    FIGI is given a nonsymmetric tridiagonal matrix such that the products
    !!----    of corresponding pairs of off-diagonal elements are all
    !!----    non-negative.  It reduces the matrix to a symmetric
    !!----    tridiagonal matrix with the same eigenvalues.  If, further,
    !!----    a zero product only occurs when both factors are zero,
    !!----    the reduced matrix is similar to the original matrix.
    !!----
    !!----  Modified:
    !!----
    !!----    08 February 2018
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input, real(kind=dp) T(N,3) contains the input matrix.  Its subdiagonal
    !!----    is stored in the last N-1 positions of the first column, its diagonal in
    !!----    the N positions of the second column, and its superdiagonal in the
    !!----    first N-1 positions of the third column.  T(1,1) and T(N,3) are arbitrary.
    !!----
    !!----    Output, real(kind=dp) D(N), the diagonal elements of the symmetric
    !!----    matrix.
    !!----
    !!----    Output, real(kind=dp) E(N), contains the subdiagonal elements of
    !!----    the symmetric matrix in E(2:N).  E(1) is not set.
    !!----
    !!----    Output, real(kind=dp) E2(N), the squares of the corresponding elements
    !!----    of E.  E2 may coincide with E if the squares are not needed.
    !!----
    !!----    Output, integer ::IERR, error flag.
    !!----    0, for normal return,
    !!----    N+I, if T(I,1) * T(I-1,3) is negative,
    !!----    -(3*N+I), if T(I,1) * T(I-1,3) is zero with one factor non-zero.  In
    !!----      this case, the eigenvectors of the symmetric matrix are not simply
    !!----      related to those of T and should not be sought.
    !!----
    Subroutine figi ( n, t, d, e, e2, ierr )
      integer,                      intent(in)  :: n
      real(kind=dp), dimension(n,3),intent(in)  :: t
      real(kind=dp), dimension(n),  intent(out) :: d,e,e2
      integer,                      intent(out) :: ierr

      integer :: i
      ierr = 0

      do i = 1, n

        if ( 1 < i ) then

          e2(i) = t(i,1) * t(i-1,3)

          if ( e2(i) < 0.0_dp ) then

            ierr = n + i
            return

          else if ( e2(i) == 0.0_dp ) then

            if ( t(i,1) /= 0.0_dp .or. t(i-1,3) /= 0.0_dp ) then
              ierr = - 3 * n - i
              return
            end if

            e(i) = 0.0_dp

          else

            e(i) = sqrt ( e2(i) )

          end if

        else

          e(i) = 0.0_dp
          e2(i) = 0.0_dp

        end if

        d(i) = t(i,2)

      end do

    End Subroutine figi



    !!----  Subroutine figi2 ( n, t, d, e, z, ierr )
    !!----
    !!----  FIGI2 transforms a real nonsymmetric tridiagonal matrix to symmetric form.
    !!----
    !!----  Discussion:
    !!----
    !!----    FIGI2 is given a nonsymmetric tridiagonal matrix such that the products
    !!----    of corresponding pairs of off-diagonal elements are all
    !!----    non-negative, and zero only when both factors are zero.
    !!----
    !!----    It reduces the matrix to a symmetric tridiagonal matrix
    !!----    using and accumulating diagonal similarity transformations.
    !!----
    !!----  Modified:
    !!----
    !!----    08 February 2018
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input, real(kind=dp) T(N,3) contains the input matrix.  Its subdiagonal
    !!----    is stored in the last N-1 positions of the first column, its diagonal in
    !!----    the N positions of the second column, and its superdiagonal in the
    !!----    first N-1 positions of the third column.  T(1,1) and T(N,3) are arbitrary.
    !!----
    !!----    Output, real(kind=dp) D(N), the diagonal elements of the symmetric
    !!----    matrix.
    !!----
    !!----    Output, real(kind=dp) E(N), contains the subdiagonal elements of the
    !!----    symmetric matrix in E(2:N).  E(1) is not set.
    !!----
    !!----    Output, real(kind=dp) Z(N,N), contains the transformation matrix
    !!----    produced in the reduction.
    !!----
    !!----    Output, integer ::IERR, error flag.
    !!----    0, for normal return,
    !!----    N+I, if T(I,1) * T(I-1,3) is negative,
    !!----    2*N+I, if T(I,1) * T(I-1,3) is zero with one factor non-zero.
    !!----
    Subroutine figi2 ( n, t, d, e, z, ierr )
      integer,                      intent(in)  :: n
      real(kind=dp), dimension(n,3),intent(in)  :: t
      real(kind=dp), dimension(n),  intent(out) :: d,e
      real(kind=dp), dimension(n,n),intent(out) :: z
      integer,                      intent(out) :: ierr


      real(kind=dp) :: h
      integer :: i

      ierr = 0
    !
    !  Initialize Z to the identity matrix.
    !
      call r8mat_identity ( n, z )

      do i = 1, n

        if ( i == 1 ) then

          e(i) = 0.0_dp

        else

          h = t(i,1) * t(i-1,3)

          if ( h < 0.0_dp ) then

            ierr = n + i
            return

          else if ( h == 0 ) then

            if ( t(i,1) /= 0.0_dp .or. t(i-1,3) /= 0.0_dp ) then
              ierr = 2 * n + i
              return
            end if

            e(i) = 0.0_dp
            z(i,i) = 1.0_dp

          else

            e(i) = sqrt ( h )
            z(i,i) = z(i-1,i-1) * e(i) / t(i-1,3)

          end if

        end if

        d(i) = t(i,2)

      end do

    End Subroutine figi2



    !!----  Subroutine hqr ( n, low, igh, h, wr, wi, ierr )
    !!----
    !!----  HQR computes all eigenvalues of a real upper Hessenberg matrix.
    !!----
    !!----  Discussion:
    !!----
    !!----    HQR finds the eigenvalues of a real
    !!----    upper Hessenberg matrix by the QR method.
    !!----
    !!----  Modified:
    !!----
    !!----    31 January 2018
    !!----
    !!----  Reference:
    !!----
    !!----    Martin, Peters, James Wilkinson,
    !!----    HQR,
    !!----    Numerische Mathematik,
    !!----    Volume 14, pages 219-231, 1970.
    !!----
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input, integer ::LOW, IGH, two integers determined by
    !!----    BALANC.  If BALANC is not used, set LOW=1, IGH=N.
    !!----
    !!----    Input/output, real(kind=dp) H(N,N), the N by N upper Hessenberg matrix.
    !!----    Information about the transformations used in the reduction to
    !!----    Hessenberg form by ELMHES or ORTHES, if performed, is stored
    !!----    in the remaining triangle under the Hessenberg matrix.
    !!----    On output, the information in H has been destroyed.
    !!----
    !!----    Output, real(kind=dp) WR(N), WI(N), the real and imaginary parts of the
    !!----    eigenvalues.  The eigenvalues are unordered, except that complex
    !!----    conjugate pairs of values appear consecutively, with the eigenvalue
    !!----    having positive imaginary part listed first.  If an error exit
    !!----    occurred, then the eigenvalues should be correct for indices
    !!----    IERR+1 through N.
    !!----
    !!----    Output, integer ::IERR, error flag.
    !!----    0, no error.
    !!----    J, the limit of 30*N iterations was reached while searching for
    !!----      the J-th eigenvalue.
    !!----
    Subroutine hqr ( n, low, igh, h, wr, wi, ierr )
      integer,                      intent(in)     :: n, low, igh
      real(kind=dp), dimension(n,n),intent(in out) :: h
      real(kind=dp), dimension(n),  intent(out)    :: wr, wi
      integer,                      intent(out)    :: ierr

      integer :: en,enm2,i,itn,its,j,k,l,m,na
      logical :: notlas
      real(kind=dp) :: norm,p,q,r,s,t,tst1,tst2,w,x,y,zz

      ierr = 0
      norm = 0.0_dp
      k = 1
      !
      !  Store roots isolated by BALANC and compute matrix norm.
      !
      do i = 1, n

        do j = k, n
          norm = norm + abs ( h(i,j) )
        end do

        k = i
        if ( i < low .or. igh < i ) then
          wr(i) = h(i,i)
          wi(i) = 0.0_dp
        end if

      end do

      en = igh
      t = 0.0_dp
      itn = 30 * n
      !
      !  Search for next eigenvalues.
      !
      if ( igh < low ) return

      its = 0
      na = igh - 1
      enm2 = igh - 2
      !
      !  Look for a single small sub-diagonal element.
      !
      do

        do l = en, low, -1
          if ( l == low ) exit
          s = abs ( h(l-1,l-1) ) + abs ( h(l,l) )
          if ( s == 0.0_dp ) then
            s = norm
          end if
          tst1 = s
          tst2 = tst1 + abs ( h(l,l-1) )
          if ( tst2 == tst1 ) exit
        end do
        !
        !  Form shift.
        !
        x = h(en,en)
        !
        !  One root found.
        !
        if ( l == en ) then
          wr(en) = x + t
          wi(en) = 0.0_dp
          en = na
          if ( en < low ) return
          its = 0
          na = en - 1
          enm2 = na - 1
          cycle
        end if

        y = h(na,na)
        w = h(en,na) * h(na,en)
        !
        !  Two roots found.
        !
        if ( l == na ) then

          p = ( y - x ) / 2.0_dp
          q = p * p + w
          zz = sqrt ( abs ( q ) )
          x = x + t
          !
          !  Real root, or complex pair.
          !
          if ( 0.0_dp <= q ) then

            zz = p + sign ( zz, p )
            wr(na) = x + zz
            if ( zz == 0.0_dp ) then
              wr(en) = wr(na)
            else
              wr(en) = x - w / zz
            end if
            wi(na) = 0.0_dp
            wi(en) = 0.0_dp

          else

            wr(na) = x + p
            wr(en) = x + p
            wi(na) = zz
            wi(en) = - zz

          end if

          en = enm2

          if ( en < low ) then
            return
          end if

          its = 0
          na = en - 1
          enm2 = na - 1
          cycle

        end if

        if ( itn == 0 ) then
          ierr = en
          return
        end if
        !
        !  Form an exceptional shift.
        !
        if ( its == 10 .or. its == 20 ) then

          t = t + x

          do i = low, en
            h(i,i) = h(i,i) - x
          end do

          s = abs ( h(en,na) ) + abs ( h(na,enm2) )
          x = 0.75_dp * s
          y = x
          w = - 0.4375_dp * s * s

        end if

        its = its + 1
        itn = itn - 1
        !
        !  Look for two consecutive small sub-diagonal elements.
        !
        do m = enm2, l, -1

          zz = h(m,m)
          r = x - zz
          s = y - zz
          p = ( r * s - w ) / h(m+1,m) + h(m,m+1)
          q = h(m+1,m+1) - zz - r - s
          r = h(m+2,m+1)
          s = abs ( p ) + abs ( q ) + abs ( r )
          p = p / s
          q = q / s
          r = r / s

          if ( m == l ) exit

          tst1 = abs ( p ) * ( abs ( h(m-1,m-1) ) + abs ( zz ) &
            + abs ( h(m+1,m+1) ) )
          tst2 = tst1 + abs ( h(m,m-1) ) * ( abs ( q ) + abs ( r ) )

          if ( tst2 == tst1 ) exit

        end do

        do i = m + 2, en
          h(i,i-2) = 0.0_dp
          if ( i /= m + 2 ) then
            h(i,i-3) = 0.0_dp
          end if
        end do
        !
        !  Double QR step involving rows l to EN and columns M to EN.
        !
        do k = m, na

          notlas = k /= na

          if ( k /= m ) then

            p = h(k,k-1)
            q = h(k+1,k-1)

            if ( notlas ) then
              r = h(k+2,k-1)
            else
              r = 0.0_dp
            end if

            x = abs ( p ) + abs ( q ) + abs ( r )

            if ( x == 0.0_dp ) cycle

            p = p / x
            q = q / x
            r = r / x

          end if

          s = sign ( sqrt ( p * p + q * q + r * r ), p )

          if ( k /= m ) then
            h(k,k-1) = - s * x
          else if ( l /= m ) then
            h(k,k-1) = - h(k,k-1)
          end if

          p = p + s
          x = p / s
          y = q / s
          zz = r / s
          q = q / p
          r = r / p

          if ( .not. notlas ) then
            !
            !  Row modification.
            !
            do j = k, n
              p = h(k,j) + q * h(k+1,j)
              h(k,j) = h(k,j) - p * x
              h(k+1,j) = h(k+1,j) - p * y
            end do

            j = min ( en, k + 3 )
            !
            !  Column modification.
            !
            do i = 1, j
              p = x * h(i,k) + y * h(i,k+1)
              h(i,k) = h(i,k) - p
              h(i,k+1) = h(i,k+1) - p * q
            end do

          else
            !
            !  Row modification.
            !
            do j = k, n
              p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
              h(k,j) = h(k,j) - p * x
              h(k+1,j) = h(k+1,j) - p * y
              h(k+2,j) = h(k+2,j) - p * zz
            end do

            j = min ( en, k + 3 )
            !
            !  Column modification.
            !
            do i = 1, j
              p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
              h(i,k) = h(i,k) - p
              h(i,k+1) = h(i,k+1) - p * q
              h(i,k+2) = h(i,k+2) - p * r
            end do

          end if

        end do

      end do

    End Subroutine hqr



    !!----  Subroutine hqr2 ( n, low, igh, h, wr, wi, z, ierr )
    !!----
    !!----  HQR2 computes eigenvalues and eigenvectors of a real upper Hessenberg matrix.
    !!----
    !!----  Discussion:
    !!----
    !!----    HQR2 finds the eigenvalues and eigenvectors of a real upper Hessenberg
    !!----    matrix by the qr method.
    !!----
    !!----    The eigenvectors of a real general matrix can also be found
    !!----    if ELMHES and ELTRAN or ORTHES and ORTRAN have
    !!----    been used to reduce this general matrix to Hessenberg form
    !!----    and to accumulate the similarity transformations.
    !!----
    !!----   THE STATEMENT BELOW CAME FROM THE ORIGINAL Eispack.f90 code
    !!----    Thanks to David Chichka, 02 May 2019, for pointing out that a previous
    !!----    version of this F90 translation of hqr2 gave erroneous results when some
    !!----    eigenvalues were complex.  I gave up my ideal of a complete rewrite of
    !!----    this function, and did a much lighter conversion of the F77 code
    !!----    for a second try.
    !!----   IN THIS VERSION ALL GOTO's AND NUMERICAL LABELS HAVE BEEN ELIMINATED
    !!----
    !!----  Modified:
    !!----
    !!----    02 May 2019
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input, integer ::LOW, IGH, determined by the balancing
    !!----    routine BALANC.  If BALANC has not been used, set LOW = 1, IGH = N.
    !!----
    !!----    Input/output, real(kind=dp) H(N,N), the N by N upper Hessenberg matrix.
    !!----    On output, the information in H has been destroyed.
    !!----
    !!----    Output, real(kind=dp) WR(N), WI(N), the real and imaginary parts of the
    !!----    eigenvalues.  The eigenvalues are unordered, except that complex
    !!----    conjugate pairs of values appear consecutively, with the eigenvalue
    !!----    having positive imaginary part listed first.  If an error exit
    !!----    occurred, then the eigenvalues should be correct for indices
    !!----    IERR+1 through N.
    !!----
    !!----    Input/output, real(kind=dp) Z(N,N).  On input, the transformation
    !!----    matrix produced by ELTRAN after the reduction by ELMHES, or by ORTRAN after
    !!----    the reduction by ORTHES, if performed.  If the eigenvectors of the
    !!----    Hessenberg matrix are desired, Z must contain the identity matrix.  On
    !!----    output, Z contains the real and imaginary parts of the eigenvectors.
    !!----    If the I-th eigenvalue is real, the I-th column of Z contains its
    !!----    eigenvector.  If the I-th eigenvalue is complex with positive imaginary
    !!----    part, the I-th and (I+1)-th columns of Z contain the real and imaginary
    !!----    parts of its eigenvector.  The eigenvectors are unnormalized.  If an
    !!----    error exit is made, none of the eigenvectors has been found.
    !!----
    !!----    Output, integer ::IERR, error flag.
    !!----    0, for normal return,
    !!----    J, if the limit of 30*N iterations is exhausted while the J-th
    !!----      eigenvalue is being sought.
    !!----

    Subroutine hqr2 ( n, low, igh, h, wr, wi, z, ierr )
      integer,                      intent(in)     :: n, low, igh
      real(kind=dp), dimension(n,n),intent(in out) :: h
      real(kind=dp), dimension(n),  intent(out)    :: wr, wi
      real(kind=dp), dimension(n,n),intent(in out) :: z
      integer,                      intent(out)    :: ierr

      integer :: i,j,k,l,m,en,ii,jj,ll,mm,na,nn
      integer :: itn,its,mp2,enm2
      real(kind=dp) :: p,q,r,s,t,w,x,y,ra,sa,vi,vr,zz,norm,tst1,tst2
      logical :: notlas,two_found !This last logical has been introduced to simplify the code (JRC)

       ierr = 0
       norm = 0.0_dp
       k = 1
       two_found=.false.
       !
       !  store roots isolated by balanc and compute matrix norm
       !
       do i = 1, n

          do j = k, n
            norm = norm + abs(h(i,j))
          end do

          k = i
          if ( i >= low .and. i <= igh ) cycle
          wr(i) = h(i,i)
          wi(i) = 0.0_dp
       end do

       en = igh
       t = 0.0_dp
       itn = 30*n
       !
       !  search for next eigenvalues
       !
       do

          if (en < low) exit
          its = 0
          na = en - 1
          enm2 = na - 1
         !
         !  look for single small sub-diagonal element
         !  for l=en step -1 until low do --
         !
         do

            do ll = low, en
               l = en + low - ll
               if (l == low) exit
               s = abs(h(l-1,l-1)) + abs(h(l,l))
               if (s == 0.0_dp) s = norm
               tst1 = s
               tst2 = tst1 + abs(h(l,l-1))
               if (tst2 == tst1) exit
            end do
            !
            !  form shift
            !
            x = h(en,en)
            if (l == en) then
              two_found=.false.
              exit
            end if
            y = h(na,na)
            w = h(en,na) * h(na,en)
            if (l == na) then
              two_found=.true.
              exit
            end if
            if (itn == 0)  then
              !
              !  set error -- all eigenvalues have not converged after 30*n iterations
              !
              ierr = en
              return
            end if
            if (.not. (its /= 10 .and. its /= 20)) then
              !
              !  form exceptional shift
              !
               t = t + x

               do i = low, en
                 h(i,i) = h(i,i) - x
               end do

               s = abs(h(en,na)) + abs(h(na,enm2))
               x = 0.75_dp * s
               y = x
               w = -0.4375_dp * s * s
            end if
            its = its + 1
            itn = itn - 1
            !
            !  look for two consecutive small sub-diagonal elements.
            !  for m=en-2 step -1 until l do --
            !
            do mm = l, enm2
               m = enm2 + l - mm
               zz = h(m,m)
               r = x - zz
               s = y - zz
               p = (r * s - w) / h(m+1,m) + h(m,m+1)
               q = h(m+1,m+1) - zz - r - s
               r = h(m+2,m+1)
               s = abs(p) + abs(q) + abs(r)
               p = p / s
               q = q / s
               r = r / s
               if (m == l) exit
               tst1 = abs(p)*(abs(h(m-1,m-1)) + abs(zz) + abs(h(m+1,m+1)))
               tst2 = tst1 + abs(h(m,m-1))*(abs(q) + abs(r))
               if (tst2 == tst1) exit
            end do

            mp2 = m + 2

            do i = mp2, en
               h(i,i-2) = 0.0_dp
               if (i == mp2) cycle
               h(i,i-3) = 0.0_dp
            end do
            !
            !  double qr step involving rows l to en and columns m to en.
            !
            do k = m, na

               notlas = k /= na
               if (k /= m) then
                 p = h(k,k-1)
                 q = h(k+1,k-1)
                 r = 0.0_dp
                 if (notlas) r = h(k+2,k-1)
                 x = abs(p) + abs(q) + abs(r)
                 if (x == 0.0_dp) cycle
                 p = p / x
                 q = q / x
                 r = r / x
               end if

               s = sign(sqrt(p*p+q*q+r*r),p)
               if (k == m) then
                 if (l /= m) h(k,k-1) = -h(k,k-1)
               else
                 h(k,k-1) = -s * x
               end if
               p = p + s
               x = p / s
               y = q / s
               zz = r / s
               q = q / p
               r = r / p
               if (notlas) then
                !
                !  row modification
                !
                 do j = k, n
                    p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
                    h(k,j) = h(k,j) - p * x
                    h(k+1,j) = h(k+1,j) - p * y
                    h(k+2,j) = h(k+2,j) - p * zz
                 end do

                 j = min(en,k+3)
                 !
                 !  column modification
                 !
                 do i = 1, j
                    p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
                    h(i,k) = h(i,k) - p
                    h(i,k+1) = h(i,k+1) - p * q
                    h(i,k+2) = h(i,k+2) - p * r
                 end do
                 !
                 !  accumulate transformations
                 !
                 do i = low, igh
                    p = x * z(i,k) + y * z(i,k+1) + zz * z(i,k+2)
                    z(i,k) = z(i,k) - p
                    z(i,k+1) = z(i,k+1) - p * q
                    z(i,k+2) = z(i,k+2) - p * r
                 end do

               else
                 !
                 !  row modification.
                 !
                 do j = k, n
                    p = h(k,j) + q * h(k+1,j)
                    h(k,j) = h(k,j) - p * x
                    h(k+1,j) = h(k+1,j) - p * y
                 end do

                 j = min(en,k+3)
                 !
                 !  column modification
                 !
                 do i = 1, j
                    p = x * h(i,k) + y * h(i,k+1)
                    h(i,k) = h(i,k) - p
                    h(i,k+1) = h(i,k+1) - p * q
                 end do
                 !
                 !  accumulate transformations
                 !
                 do i = low, igh
                    p = x * z(i,k) + y * z(i,k+1)
                    z(i,k) = z(i,k) - p
                    z(i,k+1) = z(i,k+1) - p * q
                 end do
              end if

            end do

         end do

         if(.not. two_found) then
            !
            !  one root found
            !
            h(en,en) = x + t
            wr(en) = h(en,en)
            wi(en) = 0.0_dp
            en = na
            cycle
         else
            !
            !  two roots found
            !
            p = (y - x) / 2.0_dp
            q = p * p + w
            zz = sqrt(abs(q))
            h(en,en) = x + t
            x = h(en,en)
            h(na,na) = y + t
         end if

         if (q >= 0.0_dp) then
             !
             !  real pair
             !
             zz = p + sign(zz,p)
             wr(na) = x + zz
             wr(en) = wr(na)
             if (zz /= 0.0_dp) wr(en) = x - w / zz
             wi(na) = 0.0_dp
             wi(en) = 0.0_dp
             x = h(en,na)
             s = abs(x) + abs(zz)
             p = x / s
             q = zz / s
             r = sqrt(p*p+q*q)
             p = p / r
             q = q / r
             !
             !  row modification
             !
             do j = na, n
                zz = h(na,j)
                h(na,j) = q * zz + p * h(en,j)
                h(en,j) = q * h(en,j) - p * zz
             end do
             !
             !  column modification
             !
             do i = 1, en
                zz = h(i,na)
                h(i,na) = q * zz + p * h(i,en)
                h(i,en) = q * h(i,en) - p * zz
             end do
             !
             !  accumulate transformations
             !
             do i = low, igh
                zz = z(i,na)
                z(i,na) = q * zz + p * z(i,en)
                z(i,en) = q * z(i,en) - p * zz
             end do
         else
             !
             !  complex pair
             !
             wr(na) = x + p
             wr(en) = x + p
             wi(na) = zz
             wi(en) = -zz
         end if

         en = enm2
       end do
       !
       !  all roots found.  backsubstitute to find vectors of upper triangular form
       !
       if (norm == 0.0_dp) return
       !
       !  for en=n step -1 until 1 do
       !
          do nn = 1, n  !do 800

             en = n + 1 - nn
             p = wr(en)
             q = wi(en)
             na = en - 1
             if (q > 0.0_dp) cycle
             if (q < 0.0_dp) then
                !
                !  complex vector
                !
                m = na
                !
                !  last vector component chosen imaginary so eigenvector matrix is triangular
                !
                if (abs(h(en,na)) > abs(h(na,en))) then
                  h(na,na) = q / h(en,na)
                  h(na,en) = -(h(en,en) - p) / h(en,na)
                else
                  call cdiv(0.0_dp,-h(na,en),h(na,na)-p,q,h(na,na),h(na,en))
                end if

                h(en,na) = 0.0_dp
                h(en,en) = 1.0_dp
                enm2 = na - 1
                if (enm2 == 0) cycle
                !
                !  for i=en-2 step -1 until 1 do
                !
                do ii = 1, enm2

                   i = na - ii
                   w = h(i,i) - p
                   ra = 0.0_dp
                   sa = 0.0_dp

                   do j = m, en
                      ra = ra + h(i,j) * h(j,na)
                      sa = sa + h(i,j) * h(j,en)
                   end do

                   if (wi(i) < 0.0_dp) then
                      zz = w
                      r = ra
                      s = sa
                      cycle
                   end if

                   m = i
                   if (wi(i) == 0.0_dp) then
                     call cdiv(-ra,-sa,w,q,h(i,na),h(i,en))
                   else
                     !
                     !  solve complex equations
                     !
                     x = h(i,i+1)
                     y = h(i+1,i)
                     vr = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i) - q * q
                     vi = (wr(i) - p) * 2.0_dp * q
                     if (.not. (vr /= 0.0_dp .or. vi /= 0.0_dp)) then
                        tst1 = norm * (abs(w) + abs(q) + abs(x) + abs(y) + abs(zz))
                        vr = tst1
                     else
                        do
                          vr = 0.01d0 * vr
                          tst2 = tst1 + vr
                          if (tst2 > tst1) cycle
                          exit
                        end do
                     end if

                     call cdiv(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi,h(i,na),h(i,en))
                     if (abs(x) > abs(zz) + abs(q)) then
                       h(i+1,na) = (-ra - w * h(i,na) + q * h(i,en)) / x
                       h(i+1,en) = (-sa - w * h(i,en) - q * h(i,na)) / x
                     else
                       call cdiv(-r-y*h(i,na),-s-y*h(i,en),zz,q,h(i+1,na),h(i+1,en))
                     end if
                   end if
                   !
                   !  overflow control
                   !
                   t = max(abs(h(i,na)), abs(h(i,en)))
                   if (t == 0.0_dp) cycle
                   tst1 = t
                   tst2 = tst1 + 1.0_dp/tst1
                   if (tst2 > tst1) cycle
                   do j = i, en
                      h(j,na) = h(j,na)/t
                      h(j,en) = h(j,en)/t
                   end do

                end do
                !
                !  end complex vector
                !
             else if (q == 0.0_dp) then
               !
               !  real vector
               !
               m = en
               h(en,en) = 1.0_dp
               if (na == 0) cycle
               !
               !  for i=en-1 step -1 until 1 do
               !
               do ii = 1, na

                  i = en - ii
                  w = h(i,i) - p
                  r = 0.0_dp

                  do j = m, en
                    r = r + h(i,j) * h(j,en)
                  end do

                  if (wi(i) < 0.0_dp) then
                     zz = w
                     s = r
                     cycle
                  end if

                  m = i
                  if (wi(i) == 0.0_dp) then
                    t = w
                    if (t == 0.0_dp) then
                       tst1 = norm
                       t = tst1
                       do
                         t = 0.01d0 * t
                         tst2 = norm + t
                         if (tst2 > tst1) cycle
                         exit
                       end do
                    end if
                    h(i,en) = -r / t
                  else
                    !
                    !  solve real equations
                    !
                     x = h(i,i+1)
                     y = h(i+1,i)
                     q = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i)
                     t = (x * s - zz * r) / q
                     h(i,en) = t
                     if (abs(x) > abs(zz)) then
                     h(i+1,en) = (-r - w * t) / x
                     else
                       h(i+1,en) = (-s - y * t) / zz
                     end if
                  end if
                  !
                  !  overflow control
                  !
                  t = abs(h(i,en))
                  if (t == 0.0_dp) cycle
                  tst1 = t
                  tst2 = tst1 + 1.0_dp/tst1
                  if (tst2 > tst1) cycle
                  do j = i, en
                     h(j,en) = h(j,en)/t
                  end do
               end do
               !
               !  end real vector
               !
             end if
          end do
          !
          !  end back substitution.
          !
          !  vectors of isolated roots
          !
          do i = 1, n
             if (i >= low .and. i <= igh) cycle

             do j = i, n
               z(i,j) = h(i,j)
             end do

          end do
          !
          !  multiply by transformation matrix to give vectors of original full matrix.
          !  for j=n step -1 until low do.
          !
          do jj = low, n
             j = n + low - jj
             m = min(j,igh)

             do i = low, igh
                zz = 0.0_dp
                do k = low, m
                  zz = zz + z(i,k) * h(k,j)
                end do
                z(i,j) = zz
            end do
          end do

    End Subroutine hqr2



    !  Subroutine htrib3 ( n, a, tau, m, zr, zi )
    !
    !! HTRIB3 determines eigenvectors by undoing the HTRID3 transformation.
    !
    !  Discussion:
    !
    !    HTRIB3 forms the eigenvectors of a complex hermitian
    !    matrix by back transforming those of the corresponding
    !    real symmetric tridiagonal matrix determined by HTRID3.
    !
    !  Modified:
    !
    !    05 February 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, is the order of the matrix.
    !
    !    Input, real(kind=dp) A(N,N), contains information about the unitary
    !    transformations used in the reduction by HTRID3.
    !
    !    Input, real(kind=dp) TAU(2,N), contains further information about the
    !    transformations.
    !
    !    Input, integer ::M, the number of eigenvectors to be back
    !    transformed.
    !
    !    Input/output, real(kind=dp) ZR(N,M), ZI(N,M).  On input, ZR contains
    !    the eigenvectors to be back transformed.  On output, ZR and ZI contain
    !    the real and imaginary parts of the transformed eigenvectors.
    !
    Subroutine htrib3 ( n, a, tau, m, zr, zi )
      integer,                      intent(in)     :: n
      real(kind=dp), dimension(n,n),intent(in)     :: a
      real(kind=dp), dimension(2,n),intent(in)     :: tau
      integer,                      intent(in)     :: m
      real(kind=dp), dimension(n,m),intent(in out) :: zr,zi

      integer :: i,j,k
      real(kind=dp) :: h,s,si

      if ( m == 0 ) return
    !
    !  Transform the eigenvectors of the real symmetric tridiagonal matrix
    !  to those of the hermitian tridiagonal matrix.
    !
      do k = 1, n
        do j = 1, m
          zi(k,j) = - zr(k,j) * tau(2,k)
          zr(k,j) =   zr(k,j) * tau(1,k)
        end do
      end do
    !
    !  Recover and apply the Householder matrices.
    !
      do i = 2, n

        h = a(i,i)

        if ( h /= 0.0_dp ) then

          do j = 1, m

            s = 0.0_dp
            si = 0.0_dp

            do k = 1, i - 1
              s =  s  + a(i,k) * zr(k,j) - a(k,i) * zi(k,j)
              si = si + a(i,k) * zi(k,j) + a(k,i) * zr(k,j)
            end do

            s = ( s / h ) / h
            si = ( si / h ) / h

            zr(1:i-1,j) = zr(1:i-1,j) - s * a(i,1:i-1) - si * a(1:i-1,i)
            zi(1:i-1,j) = zi(1:i-1,j) - si * a(i,1:i-1) + s * a(1:i-1,i)

          end do

        end if

      end do

    End Subroutine htrib3



    !    Subroutine htribk ( n, ar, ai, tau, m, zr, zi )
    !
    !! HTRIBK determines eigenvectors by undoing the HTRIDI transformation.
    !
    !  Discussion:
    !
    !    HTRIBK forms the eigenvectors of a complex hermitian matrix by back
    !    transforming those of the corresponding real symmetric tridiagonal
    !    matrix determined by HTRIDI.
    !
    !  Modified:
    !
    !    04 March 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input, real(kind=dp) AR(N,N), AI(N,N), contain information about
    !    the unitary transformations used in the reduction by HTRIDI in their
    !    full lower triangles, except for the diagonal of AR.
    !
    !    Input, real(kind=dp) TAU(2,N), contains further information about the
    !    transformations.
    !
    !    Input, integer ::M, the number of eigenvectors to be back
    !    transformed.
    !
    !    Input/output, real(kind=dp) ZR(N,M), ZI(N,M).  On input, ZR contains
    !    the eigenvectors to be back transformed.  On output, ZR and ZI contain
    !    the real and imaginary parts of the transformed eigenvectors.
    !
    Subroutine htribk ( n, ar, ai, tau, m, zr, zi )
      integer,                      intent(in)     :: n
      real(kind=dp), dimension(n,n),intent(in)     :: ar,ai
      real(kind=dp), dimension(2,n),intent(in)     :: tau
      integer,                      intent(in)     :: m
      real(kind=dp), dimension(n,m),intent(in out) :: zr,zi

      integer :: i,j,k,l
      real(kind=dp) :: h,s,si

      if ( m == 0 ) return
    !
    !  Transform the eigenvectors of the real symmetric tridiagonal matrix to
    !  those of the hermitian tridiagonal matrix.
    !
      do k = 1, n
        do j = 1, m
          zi(k,j) = - zr(k,j) * tau(2,k)
          zr(k,j) =   zr(k,j) * tau(1,k)
        end do
      end do
    !
    !  Recover and apply the Householder matrices.
    !
      do i = 2, n

        l = i - 1
        h = ai(i,i)

        if ( h /= 0.0_dp ) then

          do j = 1, m

            s = 0.0_dp
            si = 0.0_dp
            do k = 1, l
              s =  s  + ar(i,k) * zr(k,j) - ai(i,k) * zi(k,j)
              si = si + ar(i,k) * zi(k,j) + ai(i,k) * zr(k,j)
            end do

            s = ( s / h ) / h
            si = ( si / h ) / h

            zr(1:l,j) = zr(1:l,j) - s  * ar(i,1:l) - si * ai(i,1:l)
            zi(1:l,j) = zi(1:l,j) - si * ar(i,1:l) + s  * ai(i,1:l)

          end do

        end if

      end do

    End Subroutine htribk



    !    Subroutine htrid3 ( n, a, d, e, e2, tau )
    !
    !! HTRID3 tridiagonalizes a complex hermitian packed matrix.
    !
    !  Discussion:
    !
    !    HTRID3 reduces a complex hermitian matrix, stored as
    !    a single square array, to a real symmetric tridiagonal matrix
    !    using unitary similarity transformations.
    !
    !  Modified:
    !
    !    05 February 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input/output, real(kind=dp) A(N,N).  On input, the lower triangle of
    !    the complex hermitian input matrix.  The real parts of the matrix elements
    !    are stored in the full lower triangle of A, and the imaginary parts are
    !    stored in the transposed positions of the strict upper triangle of A.  No
    !    storage is required for the zero imaginary parts of the diagonal elements.
    !    On output, A contains information about the unitary transformations
    !    used in the reduction.
    !
    !    Output, real(kind=dp) D(N), the diagonal elements of the
    !    tridiagonal matrix.
    !
    !    Output, real(kind=dp) E(N), the subdiagonal elements of the tridiagonal
    !    matrix in E(2:N).  E(1) is set to zero.
    !
    !    Output, real(kind=dp) E2(N), the squares of the corresponding elements
    !    of E.  E2 may coincide with E if the squares are not needed.
    !
    !    Output, real(kind=dp) TAU(2,N), contains further information about the
    !    transformations.
    !
    Subroutine htrid3 ( n, a, d, e, e2, tau )
      integer,                      intent(in)     :: n
      real(kind=dp), dimension(n,n),intent(in out) :: a
      real(kind=dp), dimension(n),  intent(out)    :: d,e,e2
      real(kind=dp), dimension(2,n),intent(out)    :: tau

      real(kind=dp) :: f,fi,g,gi,h,hh,scal,si
      integer :: i,j,k

      tau(1,n) = 1.0_dp
      tau(2,n) = 0.0_dp

      do i = n, 1, -1

        h = 0.0_dp
        scal = 0.0_dp

        if ( i < 2 ) then

          e(i) = 0.0_dp
          e2(i) = 0.0_dp

        else
    !
    !  scal row.
    !
          do k = 1, i - 1
            scal = scal + abs ( a(i,k) ) + abs ( a(k,i) )
          end do

          if ( scal == 0.0_dp ) then

            tau(1,i-1) = 1.0_dp
            tau(2,i-1) = 0.0_dp
            e(i) = 0.0_dp
            e2(i) = 0.0_dp

          else

            do k = 1, i - 1
              a(i,k) = a(i,k) / scal
              a(k,i) = a(k,i) / scal
              h = h + a(i,k) * a(i,k) + a(k,i) * a(k,i)
            end do

            e2(i) = scal * scal * h
            g = sqrt ( h )
            e(i) = scal * g
            f = pythag ( a(i,i-1), a(i-1,i) )
    !
    !  Form next diagonal element of matrix T.
    !
            if ( f /= 0.0_dp ) then

              tau(1,i-1) = ( a(i-1,i) * tau(2,i) - a(i,i-1) * tau(1,i) ) / f
              si =         ( a(i,i-1) * tau(2,i) + a(i-1,i) * tau(1,i) ) / f
              h = h + f * g
              g = 1.0_dp + g / f
              a(i,i-1) = g * a(i,i-1)
              a(i-1,i) = g * a(i-1,i)

              if ( i == 2 ) then
                a(i,1:i-1) = scal * a(i,1:i-1)
                a(1:i-1,i) = scal * a(1:i-1,i)
                tau(2,i-1) = - si
                d(i) = a(i,i)
                a(i,i) = scal * sqrt ( h )
                cycle
              end if

            else

              tau(1,i-1) = - tau(1,i)
              si = tau(2,i)
              a(i,i-1) = g

            end if

            f = 0.0_dp

            do j = 1, i - 1

              g = 0.0_dp
              gi = 0.0_dp
    !
    !  Form element of A*U.
    !
              do k = 1, j - 1
                g  = g  + a(j,k) * a(i,k) + a(k,j) * a(k,i)
                gi = gi - a(j,k) * a(k,i) + a(k,j) * a(i,k)
              end do

              g  = g  + a(j,j) * a(i,j)
              gi = gi - a(j,j) * a(j,i)

              do k = j + 1, i - 1
                g  = g  + a(k,j) * a(i,k) - a(j,k) * a(k,i)
                gi = gi - a(k,j) * a(k,i) - a(j,k) * a(i,k)
              end do
    !
    !  Form element of P.
    !
              e(j) = g / h
              tau(2,j) = gi / h
              f = f + e(j) * a(i,j) - tau(2,j) * a(j,i)

            end do

            hh = f / ( h + h )
    !
    !  Form reduced A.
    !
            do j = 1, i - 1

              f = a(i,j)
              g = e(j) - hh * f
              e(j) = g
              fi = - a(j,i)
              gi = tau(2,j) - hh * fi
              tau(2,j) = - gi
              a(j,j) = a(j,j) - 2.0_dp * ( f * g + fi * gi )

              do k = 1, j - 1
                a(j,k) = a(j,k) &
                  - f * e(k) - g * a(i,k) + fi * tau(2,k) + gi * a(k,i)
                a(k,j) = a(k,j) &
                  - f * tau(2,k) - g * a(k,i) - fi * e(k) - gi * a(i,k)
              end do

            end do

            a(i,1:i-1) = scal * a(i,1:i-1)
            a(1:i-1,i) = scal * a(1:i-1,i)
            tau(2,i-1) = - si

          end if

        end if

        d(i) = a(i,i)
        a(i,i) = scal * sqrt ( h )

      end do

    End Subroutine htrid3



    !    Subroutine htridi ( n, ar, ai, d, e, e2, tau )
    !
    !! HTRIDI tridiagonalizes a complex hermitian matrix.
    !
    !  Discussion:
    !
    !    HTRIDI reduces a complex hermitian matrix to a real symmetric
    !    tridiagonal matrix using unitary similarity transformations.
    !
    !  Modified:
    !
    !    04 March 2018
    !
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input/output, real(kind=dp) AR(N,N), AI(N,N).  On input, the real
    !    and imaginary parts, respectively, of the complex hermitian input matrix.
    !    Only the lower triangle of the matrix need be supplied.
    !    On output, information about the unitary transformations used in the
    !    reduction in their full lower triangles.  Their strict upper triangles
    !    and the diagonal of AR are unaltered.
    !
    !    Output, real(kind=dp) D(N), the diagonal elements of the
    !    tridiagonal matrix.
    !
    !    Output, real(kind=dp) E(N), the subdiagonal elements of the tridiagonal
    !    matrix in its last N-1 positions.  E(1) is set to zero.
    !
    !    Output, real(kind=dp) E2(N), the squares of the corresponding elements
    !    of E.  E2 may coincide with E if the squares are not needed.
    !
    !    Output, real(kind=dp) TAU(2,N), contains further information about the
    !    transformations.
    !
    Subroutine htridi ( n, ar, ai, d, e, e2, tau )
      integer,                      intent(in)     :: n
      real(kind=dp), dimension(n,n),intent(in out) :: ar,ai
      real(kind=dp), dimension(n),  intent(out)    :: d,e,e2
      real(kind=dp), dimension(2,n),intent(out)    :: tau

      real(kind=dp) :: f,fi,g,gi,h,hh,scal,si
      integer :: i,j,k

      tau(1,n) = 1.0_dp
      tau(2,n) = 0.0_dp

      do i = 1, n
        d(i) = ar(i,i)
      end do

      do i = n, 1, -1

        h = 0.0_dp
        scal = 0.0_dp

        if ( i == 1 ) then

          e(i) = 0.0_dp
          e2(i) = 0.0_dp
          hh = d(i)
          d(i) = ar(i,i)
          ar(i,i) = hh
          ai(i,i) = scal * sqrt ( h )
    !
    !  scal row.
    !
        else

          do k = 1, i - 1
            scal = scal + abs ( ar(i,k) ) + abs ( ai(i,k) )
          end do

          if ( scal == 0.0_dp ) then
            tau(1,i-1) = 1.0_dp
            tau(2,i-1) = 0.0_dp
            e(i) = 0.0_dp
            e2(i) = 0.0_dp
            hh = d(i)
            d(i) = ar(i,i)
            ar(i,i) = hh
            ai(i,i) = scal * sqrt ( h )
            cycle
          end if

          ar(i,1:i-1) = ar(i,1:i-1) / scal
          ai(i,1:i-1) = ai(i,1:i-1) / scal

          do k = 1, i - 1
            h = h + ar(i,k) * ar(i,k) + ai(i,k) * ai(i,k)
          end do

          e2(i) = scal * scal * h
          g = sqrt ( h )
          e(i) = scal * g
          f = pythag ( ar(i,i-1), ai(i,i-1) )
    !
    !  Form next diagonal element of matrix T.
    !
          if ( f /= 0.0_dp ) then

            tau(1,i-1) = ( ai(i,i-1) * tau(2,i) - ar(i,i-1) * tau(1,i) ) / f
            si =         ( ar(i,i-1) * tau(2,i) + ai(i,i-1) * tau(1,i) ) / f
            h = h + f * g
            g = 1.0_dp + g / f
            ar(i,i-1) = g * ar(i,i-1)
            ai(i,i-1) = g * ai(i,i-1)

            if ( i == 2 ) then
              ar(i,1:i-1) = scal * ar(i,1:i-1)
              ai(i,1:i-1) = scal * ai(i,1:i-1)
              tau(2,i-1) = - si
              hh = d(i)
              d(i) = ar(i,i)
              ar(i,i) = hh
              ai(i,i) = scal * sqrt ( h )
              cycle
            end if

          else

            tau(1,i-1) = - tau(1,i)
            si = tau(2,i)
            ar(i,i-1) = g

          end if

          f = 0.0_dp

          do j = 1, i - 1

            g = 0.0_dp
            gi = 0.0_dp
    !
    !  Form element of A*U.
    !
            do k = 1, j
              g  = g  + ar(j,k) * ar(i,k) + ai(j,k) * ai(i,k)
              gi = gi - ar(j,k) * ai(i,k) + ai(j,k) * ar(i,k)
            end do

            do k = j + 1, i - 1
              g  = g  + ar(k,j) * ar(i,k) - ai(k,j) * ai(i,k)
              gi = gi - ar(k,j) * ai(i,k) - ai(k,j) * ar(i,k)
            end do
    !
    !  Form element of P.
    !
            e(j) = g / h
            tau(2,j) = gi / h
            f = f + e(j) * ar(i,j) - tau(2,j) * ai(i,j)

          end do

          hh = f / ( h + h )
    !
    !  Form the reduced A.
    !
          do j = 1, i - 1

            f = ar(i,j)
            g = e(j) - hh * f
            e(j) = g
            fi = - ai(i,j)
            gi = tau(2,j) - hh * fi
            tau(2,j) = - gi

            do k = 1, j
              ar(j,k) = ar(j,k) - f * e(k) - g * ar(i,k) + fi * tau(2,k) &
                + gi * ai(i,k)
              ai(j,k) = ai(j,k) - f * tau(2,k) - g * ai(i,k) - fi * e(k) &
                - gi * ar(i,k)
            end do

          end do

          ar(i,1:i-1) = scal * ar(i,1:i-1)
          ai(i,1:i-1) = scal * ai(i,1:i-1)
          tau(2,i-1) = - si

          hh = d(i)
          d(i) = ar(i,i)
          ar(i,i) = hh
          ai(i,i) = scal * sqrt ( h )

        end if

      end do

    End Subroutine htridi



    !    Subroutine i4vec_print ( n, a, title )
    !
    !! I4VEC_PRINT prints an I4VEC.
    !
    !  Discussion:
    !
    !    An I4VEC is a vector of I4's.
    !
    !  Modified:
    !
    !    02 May 2010
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Arguments:
    !
    !    Input, integer ::N, the number of components of the vector.
    !
    !    Input, integer ::A(N), the vector to be printed.
    !
    !    Input, character ( len = * ) TITLE, a title.
    !
    Subroutine i4vec_print ( n, a, title )
      integer,               intent(in) :: n
      integer, dimension(n), intent(in) :: a
      character ( len = * ), intent(in) :: title

      integer ::i

      if ( 0 < len_trim ( title ) ) then
        write(*, '(a)' ) ' '
        write(*, '(a)' ) trim ( title )
      end if

      write(*, '(a)' ) ' '
      do i = 1, n
        write (*, '(2x,i8,a,2x,i12)' ) i, ':', a(i)
      end do

    End Subroutine i4vec_print



    !    Subroutine Imtql1 ( N, D, E, Ierr )
    !
    !! IMTQL1 computes all eigenvalues of a symmetric tridiagonal matrix.
    !
    !  Discussion:
    !
    !    This routine finds the eigenvalues of a symmetric
    !    tridiagonal matrix by the implicit QL method.
    !
    !
    !  Modified:
    !
    !    08 February 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input/output, real(kind=dp) D(N).  On input, the diagonal elements of
    !    the matrix.  On output, the eigenvalues in ascending order.  If an error
    !    exit is made, the eigenvalues are correct and ordered for indices
    !    1,2,...IERR-1, but may not be the smallest eigenvalues.
    !
    !    Input/output, real(kind=dp) E(N).  On input, the subdiagonal elements
    !    of the matrix in its last N-1 positions.  E(1) is arbitrary.  On output,
    !    E has been overwritten.
    !
    !    Output, integer ::IERR, error flag.
    !    0, normal return,
    !    J, if the J-th eigenvalue has not been determined after 30 iterations.
    !
    Subroutine Imtql1 ( n, d, e, Ierr )
      integer,                     intent(in)     :: n
      real(kind=dp), dimension(n), intent(in out) :: d
      real(kind=dp), dimension(n), intent(in out) :: e
      integer,                     intent(out)    :: ierr

      real(kind=dp) :: b,c,f,g,p,r,s,tst1,tst2
      integer :: i,its,l,m
      logical :: skip

      ierr = 0

      if ( n == 1 ) return

      do i = 2, n
        e(i-1) = e(i)
      end do
      e(n) = 0.0_dp

      do l = 1, n

        its = 0
    !
    !  Look for a small sub-diagonal element.
    !
        do

          m = l

          do m = l, n - 1

            tst1 = abs ( d(m) ) + abs ( d(m+1) )
            tst2 = tst1 + abs ( e(m) )

            if ( tst2 == tst1 ) then
              exit
            end if

          end do
    !
    !  Order the eigenvalues.
    !
          p = d(l)

          if ( m == l ) then

            do i = l, 1, -1

              if ( i == 1 ) then
                d(i) = p
                exit
              end if

              if ( d(i-1) <= p ) then
                d(i) = p
                exit
              end if

             d(i) = d(i-1)

            end do

            exit

          else

            if ( 30 <= its ) then
              ierr = l
              return
            end if

            its = its + 1
    !
    !  Form shift.
    !
            g = ( d(l+1) - p ) / ( 2.0_dp * e(l) )
            r = pythag ( g, 1.0_dp )
            g = d(m) - p + e(l) / ( g + sign ( r, g ) )
            s = 1.0_dp
            c = 1.0_dp
            p = 0.0_dp

            skip = .false.

            do i = m - 1, l, -1

              f = s * e(i)
              b = c * e(i)
              r = pythag ( f, g )
              e(i+1) = r
    !
    !  Recover from underflow.
    !
              if ( r == 0.0_dp ) then
                d(i+1) = d(i+1) - p
                e(m) = 0.0_dp
                skip = .true.
                exit
              end if

              s = f / r
              c = g / r
              g = d(i+1) - p
              r = ( d(i) - g ) * s + 2.0_dp * c * b
              p = s * r
              d(i+1) = g + p
              g = c * r - b

            end do

            if ( .not. skip ) then
              d(l) = d(l) - p
              e(l) = g
              e(m) = 0.0_dp
            end if

          end if

        end do

      end do

    End Subroutine Imtql1



    !    Subroutine imtql2 ( n, d, e, z, ierr )
    !
    !! IMTQL2 computes all eigenvalues/vectors of a symmetric tridiagonal matrix.
    !
    !  Discussion:
    !
    !    IMTQL2 finds the eigenvalues and eigenvectors of a symmetric tridiagonal
    !    matrix by the implicit QL method.
    !
    !    The eigenvectors of a full symmetric matrix can also be found if TRED2
    !    has been used to reduce this full matrix to tridiagonal form.
    !
    !  Modified:
    !
    !    04 February 2018
    !
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input/output, real(kind=dp) D(N).  On input, the diagonal elements of
    !    the input matrix.  On output, the eigenvalues in ascending order.  If an
    !    error exit is made, the eigenvalues are correct but
    !    unordered for indices 1,2,...,IERR-1.
    !
    !    Input/output, real(kind=dp) E(N).  On input, the subdiagonal elements
    !    of the input matrix in E(2:N).  E(1) is arbitrary.  On output, E is
    !    overwritten.
    !
    !    Input/output, real(kind=dp) Z(N,N).  On input, the transformation
    !    matrix produced in the reduction by TRED2, if performed.  If the
    !    eigenvectors of the tridiagonal matrix are desired, Z must contain the
    !    identity matrix.  On output, Z contains orthonormal eigenvectors of the
    !    symmetric tridiagonal (or full) matrix.  If an error exit is made, Z
    !    contains the eigenvectors associated with the stored eigenvalues.
    !
    !    Output, integer ::IERR, error flag.
    !    0, for normal return,
    !    J, if the J-th eigenvalue has not been determined after 30 iterations.
    !
    Subroutine imtql2 ( n, d, e, z, ierr )
      integer,                     intent(in)     :: n
      real(kind=dp), dimension(n), intent(in out) :: d
      real(kind=dp), dimension(n), intent(in out) :: e
      real(kind=dp),dimension(n,n),intent(in out) :: z
      integer,                     intent(out)    :: ierr

      real(kind=dp) :: b,c,f,g,p,r,s,tst1,tst2
      integer       :: i,its,j,k,l,m
      real(kind=dp), dimension(n) :: t

      ierr = 0

      if ( n == 1 )  return

      do i = 2, n
        e(i-1) = e(i)
      end do
      e(n) = 0.0_dp

      do l = 1, n

        its = 0
    !
    !  Look for a small sub-diagonal element.
    !
        do

          m = l

          do m = l, n - 1

            tst1 = abs ( d(m) ) + abs ( d(m+1) )
            tst2 = tst1 + abs ( e(m) )

            if ( tst2 == tst1 ) then
              exit
            end if

          end do

          p = d(l)

          if ( m == l ) then
            exit
          end if

          if ( 30 <= its ) then
            ierr = l
            return
          end if

          its = its + 1
    !
    !  Form shift.
    !
          g = ( d(l+1) - p ) / ( 2.0_dp * e(l) )
          r = pythag ( g, 1.0_dp )
          g = d(m) - p + e(l) / ( g + sign ( r, g ) )
          s = 1.0_dp
          c = 1.0_dp
          p = 0.0_dp

          do i = m - 1, l, -1

            f = s * e(i)
            b = c * e(i)
            r = pythag ( f, g )
            e(i+1) = r
    !
    !  Recover from underflow.
    !
            if ( r == 0.0_dp ) then
              d(i+1) = d(i+1) - p
              e(m) = 0.0_dp
              cycle
            end if

            s = f / r
            c = g / r
            g = d(i+1) - p
            r = ( d(i) - g ) * s + 2.0_dp * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
    !
    !  Form vector.
    !
            do k = 1, n
              f = z(k,i+1)
              z(k,i+1) = s * z(k,i) + c * f
              z(k,i) = c * z(k,i) - s * f
            end do

          end do

          d(l) = d(l) - p
          e(l) = g
          e(m) = 0.0_dp

        end do

      end do
    !
    !  Order eigenvalues and eigenvectors.
    !
      do i = 1, n - 1

        k = i
        p = d(i)

        do j = i + 1, n
          if ( d(j) < p ) then
            k = j
            p = d(j)
          end if
        end do

        if ( k /= i ) then

          d(k) = d(i)
          d(i) = p

          t(1:n)   = z(1:n,i)
          z(1:n,i) = z(1:n,k)
          z(1:n,k) = t(1:n)

        end if

      end do

    End Subroutine imtql2



    !    Subroutine imtqlv ( n, d, e, e2, w, ind, ierr )
    !
    !! IMTQLV computes all eigenvalues of a real symmetric tridiagonal matrix.
    !
    !  Discussion:
    !
    !    IMTQLV finds the eigenvalues of a symmetric tridiagonal matrix by
    !    the implicit QL method and associates with them their corresponding
    !    submatrix indices.
    !
    !  Modified:
    !
    !    03 February 2018
    !
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input, real(kind=dp) D(N), the diagonal elements of the input matrix.
    !
    !    Input, real(kind=dp) E(N), the subdiagonal elements of the input matrix
    !    in E(2:N).  E(1) is arbitrary.
    !
    !    Input/output, real(kind=dp) E2(N).  On input, the squares of the
    !    corresponding elements of E.  E2(1) is arbitrary.  On output, elements of
    !    E2 corresponding to elements of E regarded as negligible have been
    !    replaced by zero, causing the matrix to split into a direct sum of
    !    submatrices.  E2(1) is also set to zero.
    !
    !    Output, real(kind=dp) W(N), the eigenvalues in ascending order.  If an
    !    error exit is made, the eigenvalues are correct and ordered for
    !    indices 1,2,...IERR-1, but may not be the smallest eigenvalues.
    !
    !    Output, integer ::IND(N), the submatrix indices associated with
    !    the corresponding eigenvalues in W: 1 for eigenvalues belonging to the
    !    first submatrix from the top, 2 for those belonging to the second
    !    submatrix, and so on.
    !
    !    Output, integer ::IERR, error flag.
    !    0, for normal return,
    !    J, if the J-th eigenvalue has not been determined after 30 iterations.
    !
    Subroutine imtqlv ( n, d, e, e2, w, ind, ierr )
      integer,                     intent(in)     :: n
      real(kind=dp), dimension(n), intent(in)     :: d
      real(kind=dp), dimension(n), intent(in)     :: e
      real(kind=dp), dimension(n), intent(in out) :: e2
      real(kind=dp), dimension(n), intent(out)    :: w
      integer,       dimension(n), intent(out)    :: ind
      integer,                     intent(out)    :: ierr

      real(kind=dp) :: b,c,f,g,p,r,s,tst1,tst2
      integer       :: i,its,k,l,m, tag
      real(kind=dp), dimension(n) :: rv1
      logical :: skip

      ierr = 0
      k = 0
      tag = 0
      w(1:n) = d(1:n)
      e2(1) = 0.0_dp
      rv1(1:n-1) = e(2:n)
      rv1(n) = 0.0_dp

      do l = 1, n

        its = 0
    !
    !  Look for a small sub-diagonal element.
    !
        do

          do m = l, n

            if ( m == n ) then
              exit
            end if

            tst1 = abs ( w(m) ) + abs ( w(m+1) )
            tst2 = tst1 + abs ( rv1(m) )

            if ( tst2 == tst1 ) then
              exit
            end if
    !
    !  Guard against underflowed element of E2.
    !
            if ( e2(m+1) == 0.0_dp ) then
              k = m
              tag = tag + 1
              exit
            end if

          end do

          if ( k < m ) then
            if ( m /= n ) then
              e2(m+1) = 0.0_dp
            end if
            k = m
            tag = tag + 1
          end if

          p = w(l)

          if ( m == l ) then

            do i = l, 1, -1
              if ( i == 1 ) then
                w(i) = p
                ind(i) = tag
              else if ( w(i-1) <= p ) then
                w(i) = p
                ind(i) = tag
                exit
              else
                w(i) = w(i-1)
                ind(i) = ind(i-1)
              end if
            end do

            exit

          else

            if ( 30 <= its ) then
              ierr = l
              return
            end if

            its = its + 1
    !
    !  Form shift.
    !
            g = ( w(l+1) - p ) / ( 2.0_dp * rv1(l) )
            r = pythag ( g, 1.0_dp )
            g = w(m) - p + rv1(l) / ( g + sign ( r, g ) )
            s = 1.0_dp
            c = 1.0_dp
            p = 0.0_dp

            skip = .false.

            do i = m - 1, l, -1
              f = s * rv1(i)
              b = c * rv1(i)
              r = pythag ( f, g )
              rv1(i+1) = r

              if ( r == 0.0_dp ) then
                w(i+1) = w(i+1) - p
                rv1(m) = 0.0_dp
                skip = .true.
                exit
              end if

              s = f / r
              c = g / r
              g = w(i+1) - p
              r = ( w(i) - g ) * s + 2.0_dp * c * b
              p = s * r
              w(i+1) = g + p
              g = c * r - b
            end do

            if ( .not. skip ) then
              w(l) = w(l) - p
              rv1(l) = g
              rv1(m) = 0.0_dp
            end if

          end if

        end do

      end do

    End Subroutine imtqlv



    !    Subroutine invit ( n, a, wr, wi, select, mm, m, z, ierr )
    !
    !! INVIT computes eigenvectors of a real upper Hessenberg matrix.
    !
    !  Discussion:
    !
    !    INVIT finds those eigenvectors of a real upper Hessenberg
    !    matrix corresponding to specified eigenvalues, using inverse iteration.
    !
    !  Modified:
    !
    !    02 March 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input, real(kind=dp) A(N,N), the Hessenberg matrix.
    !
    !    Input/output, real(kind=dp) WR(N), WI(N).  On input, the real and
    !    imaginary parts, respectively, of the eigenvalues of the matrix.  The
    !    eigenvalues must be stored in a manner identical to that of HQR,
    !    which recognizes possible splitting of the matrix.  On output,
    !    WR may have been altered since close eigenvalues are perturbed
    !    slightly in searching for independent eigenvectors.
    !
    !    Input/output, logical SELECT(N).  On input, specifies the eigenvectors
    !    to be found.  The eigenvector corresponding to the J-th eigenvalue is
    !    specified by setting SELECT(J) to TRUE.  On output, SELECT may have been
    !    altered.  If the elements corresponding to a pair of conjugate complex
    !    eigenvalues were each initially set to TRUE, the program resets the
    !    second of the two elements to FALSE.
    !
    !    Input, integer ::MM, an upper bound for the number of columns
    !    required to store the eigenvectors to be found.  Note that two columns are
    !    required to store the eigenvector corresponding to a complex eigenvalue.
    !
    !    Output, integer ::M, the number of columns actually used to store
    !    the eigenvectors.
    !
    !    Output, real(kind=dp) Z(N,MM), the real and imaginary parts of the
    !    eigenvectors.  If the next selected eigenvalue is real, the next column
    !    of Z contains its eigenvector.  If the eigenvalue is complex, the next
    !    two columns of Z contain the real and imaginary parts of its eigenvector.
    !    The eigenvectors are normalized so that the component of largest
    !    magnitude is 1.  Any vector which fails the acceptance test is set to zero.
    !
    !    Output, integer ::IERR, error flag.
    !    0, for normal return,
    !    -(2*N+1), if more than MM columns of Z are necessary to store the
    !      eigenvectors corresponding to the specified eigenvalues.
    !    -K, if the iteration corresponding to the K-th value fails,
    !    -(N+K), if both error situations occur.
    !
    Subroutine invit ( n, a, wr, wi, select, mm, m, z, ierr )
      integer,                       intent(in)     :: n
      real(kind=dp), dimension(n,n), intent(in)     :: a
      real(kind=dp), dimension(n),   intent(in out) :: wr,wi
      logical,       dimension(n),   intent(in out) :: select
      integer,                       intent(in)     :: mm
      integer,                       intent(out)    :: m
      real(kind=dp), dimension(n,mm),intent(out)    :: z
      integer,                       intent(out)    :: ierr

      real(kind=dp) :: eps3,growto,ilambd,norm,normv,rlambd,ukroot,w,x,y,t

      integer :: i,ip,its,j,k,l,mp,n1,ns,s,uk
      logical :: repeat
      real(kind=dp),dimension(n,n) :: rm1
      real(kind=dp),dimension(n)   :: rv1,rv2

      ierr = 0
      uk = 0
      s = 1
      !
      !  The value of IP is:
      !
      !   0, real eigenvalue;
      !   1, first of conjugate complex pair;
      !  -1, second of conjugate complex pair.
      !
      ip = 0
      n1 = n - 1

      do_k: do k = 1, n

        if ( wi(k) /= 0.0_dp .and. 0 <= ip ) then
          ip = 1
          if ( select(k) .and. select(k+1) ) then
            select(k+1) = .false.
          end if
        end if

        if ( .not. select(k) ) then
          if ( ip == -1 ) then
            ip = 0
          end if
          if ( ip == 1 ) then
            ip = -1
          end if
          cycle do_k
        end if

        if ( wi(k) /= 0.0_dp ) then
          s = s + 1
        end if

        if ( mm < s ) then

          if ( ierr /= 0 ) then
            ierr = ierr - n
          end if
          if ( ierr == 0 ) then
            ierr = - ( 2 * n + 1 )
          end if
          m = s - 1 - abs ( ip )
          return

        end if

        if ( uk < k ) then
          !
          !  Check for possible splitting.
          !
          do uk = k, n

            if ( uk == n ) then
              exit
            end if
            if ( a(uk+1,uk) == 0.0_dp ) then
              exit
            end if
          end do
          !
          !  Compute infinity norm of leading UK by UK (Hessenberg) matrix.
          !
          norm = 0.0_dp
          mp = 1

          do i = 1, uk

            x = sum ( abs ( a(i,mp:uk) ) )
            norm = max ( norm, x )
            mp = i

          end do
          !
          !  EPS3 replaces zero pivot in decomposition and close roots are modified
          !  by EPS3.
          !
          if ( norm == 0.0_dp ) then
            norm = 1.0_dp
          end if

          eps3 = abs ( norm ) * epsilon ( eps3 )
          !
          !  GROWTO is the criterion for the growth.
          !
          ukroot = uk
          ukroot = sqrt ( ukroot )
          growto = 0.1_dp / ukroot

        end if

        rlambd = wr(k)
        ilambd = wi(k)
        !
        !  Perturb eigenvalue if it is close to any previous eigenvalue.
        !
        if ( 1 < k ) then

          repeat = .true.

          do while ( repeat )

            repeat = .false.

            do i = k - 1, 1, -1
              if ( select(i) .and. &
                   abs ( wr(i) - rlambd ) < eps3 .and. &
                   abs ( wi(i) - ilambd ) < eps3 ) then
               rlambd = rlambd + eps3
               repeat = .true.
               exit
              end if
            end do

          end do

          wr(k) = rlambd
          wr(k+ip) = rlambd

        end if
        !
        !  Form upper Hessenberg A - rlambd*I (transposed) and initial real vector.
        !
        mp = 1

        do i = 1, uk
          rm1(mp:uk,i) = a(i,mp:uk)
          rm1(i,i) = rm1(i,i) - rlambd
          mp = i
          rv1(i) = eps3
        end do

        its = 0
        !
        !  Real eigenvalue.
        !
        !  Triangular decomposition with interchanges, replacing zero pivots by eps3.
        !
        if ( ilambd == 0.0_dp ) then

          do i = 2, uk

            mp = i - 1

            if ( abs ( rm1(mp,mp) ) < abs ( rm1(mp,i) ) ) then

              do j = mp, uk
                t         = rm1(j,i)
                rm1(j,i)  = rm1(j,mp)
                rm1(j,mp) = t
              end do

            end if

            if ( rm1(mp,mp) == 0.0_dp ) then
              rm1(mp,mp) = eps3
            end if

            x = rm1(mp,i) / rm1(mp,mp)

            if ( x /= 0.0_dp ) then
              rm1(i:uk,i) = rm1(i:uk,i) - x * rm1(i:uk,mp)
            end if

          end do

          if ( rm1(uk,uk) == 0.0_dp ) then
            rm1(uk,uk) = eps3
          end if
          !
          !  Back substitution for real vector.
          !
          do i = uk, 1, -1
            y = rv1(i)
            do j = i + 1, uk
              y = y - rm1(j,i) * rv1(j)
            end do
            rv1(i) = y / rm1(i,i)
          end do

          go to 740

        end if
        !
        !  Complex eigenvalue.
        !
        !  Triangular decomposition with interchanges,
        !  replacing zero pivots by EPS3.
        !  Store imaginary parts in upper triangle starting at (1,3)
        !
        ns = n - s
        z(1,s-1) = - ilambd
        z(1,s) = 0.0_dp

        if ( n /= 2 ) then
          rm1(1,3) = - ilambd
          z(1,s-1) = 0.0_dp
          rm1(1,4:n) = 0.0_dp
        end if

        do i = 2, uk

          mp = i - 1
          w = rm1(mp,i)

          if ( i < n ) then
            t = rm1(mp,i+1)
          else if ( i == n ) then
            t = z(mp,s-1)
          end if

          x = rm1(mp,mp) * rm1(mp,mp) + t * t

          if ( x < w * w ) then

            x = rm1(mp,mp) / w
            y = t / w
            rm1(mp,mp) = w

            if ( i < n ) then
              rm1(mp,i+1) = 0.0_dp
            else if ( i == n ) then
              z(mp,s-1) = 0.0_dp
            end if

            do j = i, uk

              w = rm1(j,i)
              rm1(j,i) = rm1(j,mp) - x * w
              rm1(j,mp) = w

              if ( n1 <= j ) then
                l = j - ns
                z(i,l) = z(mp,l) - y * w
                z(mp,l) = 0.0_dp
              else
                rm1(i,j+2) = rm1(mp,j+2) - y * w
                rm1(mp,j+2) = 0.0_dp
              end if

            end do

            rm1(i,i) = rm1(i,i) - y * ilambd

            if ( n1 <= i ) then
              l = i - ns
              z(mp,l) = -ilambd
              z(i,l) = z(i,l) + x * ilambd
            else
              rm1(mp,i+2) = -ilambd
              rm1(i,i+2) = rm1(i,i+2) + x * ilambd
            end if

          else

            if ( x == 0.0_dp ) then
              rm1(mp,mp) = eps3
              if ( i < n ) then
                rm1(mp,i+1) = 0.0_dp
              else if ( i == n ) then
                z(mp,s-1) = 0.0_dp
              end if
              t = 0.0_dp
              x = eps3 ** 2
            end if

            w = w / x
            x = rm1(mp,mp) * w
            y = - t * w

            do j = i, uk

              if ( n1 <= j ) then
                l = j - ns
                t = z(mp,l)
                z(i,l) = - x * t - y * rm1(j,mp)
              else
                t = rm1(mp,j+2)
                rm1(i,j+2) = - x * t - y * rm1(j,mp)
              end if

              rm1(j,i) = rm1(j,i) - x * rm1(j,mp) + y * t

            end do

            if ( n1 <= i ) then
              l = i - ns
              z(i,l) = z(i,l) - ilambd
            else
              rm1(i,i+2) = rm1(i,i+2) - ilambd
            end if

          end if

        end do

        if ( n1 <= uk ) then
          l = uk - ns
          t = z(uk,l)
        else
          t = rm1(uk,uk+2)
        end if

        if ( rm1(uk,uk) == 0.0_dp .and. t == 0.0_dp ) then
          rm1(uk,uk) = eps3
        end if
        !
        !  Back substitution for complex vector.
        !
    660 continue

           do i = uk, 1, -1

             x = rv1(i)
             y = 0.0_dp

             do j = i + 1, uk

               if ( n1 <= j ) then
                 t = z(i,j-ns)
               else
                 t = rm1(i,j+2)
               end if

               x = x - rm1(j,i) * rv1(j) + t * rv2(j)
               y = y - rm1(j,i) * rv2(j) - t * rv1(j)

             end do

             if ( n1 <= i ) then
               t = z(i,i-ns)
             else
               t = rm1(i,i+2)
             end if

             call cdiv ( x, y, rm1(i,i), t, rv1(i), rv2(i) )

           end do
           !
           !  Acceptance test for real or complex eigenvector and normalization.
           !
       740 continue

           its = its + 1
           norm = 0.0_dp
           normv = 0.0_dp

           do i = 1, uk

             if ( ilambd == 0.0_dp ) then
               x = abs ( rv1(i) )
             else
               x = pythag ( rv1(i), rv2(i) )
             end if

             if ( normv < x )  then
               normv = x
               j = i
             end if

             norm = norm + x

           end do
           !
           !  Choose a new starting vector.
           !
           if ( norm < growto ) then

             if ( uk <= its ) then
               j = 1
               ierr = - k
               do i = j, n
                 z(i,s) = 0.0_dp
                 if ( ilambd /= 0.0_dp ) then
                   z(i,s-1) = 0.0_dp
                 end if
               end do

               s = s + 1

               if ( ip == -1 ) then
                 ip = 0
               end if

               if ( ip == 1 ) then
                 ip = -1
               end if

               cycle do_k

             end if

             x = ukroot
             y = eps3 / ( x + 1.0_dp )

             rv1(1) = eps3
             rv1(2:uk) = y

             j = uk - its + 1
             rv1(j) = rv1(j) - eps3 * x

             if ( ilambd == 0.0_dp ) then
               do i = uk, 1, -1
                 y = rv1(i)
                 do j = i + 1, uk
                   y = y - rm1(j,i) * rv1(j)
                 end do
                 rv1(i) = y / rm1(i,i)
               end do
               go to 740
             end if

             go to 660

           end if
        !
        !  Accept vector.
        !
        x = rv1(j)
        if ( ilambd == 0.0_dp ) then
          x = 1.0_dp / x
        else
          y = rv2(j)
        end if

        do i = 1, uk
          if ( ilambd == 0.0_dp ) then
            z(i,s) = rv1(i) * x
          else
            call cdiv ( rv1(i), rv2(i), x, y, z(i,s-1), z(i,s) )
          end if
        end do

        j = uk + 1
        do i = j, n
          z(i,s) = 0.0_dp
          if ( ilambd /= 0.0_dp ) then
            z(i,s-1) = 0.0_dp
          end if
        end do

        s = s + 1

        if ( ip == -1 ) then
          ip = 0
        end if

        if ( ip == 1 ) then
          ip = -1
        end if

      end do do_k !k

      m = s - 1 - abs ( ip )

    End Subroutine invit



    !    Subroutine minfit ( nm, m, n, a, w, ip, b, ierr )
    !
    !! MINFIT: least squares problem for a real overdetermined linear system.
    !
    !  Discussion:
    !
    !    MINFIT is part of an algorithm for solving general linear
    !    systems of the form A*X=B.
    !
    !    It determines the singular value decomposition
    !      A = U * S * V'
    !    of a real M by N rectangular matrix, forming U' * B
    !    rather than U.  Householder bidiagonalization and a variant of the
    !    QR algorithm are used.
    !
    !  Modified:
    !
    !    01 February 2018
    !
    !  Arguments:
    !
    !    Input, integer ::NM, the leading dimension of the
    !    two-dimensional arrays.  NM must be at least as large as the maximum
    !    of M and N.
    !
    !    Input, integer ::M, the number of rows of A and B.
    !
    !    Input, integer ::N, the number of columns of A, and the order
    !    of V.
    !
    !    Input/output, real(kind=dp) A(NM,N). On input, the rectangular
    !    coefficient matrix.  On output, A has been overwritten by the orthogonal
    !    matrix V of the decomposition in its first N rows and columns.  If an
    !    error exit is made, the columns of V corresponding to indices of correct
    !    singular values should be correct.
    !
    !    Output, real(kind=dp) W(N), the singular values of A.  These are the
    !    diagonal elements of S.  They are unordered.  If an error exit is made, the
    !    singular values should be correct for indices IERR+1, IERR+2,...,N.
    !
    !    Input, integer ::IP, is the number of columns of B.  IP can
    !    be zero.
    !
    !    Input/output, real(kind=dp) B(NM,IP).  On input, the constant column
    !    matrix.  On output, B has been overwritten by U'*B.  If an error exit is
    !    made, the rows of U'*B corresponding to indices of correct singular values
    !    should be correct.
    !
    !    Output, integer ::IERR, error flag.
    !    0, for normal return,
    !    K, if the K-th singular value has not been determined after 30 iterations.
    !
    Subroutine minfit ( nm, m, n, a, w, ip, b, ierr )
      integer,                        intent(in)     :: nm,m,n
      real(kind=dp), dimension(nm,n), intent(in out) :: a
      real(kind=dp), dimension(n),    intent(out)    :: w
      integer,                        intent(in)     :: ip
      real(kind=dp), dimension(nm,ip),intent(in out) :: b
      integer,                        intent(out)    :: ierr

      real(kind=dp) :: c,f,g,h,x,y,z,s,tst1,tst2,scal
      integer       :: i,its,j,k,l,m1
      logical       :: skip
      real(kind=dp), dimension(n) :: rv1

      ierr = 0
    !
    !  Householder reduction to bidiagonal form.
    !
      g = 0.0_dp
      scal = 0.0_dp
      x = 0.0_dp

      do i = 1, n

        l = i + 1
        rv1(i) = scal * g
        g = 0.0_dp
        s = 0.0_dp
        scal = 0.0_dp

        if ( i <= m ) then

          scal = sum ( abs ( a(i:m,i) ) )

          if ( scal /= 0.0_dp ) then

            a(i:m,i) = a(i:m,i) / scal

            s = s + sum ( a(i:m,i) ** 2 )

            f = a(i,i)
            g = - sign ( sqrt ( s ), f )
            h = f * g - s
            a(i,i) = f - g

            do j = i + 1, n
              s = dot_product ( a(i:m,i), a(i:m,j) )
              f = s / h
              a(i:m,j) = a(i:m,j) + f * a(i:m,i)
            end do

            do j = 1, ip
              s = dot_product ( a(i:m,i), b(i:m,j) )
              b(i:m,j) = b(i:m,j) + s * a(i:m,i) / h
            end do

            a(i:m,i) = scal * a(i:m,i)

          end if

        end if

        w(i) = scal * g
        g = 0.0_dp
        s = 0.0_dp
        scal = 0.0_dp

        if ( i <= m .and. i /= n ) then

          scal = scal + sum ( abs ( a(i,i+1:n) ) )

          if ( scal /= 0.0_dp ) then

            a(i,i+1:n) = a(i,i+1:n) / scal

            s = s + sum ( a(i,i+1:n)**2 )

            f = a(i,i+1)
            g = - sign ( sqrt ( s ), f )
            h = f * g - s
            a(i,i+1) = f - g
            rv1(i+1:n) = a(i,i+1:n) / h

            do j = i + 1, m
              s = dot_product ( a(j,i+1:n), a(i,i+1:n) )
              a(j,i+1:n) = a(j,i+1:n) + s * rv1(i+1:n)
            end do

            a(i,i+1:n) = scal * a(i,i+1:n)

          end if

        end if

        x = max ( x, abs ( w(i) ) + abs ( rv1(i) ) )

      end do
    !
    !  Accumulation of right-hand transformations.
    !
      do i = n, 1, -1

        if ( i < n ) then

          if ( g /= 0.0_dp ) then

            a(i+1:n,i) = ( a(i,i+1:n) / a(i,i+1) ) / g

            do j = i + 1, n
              s = dot_product ( a(i,i+1:n), a(i+1:n,j) )
              a(i+1:n,j) = a(i+1:n,j) + s * a(i+1:n,i)
            end do

          end if

          a(i,i+1:n) = 0.0_dp
          a(i+1:n,i) = 0.0_dp

        end if

        a(i,i) = 1.0_dp
        g = rv1(i)

      end do

      if ( m < n .and. ip /= 0 ) then
        m1 = m + 1
        b(m+1:n,1:ip) = 0.0_dp
      end if
    !
    !  Diagonalization of the bidiagonal form.
    !
      tst1 = x

      do k = n, 1, -1

        its = 0

        do
    !
    !  Test for splitting.
    !
          skip = .false.

          do l = k, 1, -1

            tst2 = tst1 + abs ( rv1(l) )

            if ( tst2 == tst1 ) then
              skip = .true.
              exit
            end if

            tst2 = tst1 + abs ( w(l-1) )

            if ( tst2 == tst1 ) then
              exit
            end if

          end do
    !
    !  Cancellation of RV1(l) if l greater than 1.
    !
          if ( .not. skip ) then

            c = 0.0_dp
            s = 1.0_dp

            do i = l, k

              f = s * rv1(i)
              rv1(i) = c * rv1(i)
              tst2 = tst1 + abs ( f )

              if ( tst2 == tst1 ) then
                exit
              end if

              g = w(i)
              h = pythag ( f, g )
              w(i) = h
              c = g / h
              s = - f / h

              do j = 1, ip
                y = b(l-1,j)
                z = b(i,j)
                b(l-1,j) =  y * c + z * s
                b(i,j) = - y * s + z * c
              end do

            end do

          end if
    !
    !  Test for convergence.
    !
          z = w(k)

          if ( l == k ) then

            if ( z < 0.0_dp ) then
              w(k) = - z
              a(1:n,k) = - a(1:n,k)
            end if

            return

          end if
    !
    !  Shift from bottom 2 by 2 minor.
    !
          if ( 30 <= its ) then
            ierr = k
            return
          end if

          its = its + 1
          x = w(l)
          y = w(k-1)
          g = rv1(k-1)
          h = rv1(k)
          f = 0.5_dp * ( ( ( g + z ) / h ) * ( ( g - z ) / y ) + y / h - h / y )
          g = pythag ( f, 1.0_dp )
          f = x - ( z / x ) * z + ( h / x ) * ( y / ( f + sign ( g, f ) ) - h )
    !
    !  Next QR transformation.
    !
          c = 1.0_dp
          s = 1.0_dp

          do i = l + 1, k

            g = rv1(i)
            y = w(i)
            h = s * g
            g = c * g
            z = pythag ( f, h )
            rv1(i-1) = z
            c = f / z
            s = h / z
            f =   x * c + g * s
            g = - x * s + g * c
            h = y * s
            y = y * c

            do j = 1, n
              x = a(j,i-1)
              z = a(j,i)
              a(j,i-1) =  x * c + z * s
              a(j,i) =  - x * s + z * c
            end do

            z = pythag ( f, h )
            w(i-1) = z

            if ( z /= 0.0_dp ) then
              c = f / z
              s = h / z
            end if

            f =   c * g + s * y
            x = - s * g + c * y

            do j = 1, ip
              y = b(i-1,j)
              z = b(i,j)
              b(i-1,j) =  y * c + z * s
              b(i,j) =  - y * s + z * c
            end do

          end do

          rv1(l) = 0.0_dp
          rv1(k) = f
          w(k) = x

        end do

      end do

    End Subroutine minfit



    !    Subroutine ortbak ( n, low, igh, a, ort, m, z )
    !
    !! ORTBAK determines eigenvectors by undoing the ORTHES transformation.
    !
    !  Discussion:
    !
    !    ORTBAK forms the eigenvectors of a real general matrix by back
    !    transforming those of the corresponding upper Hessenberg matrix
    !    determined by ORTHES.
    !
    !
    !  Modified:
    !
    !    04 March 2018
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    !    Klema, Moler.
    !    FORTRAN90 version by John Burkardt
    !
    !  Reference:
    !
    !    James Wilkinson, Christian Reinsch,
    !    Handbook for Automatic Computation,
    !    Volume II, Linear Algebra, Part 2,
    !    Springer, 1971,
    !    ISBN: 0387054146,
    !    LC: QA251.W67.
    !
    !    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    !    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    !    Matrix Eigensystem Routines, EISPACK Guide,
    !    Lecture Notes in Computer Science, Volume 6,
    !    Springer Verlag, 1976,
    !    ISBN13: 978-3540075462,
    !    LC: QA193.M37.
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input, integer ::LOW, IGH, are determined by the balancing
    !    routine BALANC.  If BALANC has not been used, set LOW = 1 and IGH equal
    !    to the order of the matrix.
    !
    !    Input, real(kind=dp) A(N,IGH), contains information about the
    !    orthogonal transformations used in the reduction by ORTHES in its strict
    !    lower triangle.
    !
    !    Input/output, real(kind=dp) ORT(IGH), contains further information
    !    about the transformations used in the reduction by ORTHES.  On output, ORT
    !    has been altered.
    !
    !    Input, integer ::M, the number of columns of Z to be back
    !    transformed.
    !
    !    Input/output, real(kind=dp) Z(N,M).  On input, the real and imaginary
    !    parts of the eigenvectors to be back transformed in the first M columns.
    !    On output, the real and imaginary parts of the transformed eigenvectors.
    !
    Subroutine ortbak ( n, low, igh, a, ort, m, z )
      integer,                        intent(in)     :: n, low, igh,m
      real(kind=dp), dimension(n,igh),intent(in out) :: a
      real(kind=dp), dimension(igh),  intent(in out) :: ort
      real(kind=dp), dimension(n,m),  intent(in out) :: z


      real(kind=dp) :: g
      integer :: i,j,mp

      if ( m == 0 ) return

      do mp = igh - 1, low + 1, -1

        if ( a(mp,mp-1) /= 0.0_dp ) then

          ort(mp+1:igh) = a(mp+1:igh,mp-1)

          do j = 1, m

            g = dot_product ( ort(mp:igh), z(mp:igh,j) )

            g = ( g / ort(mp) ) / a(mp,mp-1)

            do i = mp, igh
              z(i,j) = z(i,j) + g * ort(i)
            end do

          end do

        end if

      end do

    End Subroutine ortbak



    !    Subroutine orthes ( n, low, igh, a, ort )
    !
    !! ORTHES transforms a real general matrix to upper Hessenberg form.
    !
    !  Discussion:
    !
    !    ORTHES is given a real general matrix, and reduces a submatrix
    !    situated in rows and columns LOW through IGH to upper Hessenberg form by
    !    orthogonal similarity transformations.
    !
    !  Modified:
    !
    !    04 March 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input, integer ::LOW, IGH, are determined by the balancing
    !    routine BALANC.  If BALANC has not been used, set LOW = 1 and IGH = N.
    !
    !    Input/output, real(kind=dp) A(N,N).  On input, the matrix.  On output,
    !    the Hessenberg matrix.  Information about the orthogonal transformations
    !    used in the reduction is stored in the remaining triangle under the
    !    Hessenberg matrix.
    !
    !    Output, real(kind=dp) ORT(IGH), contains further information about the
    !    transformations.
    !
    Subroutine orthes ( n, low, igh, a, ort )
      integer,                        intent(in)     :: n, low, igh
      real(kind=dp), dimension(n,n),  intent(in out) :: a
      real(kind=dp), dimension(igh),  intent(out)    :: ort

      real(kind=dp) :: f,g,h,scal
      integer :: i,j,m

      do m = low + 1, igh - 1

        h = 0.0_dp
        ort(m) = 0.0_dp
        scal = 0.0_dp
    !
    !  scal the column.
    !
        do i = m, igh
          scal = scal + abs ( a(i,m-1) )
        end do

        if ( scal /= 0.0_dp ) then

          do i = igh, m, -1
            ort(i) = a(i,m-1) / scal
            h = h + ort(i) * ort(i)
          end do

          g = - sign ( sqrt ( h ), ort(m) )
          h = h - ort(m) * g
          ort(m) = ort(m) - g
    !
    !  Form (I-(U*Ut)/h) * A.
    !
          do j = m, n

            f = 0.0_dp
            do i = igh, m, -1
              f = f + ort(i) * a(i,j)
            end do
            f = f / h

            do i = m, igh
              a(i,j) = a(i,j) - f * ort(i)
            end do

          end do
    !
    !  Form (I-(u*ut)/h) * A * (I-(u*ut)/h).
    !
          do i = 1, igh

            f = 0.0_dp
            do j = igh, m, -1
              f = f + ort(j) * a(i,j)
            end do

            a(i,m:igh) = a(i,m:igh) - f * ort(m:igh) / h

          end do

          ort(m) = scal * ort(m)
          a(m,m-1) = scal * g

        end if

      end do

    End Subroutine orthes


    !!----  Subroutine Ortran(n,low,igh,a,ort,z)
    !!----   integer,                          intent(in)     :: n,low,igh
    !!----   real(kind = dp),dimension(n,igh), intent(in)     :: a
    !!----   real(kind = dp),dimension(igh),   intent(in out) :: ort
    !!----   real(kind = dp),dimension(n,n),   intent(out)    :: z
    !!----
    !!----  Subroutine adapted from EisPack, originally called RG_ORT
    !!----  ORTRAN accumulates similarity transformations generated by ORTHES.
    !!----
    !!----  Discussion:
    !!----
    !!----    ORTRAN accumulates the orthogonal similarity transformations used in
    !!----    the reduction of a real general matrix to upper Hessenberg form by ORTHES.
    !!----
    !!----  Modified:
    !!----
    !!----    01 February 2018, 29 November 2019 (JRC)
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input, integer ::LOW, IGH, are determined by the balancing
    !!----    routine BALANC.  If BALANC has not been used, set LOW = 1, IGH = N.
    !!----
    !!----    Input, real(kind=dp) A(N,IGH), contains information about the
    !!----    orthogonal transformations used in the reduction by ORTHES in its strict
    !!----    lower triangle.
    !!----
    !!----    Input/output, real(kind=dp) ORT(IGH), contains further information
    !!----    about the transformations used in the reduction by ORTHES.  On output, ORT
    !!----    has been further altered.
    !!----
    !!----    Output, real(kind=dp) Z(N,N), contains the transformation matrix
    !!----    produced in the reduction by ORTHES.
    !!----
    Subroutine Ortran(n,low,igh,a,ort,z)
      integer,                          intent(in)     :: n,low,igh
      real(kind = dp),dimension(n,igh), intent(in)     :: a
      real(kind = dp),dimension(igh),   intent(in out) :: ort
      real(kind = dp),dimension(n,n),   intent(out)    :: z

      real(kind = dp):: g
      integer :: j,mp
      !
      !  Initialize Z to the identity matrix.
      !
      call r8mat_identity ( n, z )
      if ( igh - low < 2 ) then
        return
      end if
      do mp = igh - 1, low + 1, -1
        if ( a(mp,mp-1) /= 0.0_dp ) then
          ort(mp+1:igh) = a(mp+1:igh,mp-1)
          do j = mp, igh
            g = dot_product ( ort(mp:igh), z(mp:igh,j) )
            g = ( g / ort(mp) ) / a(mp,mp-1)
            z(mp:igh,j) = z(mp:igh,j) + g * ort(mp:igh)
          end do
        end if
      end do

    End Subroutine Ortran



    !    function pythag ( a, b ) result(pythag_val)
    !
    !! PYTHAG computes SQRT ( A * A + B * B ) carefully.
    !
    !  Discussion:
    !
    !    The formula
    !
    !      PYTHAG = sqrt ( A * A + B * B )
    !
    !    is reasonably accurate, but can fail if, for example, A^2 is larger
    !    than the machine overflow.  The formula can lose most of its accuracy
    !    if the sum of the squares is very large or very small.
    !
    !  Modified:
    !
    !    18 October 2009
    !
    !  Arguments:
    !
    !    Input, real(kind=dp) A, B, the two legs of a right triangle.
    !
    !    Output, real(kind=dp) PYTHAG, the length of the hypotenuse.
    !
    function pythag ( a, b ) result(pythag_val)
      real(kind=dp), intent(in) :: a,b
      real(kind=dp) :: pythag_val

      real(kind=dp) :: p,r,s,t,u

      p = max ( abs ( a ), abs ( b ) )

      if ( p /= 0.0_dp ) then

        r = ( min ( abs ( a ), abs ( b ) ) / p )**2

        do

          t = 4.0_dp + r

          if ( t == 4.0_dp ) then
            exit
          end if

          s = r / t
          u = 1.0_dp + 2.0_dp * s
          p = u * p
          r = ( s / u )**2 * r

        end do

      end if

      pythag_val = p

    End Function pythag



    !    Subroutine qzhes ( n, a, b, matz, z )
    !
    !! QZHES carries out transformations for a generalized eigenvalue problem.
    !
    !  Discussion:
    !
    !    QZHES is the first step of the QZ algorithm
    !    for solving generalized matrix eigenvalue problems.
    !
    !    QZHES accepts a pair of real general matrices and
    !    reduces one of them to upper Hessenberg form and the other
    !    to upper triangular form using orthogonal transformations.
    !    it is usually followed by QZIT, QZVAL and, possibly, QZVEC.
    !
    !  Modified:
    !
    !    10 February 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrices.
    !
    !    Input/output, real(kind=dp) A(N,N).  On input, the first real general
    !    matrix.  On output, A has been reduced to upper Hessenberg form.  The
    !    elements below the first subdiagonal have been set to zero.
    !
    !    Input/output, real(kind=dp) B(N,N).  On input, a real general matrix.
    !    On output, B has been reduced to upper triangular form.  The elements
    !    below the main diagonal have been set to zero.
    !
    !    Input, logical MATZ, should be true if the right hand transformations
    !    are to be accumulated for later use in computing eigenvectors.
    !
    !    Output, real(kind=dp) Z(N,N), contains the product of the right hand
    !    transformations if MATZ is true.
    !
    Subroutine qzhes ( n, a, b, matz, z )
      integer,                       intent(in)     :: n
      real(kind=dp), dimension(n,n), intent(in out) :: a,b
      logical,                       intent(in)     :: matz
      real(kind=dp), dimension(n,n), intent(out)    :: z

      integer       :: i,j,k,l
      real(kind=dp) :: r,rho,s,t,u1,u2,v1,v2
    !
    !  Set Z to the identity matrix.
    !
      if ( matz ) then
        call r8mat_identity ( n, z )
      end if
    !
    !  Reduce B to upper triangular form.
    !
      if ( n <= 1 ) return

      do l = 1, n - 1

        s = sum ( abs ( b(l+1:n,l) ) )

        if ( s /= 0.0_dp ) then

          s = s + abs ( b(l,l) )
          b(l:n,l) = b(l:n,l) / s

          r = sqrt ( sum ( b(l:n,l)**2 ) )
          r = sign ( r, b(l,l) )
          b(l,l) = b(l,l) + r
          rho = r * b(l,l)

          do j = l + 1, n

            t = dot_product ( b(l:n,l), b(l:n,j) )

            b(l:n,j) = b(l:n,j) - t * b(l:n,l) / rho

          end do

          do j = 1, n

            t = dot_product ( b(l:n,l), a(l:n,j) )

            a(l:n,j) = a(l:n,j) - t * b(l:n,l) / rho

          end do

          b(l,l) = - s * r
          b(l+1:n,l) = 0.0_dp

        end if

      end do
    !
    !  Reduce A to upper Hessenberg form, while keeping B triangular.
    !
      do k = 1, n - 2

        do l = n - 1, k + 1, -1
    !
    !  Zero A(l+1,k).
    !
          s = abs ( a(l,k) ) + abs ( a(l+1,k) )

          if ( s /= 0.0_dp ) then

            u1 = a(l,k) / s
            u2 = a(l+1,k) / s
            r = sign ( sqrt ( u1 * u1 + u2 * u2 ), u1 )
            v1 = - ( u1 + r ) / r
            v2 = - u2 / r
            u2 = v2 / v1

            do j = k, n
              t = a(l,j) + u2 * a(l+1,j)
              a(l,j) = a(l,j) + t * v1
              a(l+1,j) = a(l+1,j) + t * v2
            end do

            a(l+1,k) = 0.0_dp

            do j = l, n
              t = b(l,j) + u2 * b(l+1,j)
              b(l,j) = b(l,j) + t * v1
              b(l+1,j) = b(l+1,j) + t * v2
            end do
    !
    !  Zero B(l+1,l).
    !
            s = abs ( b(l+1,l+1) ) + abs ( b(l+1,l) )

            if ( s /= 0.0 ) then

              u1 = b(l+1,l+1) / s
              u2 = b(l+1,l) / s
              r = sign ( sqrt ( u1 * u1 + u2 * u2 ), u1 )
              v1 =  - ( u1 + r ) / r
              v2 = - u2 / r
              u2 = v2 / v1

              do i = 1, l + 1
                t = b(i,l+1) + u2 * b(i,l)
                b(i,l+1) = b(i,l+1) + t * v1
                b(i,l) = b(i,l) + t * v2
              end do

              b(l+1,l) = 0.0_dp

              do i = 1, n
                t = a(i,l+1) + u2 * a(i,l)
                a(i,l+1) = a(i,l+1) + t * v1
                a(i,l) = a(i,l) + t * v2
              end do

              if ( matz ) then

                do i = 1, n
                  t = z(i,l+1) + u2 * z(i,l)
                  z(i,l+1) = z(i,l+1) + t * v1
                  z(i,l) = z(i,l) + t * v2
                end do

              end if

            end if

          end if

        end do

      end do

    End Subroutine qzhes



    !    Subroutine qzit ( n, a, b, eps1, matz, z, ierr )
    !
    !! QZIT carries out iterations to solve a generalized eigenvalue problem.
    !
    !  Discussion:
    !
    !    QZIT is the second step of the QZ algorithm for generalized matrix
    !    eigenvalue problems.
    !
    !    QZIT accepts a pair of real matrices, one of them
    !    in upper Hessenberg form and the other in upper triangular form.
    !    It reduces the Hessenberg matrix to quasi-triangular form using
    !    orthogonal transformations while maintaining the triangular form
    !    of the other matrix.  It is usually preceded by QZHES and
    !    followed by QZVAL and, possibly, QZVEC.
    !
    !
    !  Modified:
    !
    !    27 February 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrices.
    !
    !    Input/output, real(kind=dp) A(N,N).  On input, a real upper Hessenberg
    !    matrix.  On output, A has been reduced to quasi-triangular form.  The
    !    elements below the first subdiagonal are still zero and no two consecutive
    !    subdiagonal elements are nonzero.
    !
    !    Input/output, real(kind=dp) B(N,N).  On input, a real upper triangular
    !    matrix.  On output, B is still in upper triangular form, although its
    !    elements have been altered.  The location B(N,1) is used to store EPS1
    !    times the norm of B for later use by QZVAL and QZVEC.
    !
    !    Input, real(kind=dp) EPS1, a tolerance used to determine negligible
    !    elements.  EPS1 = 0.0 (or negative) may be input, in which case an element
    !    will be neglected only if it is less than roundoff error times the
    !    norm of its matrix.  If the input EPS1 is positive, then an element
    !    will be considered negligible if it is less than EPS1 times the norm
    !    of its matrix.  A positive value of EPS1 may result in faster execution,
    !    but less accurate results.
    !
    !    Input, logical MATZ, should be true if the right hand transformations
    !    are to be accumulated for later use in computing eigenvectors.
    !
    !    Input/output, real(kind=dp) Z(N,N).  If MATZ is false, Z is not
    !    referenced.  Otherwise, on input, the transformation matrix produced in the
    !    reduction by QZHES, if performed, or else the identity matrix.  On output,
    !    Z contains the product of the right hand transformations for both steps.
    !
    !    Output, integer ::IERR, error flag.
    !    0, for normal return,
    !    J, if the limit of 30*N iterations is exhausted while the J-th
    !      eigenvalue is being sought.
    !
    Subroutine qzit ( n, a, b, eps1, matz, z, ierr )
      integer,                       intent(in)     :: n
      real(kind=dp), dimension(n,n), intent(in out) :: a,b
      real(kind=dp),                 intent(in)     :: eps1
      logical,                       intent(in)     :: matz
      real(kind=dp), dimension(n,n), intent(in out) :: z
      integer,                       intent(out)    :: ierr


      real(kind=dp) :: a1,a11,a12,a2,a21,a22,a3,a33,a34,a43,a44,ani,     &
                       anorm ,b11,b12,b22,b33,b34,b44,bni,bnorm,ep,epsa, &
                       epsb,r,s,sh,t,u2,u3,v1,v2,v3,u1
      integer :: i,ish,itn,its,en,enm2,j,k,k1,k2,km1,l,l1,ld,ll,lm1,lor1,na,enorn

      logical :: notlas,skip



      ierr = 0
    !
    !  Compute EPSA and EPSB.
    !
      anorm = 0.0_dp
      bnorm = 0.0_dp

      do i = 1, n

        if ( i == 1 ) then
          ani = 0.0_dp
        else
          ani = abs ( a(i,i-1) )
        end if

        bni = 0.0_dp

        do j = i, n
          ani = ani + abs ( a(i,j) )
          bni = bni + abs ( b(i,j) )
        end do

        anorm = max ( anorm, ani )
        bnorm = max ( bnorm, bni )

      end do

      if ( anorm == 0.0_dp ) then
        anorm = 1.0_dp
      end if

      if ( bnorm == 0.0_dp ) then
        bnorm = 1.0_dp
      end if

      ep = eps1

      if ( ep <= 0.0_dp ) then
        ep = epsilon ( ep )
      end if

      epsa = ep * anorm
      epsb = ep * bnorm
    !
    !  Reduce A to quasi-triangular form, while keeping B triangular.
    !
      lor1 = 1
      enorn = n
      en = n
      itn = 30 * n
    !
    !  Begin QZ step.
    !
      if ( en <= 2 ) then
        if ( 1 < n ) then
          b(n,1) = epsb
        end if
        return
      end if

      if ( .not. matz ) then
        enorn = en
      end if

      its = 0
      na = en - 1
      enm2 = na - 1

      do

        ish = 2
    !
    !  Check for convergence or reducibility.
    !
        skip = .false.

        do l = en, 1, -1
          lm1 = l - 1
          if ( l == 1 ) then
            ld = l
            skip = .true.
            exit
          end if
          if ( abs ( a(l,l-1) ) <= epsa ) then
            exit
          end if
        end do

    90  continue

        if ( .not. skip ) then

          a(l,l-1) = 0.0_dp
    !
    !  1-by-1 or 2-by-2 block isolated.
    !
          if ( na <= l ) then
            en = l - 1
            if ( en <= 2 ) then
              if ( 1 < n ) then
                b(n,1) = epsb
              end if
              return
            end if

            if ( .not. matz ) then
              enorn = en
            end if

            its = 0
            na = en - 1
            enm2 = na - 1
            cycle
          end if

          ld = l

        else

          skip = .false.

        end if
    !
    !  Check for small top of B.
    !
    100 continue

        l1 = l + 1
        b11 = b(l,l)

        if ( abs ( b11 ) <= epsb ) then

          b(l,l) = 0.0_dp
          s = abs ( a(l,l) ) + abs ( a(l+1,l) )
          u1 = a(l,l) / s
          u2 = a(l+1,l) / s
          r = sign ( sqrt ( u1 * u1 + u2 * u2 ), u1 )
          v1 = - ( u1 + r ) / r
          v2 = - u2 / r
          u2 = v2 / v1

          do j = l, enorn
            t = a(l,j) + u2 * a(l+1,j)
            a(l,j) = a(l,j) + t * v1
            a(l+1,j) = a(l+1,j) + t * v2
            t = b(l,j) + u2 * b(l+1,j)
            b(l,j) = b(l,j) + t * v1
            b(l+1,j) = b(l+1,j) + t * v2
          end do

          if ( 1 < l ) then
            a(l,lm1) = - a(l,lm1)
          end if
          lm1 = l
          l = l1
          go to 90

        end if

        a11 = a(l,l) / b11
        a21 = a(l1,l) / b11

        if ( ish == 1 ) then

          a1 = a11 - sh
          a2 = a21
          if ( l /= ld ) then
            a(l,lm1) = - a(l,lm1)
          end if
    !
    !  Set error: not all eigenvalues have converged after 30*N iterations.
    !  Save EPSB for use by QZVAL and QZVEC.
    !
        else if ( itn == 0 ) then

          ierr = en
          if ( 1 < n ) then
            b(n,1) = epsb
          end if
          return

        else if ( its == 10 ) then

          a1 = 0.0_dp
          a2 = 1.0_dp
          a3 = 1.1605_dp

        else
    !
    !  Determine type of shift.
    !
          b22 = b(l1,l1)
          if ( abs ( b22 ) < epsb ) then
            b22 = epsb
          end if

          b33 = b(na,na)
          if ( abs ( b33 ) < epsb ) then
            b33 = epsb
          end if

          b44 = b(en,en)
          if ( abs ( b44 ) < epsb ) then
            b44 = epsb
          end if

          a33 = a(na,na) / b33
          a34 = a(na,en) / b44
          a43 = a(en,na) / b33
          a44 = a(en,en) / b44
          b34 = b(na,en) / b44
          t = 0.5_dp * ( a43 * b34 - a33 - a44 )
          r = t * t + a34 * a43 - a33 * a44
    !
    !  Determine double shift zeroth column of A.
    !
          if ( r < 0.0_dp ) then

            a12 = a(l,l1) / b22
            a22 = a(l1,l1) / b22
            b12 = b(l,l1) / b22
            a1 = ( ( a33 - a11 ) * ( a44 - a11 ) - a34 * a43 + a43 * b34 * a11 ) &
              / a21 + a12 - a11 * b12
            a2 = ( a22 - a11 ) - a21 * b12 - ( a33 - a11 ) &
              - ( a44 - a11 ) + a43 * b34
            a3 = a(l1+1,l1) / b22

          else
    !
    !  Determine single shift zeroth column of A.
    !
            ish = 1
            r = sqrt ( r )
            sh = - t + r
            s = - t - r
            if ( abs ( s - a44 ) < abs ( sh - a44 ) ) then
              sh = s
            end if
    !
    !  Look for two consecutive small sub-diagonal elements of A.
    !
            do l = enm2, ld, -1
              if ( l == ld ) then
                exit
              end if
              lm1 = l - 1
              l1 = l + 1
              t = a(l,l)

              if ( epsb < abs ( b(l,l) ) ) then
                t = t - sh * b(l,l)
              end if

              if ( abs ( a(l,l-1) ) <= abs ( t / a(l+1,l) ) * epsa ) then
                go to 100
              end if

            end do

            a1 = a11 - sh
            a2 = a21

            if ( l /= ld ) then
              a(l,lm1) = - a(l,lm1)
            end if

          end if

        end if

        its = its + 1
        itn = itn - 1
        if ( .not. matz ) then
          lor1 = ld
        end if
    !
    !  Main loop.
    !
        do k = l, na

          notlas = ( k /= na .and. ish == 2 )
          k1 = k + 1
          k2 = k + 2
          km1 = max ( k - 1, l )
          ll = min ( en, k1 + ish )

          if ( .not. notlas ) then
    !
    !  Zero A(k+1,k-1).
    !
            if ( k /= l ) then
              a1 = a(k,km1)
              a2 = a(k1,km1)
            end if

            s = abs ( a1 ) + abs ( a2 )

            if ( s == 0.0_dp ) then
              exit
            end if

            u1 = a1 / s
            u2 = a2 / s
            r = sign ( sqrt ( u1 * u1 + u2 * u2 ), u1 )
            v1 = - ( u1 + r ) / r
            v2 = - u2 / r
            u2 = v2 / v1

            do j = km1, enorn
              t = a(k,j) + u2 * a(k1,j)
              a(k,j) = a(k,j) + t * v1
              a(k1,j) = a(k1,j) + t * v2
              t = b(k,j) + u2 * b(k1,j)
              b(k,j) = b(k,j) + t * v1
              b(k1,j) = b(k1,j) + t * v2
            end do

            if ( k /= l ) then
              a(k1,km1) = 0.0_dp
            end if
    !
    !  Zero A(k+1,k-1) and A(k+2,k-1).
    !
          else

            if ( k /= l ) then
              a1 = a(k,km1)
              a2 = a(k1,km1)
              a3 = a(k2,km1)
            end if

            s = abs ( a1 ) + abs ( a2 ) + abs ( a3 )

            if ( s == 0.0_dp ) then
              cycle
            end if

            u1 = a1 / s
            u2 = a2 / s
            u3 = a3 / s
            r = sign ( sqrt ( u1 * u1 + u2 * u2 + u3 * u3 ), u1 )
            v1 = - ( u1 + r ) / r
            v2 = - u2 / r
            v3 = - u3 / r
            u2 = v2 / v1
            u3 = v3 / v1

            do j = km1, enorn
              t = a(k,j) + u2 * a(k1,j) + u3 * a(k2,j)
              a(k,j) = a(k,j) + t * v1
              a(k1,j) = a(k1,j) + t * v2
              a(k2,j) = a(k2,j) + t * v3
              t = b(k,j) + u2 * b(k1,j) + u3 * b(k2,j)
              b(k,j) = b(k,j) + t * v1
              b(k1,j) = b(k1,j) + t * v2
              b(k2,j) = b(k2,j) + t * v3
            end do

            if ( k /= l ) then
              a(k1,km1) = 0.0_dp
              a(k2,km1) = 0.0_dp
            end if
    !
    !  Zero B(k+2,k+1) and B(k+2,k).
    !
            s = abs ( b(k2,k2) ) + abs ( b(k2,k1) ) + abs ( b(k2,k) )

            if ( s /= 0.0_dp ) then

              u1 = b(k2,k2) / s
              u2 = b(k2,k1) / s
              u3 = b(k2,k) / s
              r = sign ( sqrt ( u1 * u1 + u2 * u2 + u3 * u3 ), u1 )
              v1 = - ( u1 + r ) / r
              v2 = - u2 / r
              v3 = - u3 / r
              u2 = v2 / v1
              u3 = v3 / v1

              do i = lor1, ll
                t = a(i,k2) + u2 * a(i,k1) + u3 * a(i,k)
                a(i,k2) = a(i,k2) + t * v1
                a(i,k1) = a(i,k1) + t * v2
                a(i,k) = a(i,k) + t * v3
                t = b(i,k2) + u2 * b(i,k1) + u3 * b(i,k)
                b(i,k2) = b(i,k2) + t * v1
                b(i,k1) = b(i,k1) + t * v2
                b(i,k) = b(i,k) + t * v3
              end do

              b(k2,k) = 0.0_dp
              b(k2,k1) = 0.0_dp

              if ( matz ) then

                do i = 1, n
                  t = z(i,k2) + u2 * z(i,k1) + u3 * z(i,k)
                  z(i,k2) = z(i,k2) + t * v1
                  z(i,k1) = z(i,k1) + t * v2
                  z(i,k) = z(i,k) + t * v3
                end do

              end if

            end if

          end if
    !
    !  Zero B(k+1,k).
    !
          s = abs ( b(k1,k1) ) + abs ( b(k1,k) )

          if ( s /= 0.0_dp ) then

            u1 = b(k1,k1) / s
            u2 = b(k1,k) / s
            r = sign ( sqrt ( u1 * u1 + u2 * u2 ), u1 )
            v1 = - ( u1 + r ) / r
            v2 = - u2 / r
            u2 = v2 / v1

            do i = lor1, ll
              t = a(i,k1) + u2 * a(i,k)
              a(i,k1) = a(i,k1) + t * v1
              a(i,k) = a(i,k) + t * v2
              t = b(i,k1) + u2 * b(i,k)
              b(i,k1) = b(i,k1) + t * v1
              b(i,k) = b(i,k) + t * v2
            end do

            b(k1,k) = 0.0_dp

            if ( matz ) then

              do i = 1, n
                t = z(i,k1) + u2 * z(i,k)
                z(i,k1) = z(i,k1) + t * v1
                z(i,k) = z(i,k) + t * v2
              end do

            end if

          end if

        end do

      end do

    End Subroutine qzit



    !    Subroutine qzval ( n, a, b, alfr, alfi, beta, matz, z )
    !
    !! QZVAL computes eigenvalues for a generalized eigenvalue problem.
    !
    !  Discussion:
    !
    !    QZVAL is the third step of the QZ algorithm for solving generalized
    !    matrix eigenvalue problems.
    !
    !    QZVAL accepts a pair of real matrices, one of them in quasi-triangular
    !    form and the other in upper triangular form.  It reduces the
    !    quasi-triangular matrix further, so that any remaining 2-by-2 blocks
    !    correspond to pairs of complex eigenvalues, and returns quantities whose
    !    ratios give the generalized eigenvalues.  It is usually preceded by QZHES
    !    and QZIT and may be followed by QZVEC.
    !
    !  Modified:
    !
    !    26 February 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrices.
    !
    !    Input/output, real(kind=dp) A(N,N).  On input, a real upper
    !    quasi-triangular matrix.  On output, A has been reduced further to a
    !    quasi-triangular matrix in which all nonzero subdiagonal elements
    !    correspond to pairs of complex eigenvalues.
    !
    !    Input/output, real(kind=dp) B(N,N).  On input, a real upper triangular
    !    matrix.  In addition, location B(n,1) contains the tolerance quantity EPSB
    !    computed and saved in QZIT.  On output, B is still in upper triangular
    !    form, although its elements have been altered.  B(N,1) is unaltered.
    !
    !    Output, real(kind=dp) ALFR(N), ALFI(N), the real and imaginary parts of
    !    the diagonal elements of the triangular matrix that would be obtained
    !    if A were reduced completely to triangular form by unitary
    !    transformations.  Non-zero values of ALFI occur in pairs, the first
    !    member positive and the second negative.
    !
    !    Output, real(kind=dp) BETA(N), the diagonal elements of the
    !    corresponding B, normalized to be real and non-negative.  The generalized
    !    eigenvalues are then the ratios (ALFR + I * ALFI) / BETA.
    !
    !    Input, logical MATZ, should be true if the right hand transformations
    !    are to be accumulated for later use in computing eigenvectors, and
    !    false otherwise.
    !
    !    Input/output, real(kind=dp) Z(N,N), is only used if MATZ is true.
    !    On input, the transformation matrix produced in the reductions by QZHES
    !    and QZIT, if performed, or else the identity matrix.  On output,
    !    the product of the right hand transformations for all three steps.
    !
    Subroutine qzval ( n, a, b, alfr, alfi, beta, matz, z )
      integer,                       intent(in)     :: n
      real(kind=dp), dimension(n,n), intent(in out) :: a,b
      real(kind=dp), dimension(n),   intent(out)    :: alfr, alfi, beta
      logical,                       intent(in)     :: matz
      real(kind=dp), dimension(n,n), intent(in out) :: z


      real(kind=dp) :: a1,a11,a11i,a11r,a12,a12i,a12r,a1i,a2,a21,a22,a22i,   &
                       a22r,a2i,an,b11,b12,b22,bn,c,cq,cz,d,di,dr,e,ei,epsb, &
                       r,s,sqi,sqr,ssi,ssr,szi,szr,t,ti,tr,u1,u2,v1,v2
      integer :: en,i,isw,j,na

      epsb = b(n,1)
      isw = 1
    !
    !  Find eigenvalues of quasi-triangular matrices.
    !
      do en = n, 1, -1

        na = en - 1

        if ( isw == 2 ) then
          isw = 3 - isw
          cycle
        end if
    !
    !  1-by-1 block, one real root.
    !
        if ( en == 1 .or. a(en,na) == 0.0_dp ) then
          alfr(en) = a(en,en)
          if ( b(en,en) < 0.0_dp ) then
            alfr(en) = - alfr(en)
          end if
          beta(en) = abs ( b(en,en) )
          alfi(en) = 0.0_dp
          cycle
        end if
    !
    !  2-by-2 block.
    !
        if ( abs ( b(na,na) ) <= epsb ) then

          a1 = a(na,na)
          a2 = a(en,na)

        else

          if ( abs ( b(en,en) ) <= epsb ) then

            a1 = a(en,en)
            a2 = a(en,na)
            bn = 0.0_dp

          else

            an = abs ( a(na,na) ) + abs ( a(na,en) ) + abs ( a(en,na) ) &
              + abs ( a(en,en) )
            bn = abs ( b(na,na) ) + abs ( b(na,en) ) + abs ( b(en,en) )
            a11 = a(na,na) / an
            a12 = a(na,en) / an
            a21 = a(en,na) / an
            a22 = a(en,en) / an
            b11 = b(na,na) / bn
            b12 = b(na,en) / bn
            b22 = b(en,en) / bn
            e = a11 / b11
            ei = a22 / b22
            s = a21 / ( b11 * b22 )
            t = ( a22 - e * b22 ) / b22

            if ( abs ( ei ) < abs ( e ) ) then
              e = ei
              t = ( a11 - e * b11 ) / b11
            end if

            c = 0.5_dp * ( t - s * b12 )
            d = c * c + s * ( a12 - e * b12 )
    !
    !  Two complex roots.
    !
            if ( d < 0.0_dp ) then

              e = e + c
              ei = sqrt ( - d )
              a11r = a11 - e * b11
              a11i = ei * b11
              a12r = a12 - e * b12
              a12i = ei * b12
              a22r = a22 - e * b22
              a22i = ei * b22

              if ( abs ( a21 ) + abs ( a22r ) + abs ( a22i ) <= &
                   abs ( a11r ) + abs ( a11i ) + abs ( a12r ) + abs ( a12i ) ) then
                a1 = a12r
                a1i = a12i
                a2 = - a11r
                a2i = - a11i
              else
                a1 = a22r
                a1i = a22i
                a2 = - a21
                a2i = 0.0_dp
              end if
    !
    !  Choose complex Z.
    !
              cz = sqrt ( a1 * a1 + a1i * a1i )

              if ( cz /= 0.0_dp ) then
                szr = ( a1 * a2 + a1i * a2i ) / cz
                szi = ( a1 * a2i - a1i * a2 ) / cz
                r = sqrt ( cz * cz + szr * szr + szi * szi )
                cz = cz / r
                szr = szr / r
                szi = szi / r
              else
                szr = 1.0_dp
                szi = 0.0_dp
              end if

              if ( ( abs ( e ) + ei ) * bn <= an ) then
                a1 = cz * b11 + szr * b12
                a1i = szi * b12
                a2 = szr * b22
                a2i = szi * b22
              else
                a1 = cz * a11 + szr * a12
                a1i = szi * a12
                a2 = cz * a21 + szr * a22
                a2i = szi * a22
              end if
    !
    !  Choose complex Q.
    !
              cq = sqrt ( a1 * a1 + a1i * a1i )

              if ( cq /= 0.0_dp ) then
                sqr = ( a1 * a2 + a1i * a2i ) / cq
                sqi = ( a1 * a2i - a1i * a2 ) / cq
                r = sqrt ( cq * cq + sqr * sqr + sqi * sqi )
                cq = cq / r
                sqr = sqr / r
                sqi = sqi / r
              else
                sqr = 1.0_dp
                sqi = 0.0_dp
              end if
    !
    !  Compute diagonal elements that would result if transformations were applied.
    !
              ssr = sqr * szr + sqi * szi
              ssi = sqr * szi - sqi * szr
              i = 1
              tr = cq * cz * a11 + cq * szr * a12 + sqr * cz * a21 + ssr * a22
              ti = cq * szi * a12 - sqi * cz * a21 + ssi * a22
              dr = cq * cz * b11 + cq * szr * b12 + ssr * b22
              di = cq * szi * b12 + ssi * b22

              do

                t = ti * dr - tr * di

                if ( t < 0.0_dp ) then
                  j = en
                else
                  j = na
                end if

                r = sqrt ( dr * dr + di * di )
                beta(j) = bn * r
                alfr(j) = an * ( tr * dr + ti * di ) / r
                alfi(j) = an * t / r

                if ( i /= 1 ) then
                  exit
                end if

                i = 2
                tr =   ssr * a11 - sqr * cz * a12 - cq * szr * a21 + cq * cz * a22
                ti = - ssi * a11 - sqi * cz * a12 + cq * szi * a21
                dr =   ssr * b11 - sqr * cz * b12 + cq * cz * b22
                di = - ssi * b11 - sqi * cz * b12

              end do

              isw = 3 - isw

              cycle
            end if
    !
    !  Two real roots.
    !  Zero both A(EN,NA) and B(EN,NA).
    !
            e = e + ( c + sign ( sqrt ( d ), c ) )
            a11 = a11 - e * b11
            a12 = a12 - e * b12
            a22 = a22 - e * b22

            if ( abs ( a21 ) + abs ( a22 ) <= abs ( a11 ) + abs ( a12 ) ) then
              a1 = a12
              a2 = a11
            else
              a1 = a22
              a2 = a21
            end if

          end if
    !
    !  Choose and apply real Z.
    !
          s = abs ( a1 ) + abs ( a2 )
          u1 = a1 / s
          u2 = a2 / s
          r = sign ( sqrt ( u1 * u1 + u2 * u2 ), u1 )
          v1 = - ( u1 + r ) / r
          v2 = - u2 / r
          u2 = v2 / v1

          do i = 1, en
            t = a(i,en) + u2 * a(i,na)
            a(i,en) = a(i,en) + t * v1
            a(i,na) = a(i,na) + t * v2
            t = b(i,en) + u2 * b(i,na)
            b(i,en) = b(i,en) + t * v1
            b(i,na) = b(i,na) + t * v2
          end do

          if ( matz ) then

            do i = 1, n
              t = z(i,en) + u2 * z(i,na)
              z(i,en) = z(i,en) + t * v1
              z(i,na) = z(i,na) + t * v2
            end do

          end if

          if ( bn == 0.0_dp ) then
            a(en,na) = 0.0_dp
            b(en,na) = 0.0_dp
            alfr(na) = a(na,na)
            alfr(en) = a(en,en)
            if ( b(na,na) < 0.0_dp ) then
              alfr(na) = - alfr(na)
            end if

            if ( b(en,en) < 0.0_dp ) then
              alfr(en) = - alfr(en)
            end if

            beta(na) = abs ( b(na,na) )
            beta(en) = abs ( b(en,en) )
            alfi(en) = 0.0_dp
            alfi(na) = 0.0_dp
            isw = 3 - isw
            cycle
          end if

          if ( abs ( e ) * bn <= an ) then
            a1 = b(na,na)
            a2 = b(en,na)
          else
            a1 = a(na,na)
            a2 = a(en,na)
          end if

        end if
    !
    !  Choose and apply real Q.
    !
        s = abs ( a1 ) + abs ( a2 )

        if ( s /= 0.0_dp ) then

          u1 = a1 / s
          u2 = a2 / s
          r = sign ( sqrt ( u1 * u1 + u2 * u2 ), u1 )
          v1 = - ( u1 + r ) / r
          v2 = - u2 / r
          u2 = v2 / v1

          do j = na, n
            t = a(na,j) + u2 * a(en,j)
            a(na,j) = a(na,j) + t * v1
            a(en,j) = a(en,j) + t * v2
            t = b(na,j) + u2 * b(en,j)
            b(na,j) = b(na,j) + t * v1
            b(en,j) = b(en,j) + t * v2
          end do

        end if

        a(en,na) = 0.0_dp
        b(en,na) = 0.0_dp
        alfr(na) = a(na,na)
        alfr(en) = a(en,en)
        if ( b(na,na) < 0.0_dp ) then
          alfr(na) = - alfr(na)
        end if

        if ( b(en,en) < 0.0_dp ) then
          alfr(en) = - alfr(en)
        end if

        beta(na) = abs ( b(na,na) )
        beta(en) = abs ( b(en,en) )
        alfi(en) = 0.0_dp
        alfi(na) = 0.0_dp
        isw = 3 - isw

      end do

      b(n,1) = epsb

    End Subroutine qzval



    !    Subroutine qzvec ( n, a, b, alfr, alfi, beta, z )
    !
    !! QZVEC computes eigenvectors for a generalized eigenvalue problem.
    !
    !  Discussion:
    !
    !    QZVEC is the optional fourth step of the QZ algorithm
    !    for solving generalized matrix eigenvalue problems.
    !
    !    QZVEC accepts a pair of real matrices, one of them in
    !    quasi-triangular form (in which each 2-by-2 block corresponds to
    !    a pair of complex eigenvalues) and the other in upper triangular
    !    form.  It computes the eigenvectors of the triangular problem and
    !    transforms the results back to the original coordinate system.
    !    it is usually preceded by QZHES, QZIT, and QZVAL.
    !
    !  Modified:
    !
    !    25 February 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrices.
    !
    !    Input, real(kind=dp) A(N,N), contains a real upper quasi-triangular
    !    matrix.  Its subdiagonal elements provide information about the storage of
    !    the complex eigenvectors.
    !
    !    Input/output, real(kind=dp) B(N,N).  On input, a real upper triangular
    !    matrix.  In addition, location B(N,1) contains the tolerance quantity EPSB
    !    computed and saved in QZIT.  On output, B has been destroyed.
    !
    !    Input, real(kind=dp) ALFR(N), ALFI(N), BETA(N), vectors whose ratios
    !      ( ALFR + I * ALFI ) / BETA
    !    are the generalized eigenvalues.  They are usually obtained from QZVAL.
    !
    !    Input/output, real(kind=dp) Z(N,N).  On input, the transformation
    !    matrix produced in the reductions by QZHES, QZIT, and QZVAL, if performed.
    !    If the eigenvectors of the triangular problem are desired, Z must contain
    !    the identity matrix.  On output, Z contains the real and imaginary parts of
    !    the eigenvectors:
    !    If ALFI(I) == 0.0, the I-th eigenvalue is real and the I-th column of Z
    !    contains its eigenvector.
    !    If 0.0 < ALFI(I), the eigenvalue is the first of a complex pair and the
    !    I-th and (I+1)-th columns of Z contain its eigenvector.
    !    If ALFI(I) < 0.0, the eigenvalue is the second of a complex pair and the
    !    (I-1)-th and I-th columns of Z contain the conjugate of its eigenvector.
    !    Each eigenvector is normalized so that the modulus of its largest
    !    component is 1.0_dp .
    !
    Subroutine qzvec ( n, a, b, alfr, alfi, beta, z )
      integer,                       intent(in)     :: n
      real(kind=dp), dimension(n,n), intent(in)     :: a
      real(kind=dp), dimension(n,n), intent(in out) :: b
      real(kind=dp), dimension(n),   intent(in)     :: alfr, alfi, beta
      real(kind=dp), dimension(n,n), intent(in out) :: z


      real(kind=dp) :: alfm,almi,almr,betm,d,di,dr,epsb,q,r,ra,rr,s,sa,t,t1,t2, &
                       ti,tr,w,w1,x,x1,y,z1,zz
      integer       :: en,enm2,i,isw,j,k,m,na

      epsb = b(n,1)
      isw = 1

      do en = n, 1, -1

        na = en - 1

        if ( isw == 2 ) then
          isw = 3 - isw
          go to 800
        end if

        if ( alfi(en) /= 0.0_dp ) then
          go to 710
        end if
    !
    !  Real vector.
    !
        m = en
        b(en,en) = 1.0_dp

        if ( na == 0 ) then
          go to 800
        end if

        alfm = alfr(m)
        betm = beta(m)

        do i = en - 1, en - na, -1

          w = betm * a(i,i) - alfm * b(i,i)

          r = 0.0_dp
          do j = m, en
            r = r + ( betm * a(i,j) - alfm * b(i,j) ) * b(j,en)
          end do

          if ( 1 < i .and. isw /= 2 ) then

            if ( betm * a(i,i-1) /= 0.0_dp ) then
              zz = w
              s = r
              isw = 3 - isw
              cycle
            end if

          end if

          m = i
    !
    !  Real 1-by-1 block.
    !
          if ( isw /= 2 ) then

            if ( w == 0.0_dp ) then
              t = epsb
            else
              t = w
            end if

            b(i,en) = - r / t
            cycle

          end if
    !
    !  Real 2-by-2 block.
    !
          x = betm * a(i,i+1) - alfm * b(i,i+1)
          y = betm * a(i+1,i)
          q = w * zz - x * y
          t = ( x * s - zz * r ) / q
          b(i,en) = t

          if ( abs ( x ) <= abs ( zz ) ) then
            b(i+1,en) = ( - s - y * t ) / zz
          else
            b(i+1,en) = ( - r - w * t ) / x
          end if

          isw = 3 - isw

        end do
    !
    !  End real vector.
    !
        go to 800
    !
    !  Complex vector.
    !
    710 continue

        m = na
        almr = alfr(m)
        almi = alfi(m)
        betm = beta(m)
    !
    !  Last vector component chosen imaginary so eigenvector matrix is triangular.
    !
        y = betm * a(en,na)
        b(na,na) = - almi * b(en,en) / y
        b(na,en) = ( almr * b(en,en) - betm * a(en,en) ) / y
        b(en,na) = 0.0_dp
        b(en,en) = 1.0_dp
        enm2 = na - 1

        do i = na - 1, na - enm2, -1

          w = betm * a(i,i) - almr * b(i,i)
          w1 = - almi * b(i,i)
          ra = 0.0_dp
          sa = 0.0_dp

          do j = m, en
            x = betm * a(i,j) - almr * b(i,j)
            x1 = - almi * b(i,j)
            ra = ra + x * b(j,na) - x1 * b(j,en)
            sa = sa + x * b(j,en) + x1 * b(j,na)
          end do

          if ( i /= 1 .and. isw /= 2 ) then

            if ( betm * a(i,i-1) /= 0.0_dp ) then
              zz = w
              z1 = w1
              r = ra
              s = sa
              isw = 2
              cycle
            end if

          end if

          m = i
    !
    !  Complex 2-by-2 block.
    !
          if ( isw == 2 ) then

            x = betm * a(i,i+1) - almr * b(i,i+1)
            x1 = - almi * b(i,i+1)
            y = betm * a(i+1,i)
            tr = y * ra - w * r + w1 * s
            ti = y * sa - w * s - w1 * r
            dr = w * zz - w1 * z1 - x * y
            di = w * z1 + w1 * zz - x1 * y
            if ( dr == 0.0_dp .and. di == 0.0_dp ) then
              dr = epsb
            end if
    !
    !  Complex 1-by-1 block.
    !
          else

            tr = - ra
            ti = - sa
            dr = w
            di = w1

          end if
    !
    !  Complex divide (t1,t2) = (tr,ti) / (dr,di),
    !
    775   continue

          if ( abs ( di ) <= abs ( dr ) ) then
            rr = di / dr
            d = dr + di * rr
            t1 = ( tr + ti * rr ) / d
            t2 = ( ti - tr * rr ) / d
          else
            rr = dr / di
            d = dr * rr + di
            t1 = ( tr * rr + ti ) / d
            t2 = ( ti * rr - tr ) / d
          end if

          if ( isw == 2 ) then

            b(i+1,na) = t1
            b(i+1,en) = t2
            isw = 1

            if ( abs ( y ) <= abs ( w ) + abs ( w1 ) ) then
              tr = - ra - x * b(i+1,na) + x1 * b(i+1,en)
              ti = - sa - x * b(i+1,en) - x1 * b(i+1,na)
              dr = w
              di = w1
              go to 775
            end if

            t1 = ( - r - zz * b(i+1,na) + z1 * b(i+1,en) ) / y
            t2 = ( - s - zz * b(i+1,en) - z1 * b(i+1,na) ) / y

          end if

          b(i,na) = t1
          b(i,en) = t2

        end do
    !
    !  End complex vector.
    !
        isw = 3 - isw

    800 continue

      end do
    !
    !  End back substitution.
    !  Transform to original coordinate system.
    !
      do j = n, 1, -1
        do i = 1, n
          zz = 0.0_dp
          do k = 1, j
            zz = zz + z(i,k) * b(k,j)
          end do
          z(i,j) = zz
        end do
      end do
    !
    !  Normalize so that modulus of largest component of each vector is 1.
    !  (ISW is 1 initially from before).
    !
      do j = 1, n

        d = 0.0_dp

        if ( isw /= 2 ) then

          if ( alfi(j) /= 0.0_dp ) then
            isw = 3 - isw
            cycle
          end if

          do i = 1, n
            d = max ( d, abs ( z(i,j) ) )
          end do

          z(1:n,j) = z(1:n,j) / d

          cycle

        else

          do i = 1, n
            r = abs ( z(i,j-1) ) + abs ( z(i,j) )
            if ( r /= 0.0_dp ) then
              r = r * sqrt ( ( z(i,j-1) / r )**2 + ( z(i,j) / r )**2 )
            end if
            d = max ( d, r )
          end do

          z(1:n,j-1) = z(1:n,j-1) / d
          z(1:n,j) = z(1:n,j) / d

          isw = 3 - isw

        end if

      end do

    End Subroutine qzvec



    !     Subroutine r8mat_identity ( n, a )
    !
    !! R8MAT_IDENTITY stores the identity matrix in an R8MAT.
    !
    !  Discussion:
    !
    !    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
    !
    !  Modified:
    !
    !    24 March 2000
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of A.
    !
    !    Output, real(kind=dp) A(N,N), the N by N identity matrix.
    !
    Subroutine r8mat_identity ( n, a )
      integer, intent(in) :: n

      real(kind=dp), dimension(n,n), intent(out) :: a
      integer ::i

      a(1:n,1:n) = 0.0_dp
      do i = 1, n
        a(i,i) = 1.0_dp
      end do

    End Subroutine r8mat_identity

    !    Subroutine r8mat_print ( m, n, a, title )
    !
    !! R8MAT_PRINT prints an R8MAT.
    !
    !  Discussion:
    !
    !    An R8MAT is an array of R8 values.
    !
    !  Modified:
    !
    !    12 September 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Arguments:
    !
    !    Input, integer ::M, the number of rows in A.
    !
    !    Input, integer ::N, the number of columns in A.
    !
    !    Input, real(kind=dp) A(M,N), the matrix.
    !
    !    Input, character ( len = * ) TITLE, a title.
    !
    Subroutine r8mat_print ( m, n, a, title )
      integer,                       intent(in) :: m, n
      real(kind=dp), dimension(m,n), intent(in) :: a
      character ( len = * ),         intent(in) :: title

      call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

    End Subroutine r8mat_print



    !    Subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )
    !
    !! R8MAT_PRINT_SOME prints some of an R8MAT.
    !
    !  Discussion:
    !
    !    An R8MAT is an array of R8 values.
    !
    !  Modified:
    !
    !    26 March 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Arguments:
    !
    !    Input, integer ::M, N, the number of rows and columns.
    !
    !    Input, real(kind=dp) A(M,N), an M by N matrix to be printed.
    !
    !    Input, integer ::ILO, JLO, the first row and column to print.
    !
    !    Input, integer ::IHI, JHI, the last row and column to print.
    !
    !    Input, character ( len = * ) TITLE, a title.
    !
    Subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )
      integer,                       intent(in) :: m, n, ilo, jlo, ihi, jhi
      real(kind=dp), dimension(m,n), intent(in) :: a
      character ( len = * ),         intent(in) :: title

      integer, parameter :: incx = 5
      character ( len = 14 ), dimension(incx) :: ctemp
      integer :: i,i2hi,i2lo,inc,j,j2,j2hi,j2lo

      write(*, '(a)' ) ''
      write(*, '(a)' ) trim ( title )

      do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write(*, '(a)' ) ''

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write(ctemp(j2), '(i8,6x)' ) j
        end do

        write(*, '(''  Col   '',5a14)' ) ctemp(1:inc)
        write(*, '(a)' ) '  Row'
        write(*, '(a)' ) ''

        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi

          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
              write(ctemp(j2), '(f8.0,6x)' ) a(i,j)
            else
              write(ctemp(j2), '(g14.6)' ) a(i,j)
            end if

          end do

          write(*, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

        end do

      end do

    End Subroutine r8mat_print_some



    !    Subroutine r8mat_uniform_01 ( m, n, seed, r )
    !
    !! R8MAT_UNIFORM_01 fills an R8MAT with unit pseudorandom numbers.
    !
    !  Discussion:
    !
    !    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
    !
    !  Modified:
    !
    !    11 August 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Paul Bratley, Bennett Fox, Linus Schrage,
    !    A Guide to Simulation,
    !    Springer Verlag, pages 201-202, 1983.
    !
    !    Bennett Fox,
    !    Algorithm 647:
    !    Implementation and Relative Efficiency of Quasirandom
    !    Sequence Generators,
    !    ACM Transactions on Mathematical Software,
    !    Volume 12, Number 4, pages 362-376, 1986.
    !
    !    Peter Lewis, Allen Goodman, James Miller,
    !    A Pseudo-Random Number Generator for the System/360,
    !    IBM Systems Journal,
    !    Volume 8, pages 136-143, 1969.
    !
    !  Arguments:
    !
    !    Input, integer ::M, N, the number of rows and columns in
    !    the array.
    !
    !    Input/output, integer ::SEED, the "seed" value, which
    !    should NOT be 0.  On output, SEED has been updated.
    !
    !    Output, real(kind=dp) R(M,N), the array of pseudorandom values.
    !
    Subroutine r8mat_uniform_01 ( m, n, seed, r )
      integer,                       intent(in)     :: m, n
      integer,                       intent(in out) :: seed
      real(kind=dp), dimension(m,n), intent(out)    :: r

      integer, parameter :: i4_huge = 2147483647
      integer ::i,j,k


      do j = 1, n

        do i = 1, m

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed < 0 ) then
            seed = seed + i4_huge
          end if

          r(i,j) = real ( seed, kind = 8 ) * 4.656612875D-10

        end do
      end do

    End Subroutine r8mat_uniform_01



    !    Subroutine r8vec_print ( n, a, title )
    !
    !! R8VEC_PRINT prints an R8VEC.
    !
    !  Discussion:
    !
    !    An R8VEC is a vector of R8 values.
    !
    !  Modified:
    !
    !    22 August 2000
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Arguments:
    !
    !    Input, integer ::N, the number of components of the vector.
    !
    !    Input, real(kind=dp) A(N), the vector to be printed.
    !
    !    Input, character ( len = * ) TITLE, a title.
    !
    Subroutine r8vec_print ( n, a, title )
      integer,                       intent(in) :: n
      real(kind=dp), dimension(n),   intent(in) :: a
      character ( len = * ),         intent(in) :: title

      integer ::i

      write(*, '(a)' ) ''
      write(*, '(a)' ) trim ( title )
      write(*, '(a)' ) ''
      do i = 1, n
        write(*, '(2x,i8,2x,g16.8)' ) i, a(i)
      end do

    End Subroutine r8vec_print



    !    Subroutine r8vec2_print ( n, a1, a2, title )
    !
    !! R8VEC2_PRINT prints an R8VEC2.
    !
    !  Discussion:
    !
    !    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
    !    as two separate vectors A1 and A2.
    !
    !  Modified:
    !
    !    13 December 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Arguments:
    !
    !    Input, integer ::N, the number of components of the vector.
    !
    !    Input, real(kind=dp) A1(N), A2(N), the vectors to be printed.
    !
    !    Input, character ( len = * ) TITLE, a title.
    !
    Subroutine r8vec2_print ( n, a1, a2, title )
      integer,                       intent(in) :: n
      real(kind=dp), dimension(n),   intent(in) :: a1,a2
      character ( len = * ),         intent(in) :: title

      integer ::i

      write(*, '(a)' ) ''
      write(*, '(a)' ) trim ( title )
      write(*, '(a)' ) ''

      if ( all ( a1(1:n) == aint ( a1(1:n) ) ) .and. &
           all ( a2(1:n) == aint ( a2(1:n) ) ) ) then
        do i = 1, n
          write(*, '(i8,2i8)' ) i, int ( a1(i) ), int ( a2(i) )
        end do
      else if ( all ( abs ( a1(1:n) ) < 1000000.0_dp ) .and. &
                all ( abs ( a2(1:n) ) < 1000000.0_dp ) ) then
        do i = 1, n
          write(*, '(i8,2f14.6)' ) i, a1(i), a2(i)
        end do
      else
        do i = 1, n
          write(*, '(i8,2g14.6)' ) i, a1(i), a2(i)
        end do
      end if

    End Subroutine r8vec2_print



    !    Subroutine ratqr ( n, eps1, d, e, e2, m, w, ind, bd, typ, idef, ierr )
    !
    !! RATQR computes selected eigenvalues of a real symmetric tridiagonal matrix.
    !
    !  Discussion:
    !
    !    RATQR finds the algebraically smallest or largest eigenvalues of a
    !    symmetric tridiagonal matrix by the rational QR method with Newton
    !    corrections.
    !
    !  Modified:
    !
    !    03 February 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input/output, real(kind=dp) EPS1.  On input, a theoretical absolute
    !    error tolerance for the computed eigenvalues.  If the input EPS1 is
    !    non-positive, or indeed smaller than its default value, it is reset at
    !    each iteration to the respective default value, namely, the product of
    !    the relative machine precision and the magnitude of the current eigenvalue
    !    iterate.  The theoretical absolute error in the K-th eigenvalue is usually
    !    not greater than K times EPS1.  On output, EPS1 is unaltered unless it has
    !    been reset to its (last) default value.
    !
    !    Input, real(kind=dp) D(N), the diagonal elements of the input matrix.
    !
    !    Input, real(kind=dp) E(N), the subdiagonal elements of the input matrix
    !    in E(2:N).  E(1) is arbitrary.
    !
    !    Input/output, real(kind=dp) E2(N).  On input, E2(2:N-1) contains the
    !    squares of the corresponding elements of E, and E2(1) is arbitrary.  On
    !    output, elements of E2 corresponding to elements of E regarded as
    !    negligible have been replaced by zero, causing the matrix to split into
    !    a direct sum of submatrices.  E2(1) is set to 0.0_dp if the smallest
    !    eigenvalues have been found, and to 2.0_dp if the largest eigenvalues
    !    have been found.  E2 is otherwise unaltered (unless overwritten by BD).
    !
    !    Input, integer ::M, the number of eigenvalues to be found.
    !
    !    Output, real(kind=dp) W(M), the M algebraically smallest eigenvalues in
    !    ascending order, or the M largest eigenvalues in descending order.
    !    If an error exit is made because of an incorrect specification of IDEF,
    !    no eigenvalues are found.  If the Newton iterates for a particular
    !    eigenvalue are not monotone, the best estimate obtained is returned
    !    and IERR is set.  W may coincide with D.
    !
    !    Output, integer IND(N), contains in its first M positions the submatrix
    !    indices associated with the corresponding eigenvalues in W:
    !    1 for eigenvalues belonging to the first submatrix from the top, 2 for
    !    those belonging to the second submatrix, and so on.
    !
    !    Output, real(kind=dp) BD(N), contains refined bounds for the
    !    theoretical errors of the corresponding eigenvalues in W.  These bounds
    !    are usually within the tolerance specified by EPS1.  BD may coincide
    !    with E2.
    !
    !    Input, logical typ, should be set to TRUE if the smallest eigenvalues
    !    are to be found, and to FALSE if the largest eigenvalues are to be found.
    !
    !    Input, integer ::IDEF, should be set to 1 if the input matrix
    !    is known to be positive definite, to -1 if the input matrix is known to
    !    be negative  definite, and to 0 otherwise.
    !
    !    Output, integer ::IERR, error flag.
    !    0, for normal return,
    !    6*N+1, if IDEF is set to 1 and typ to .true. when the matrix is not
    !      positive definite, or if IDEF is set to -1 and typ to .false.
    !      when the matrix is not negative definite,
    !    5*N+K, if successive iterates to the K-th eigenvalue are not monotone
    !      increasing, where K refers to the last such occurrence.
    !
    Subroutine ratqr ( n, eps1, d, e, e2, m, w, ind, bd, typ, idef, ierr )
      integer,                       intent(in)     :: n
      real(kind=dp),                 intent(in out) :: eps1
      real(kind=dp), dimension(n),   intent(in)     :: d,e
      real(kind=dp), dimension(n),   intent(in out) :: e2
      integer,                       intent(in)     :: m
      real(kind=dp), dimension(m),   intent(out)    :: w
      integer,       dimension(n),   intent(out)    :: ind
      real(kind=dp), dimension(n),   intent(out)    :: bd
      logical,                       intent(in)     :: typ
      integer,                       intent(in)     :: idef
      integer,                       intent(out)    :: ierr

      real(kind=dp) :: delta ,ep,errv,f,p,q,qp,r,s,tot

      integer :: i,ii,j,jdef,k
      logical :: irreg

      ierr = 0
      jdef = idef
      w(1:n) = d(1:n)

      if ( .not. typ ) then
        j = 1
        w(1:n) = - w(1:n)
        jdef = - jdef
      end if

      errv = 0.0_dp
      s = 0.0_dp
    !
    !  Look for small sub-diagonal entries and define initial shift
    !  from lower Gerschgorin bound.
    !
    !  Copy E2 array into BD.
    !
      tot = w(1)
      q = 0.0_dp
      j = 0

      do i = 1, n

        p = q

        if ( 1 == i ) then
          e2(i) = 0.0_dp
        else if ( p <= ( abs ( d(i) ) + abs (  d(i-1) ) ) * epsilon ( p ) ) then
          e2(i) = 0.0_dp
        end if

        bd(i) = e2(i)
    !
    !  Count also if element of E2 has underflowed.
    !
        if ( e2(i) == 0.0_dp ) then
          j = j + 1
        end if

        ind(i) = j
        q = 0.0_dp
        if ( i < n ) then
          q = abs ( e(i+1) )
        end if

        tot = min ( w(i) - p - q, tot )

      end do

      if ( jdef == 1 .and. tot < 0.0_dp ) then
        tot = 0.0_dp
      else
        w(1:n) = w(1:n) - tot
      end if

      do k = 1, m
    !
    !  Next QR transformation.
    !
        irreg = .true.

        do

          tot = tot + s
          delta = w(n) - s
          i = n
          f = abs ( tot ) * epsilon ( f )
          if ( eps1 < f ) then
            eps1 = f
          end if

          if ( delta <= eps1 ) then

            if ( delta < - eps1 ) then
              ierr = 6 * n + 1
              return
            end if

            irreg = .false.
            exit

          end if
    !
    !  Replace small sub-diagonal squares by zero to reduce the incidence of
    !  underflows.
    !
          do j = k + 1, n
            if ( bd(j) <= epsilon ( bd(j) ) ** 2 ) then
              bd(j) = 0.0_dp
            end if
          end do

          f = bd(n) / delta
          qp = delta + f
          p = 1.0_dp

          do i = n - 1, k, -1
            q = w(i) - s - f
            r = q / qp
            p = p * r + 1.0_dp
            ep = f * r
            w(i+1) = qp + ep
            delta = q - ep

            if ( delta <= eps1 ) then

              if ( delta < - eps1 ) then
                ierr = 6 * n + 1
                return
              end if

              irreg = .false.
              exit

            end if

            f = bd(i) / q
            qp = delta + f
            bd(i+1) = qp * ep

          end do

          if ( .not. irreg ) then
            exit
          end if

          w(k) = qp
          s = qp / p

          if ( tot + s <= tot ) then
            exit
          end if

        end do
    !
    !  Set error: irregular end of iteration.
    !  Deflate minimum diagonal element.
    !
        if ( irreg ) then

          ierr = 5 * n + k
          s = 0.0_dp
          delta = qp

          do j = k, n
            if ( w(j) <= delta ) then
              i = j
              delta = w(j)
            end if
          end do

        end if
    !
    !  Convergence.
    !
        if ( i < n ) then
          bd(i+1) = bd(i) * f / qp
        end if

        ii = ind(i)

        do j = i - 1, k, -1
          w(j+1) = w(j) - s
          bd(j+1) = bd(j)
          ind(j+1) = ind(j)
        end do

        w(k) = tot
        errv = errv + abs ( delta)
        bd(k) = errv
        ind(k) = ii

      end do

      if ( .not. typ ) then
        f = bd(1)
        e2(1) = 2.0_dp
        bd(1) = f
        w(1:n) = - w(1:n)
        jdef = - jdef
      end if

    End Subroutine ratqr



    !    Subroutine rebak ( n, b, dl, m, z )
    !
    !! REBAK determines eigenvectors by undoing the REDUC transformation.
    !
    !  Discussion:
    !
    !    REBAK forms the eigenvectors of a generalized
    !    symmetric eigensystem by back transforming those of the
    !    derived symmetric matrix determined by REDUC.
    !
    !  Modified:
    !
    !    04 March 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input, real(kind=dp) B(N,N), contains information about the similarity
    !    transformation (Cholesky decomposition) used in the reduction by REDUC
    !    in its strict lower triangle.
    !
    !    Input, real(kind=dp) DL(N), further information about the
    !    transformation.
    !
    !    Input, integer ::M, the number of eigenvectors to be back
    !    transformed.
    !
    !    Input/output, real(kind=dp) Z(N,M).  On input, the eigenvectors to be
    !    back transformed in its first M columns.  On output, the transformed
    !    eigenvectors.
    !
    Subroutine rebak ( n, b, dl, m, z )
      integer,                       intent(in)     :: n
      real(kind=dp), dimension(n,n), intent(in)     :: b
      real(kind=dp), dimension(n),   intent(in)     :: dl
      integer,                       intent(in)     :: m
      real(kind=dp), dimension(n,m), intent(in out) :: z

      integer ::i,j

      do j = 1, m
        do i = n, 1, -1
          z(i,j) = ( z(i,j) - dot_product ( b(i+1:n,i), z(i+1:n,j) ) ) / dl(i)
        end do
      end do

    End Subroutine rebak



    !    Subroutine rebakb ( n, b, dl, m, z )
    !
    !! REBAKB determines eigenvectors by undoing the REDUC2 transformation.
    !
    !  Discussion:
    !
    !    REBAKB forms the eigenvectors of a generalized symmetric eigensystem by
    !    back transforming those of the derived symmetric matrix determined
    !    by REDUC2.
    !
    !  Modified:
    !
    !    04 March 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input, real(kind=dp) B(N,N), contains information about the similarity
    !    transformation (Cholesky decomposition) used in the reduction by REDUC2
    !    in its strict lower triangle.
    !
    !    Input, real(kind=dp) DL(N), further information about the
    !    transformation.
    !
    !    Input, integer ::M, the number of eigenvectors to be back
    !    transformed.
    !
    !    Input/output, real(kind=dp) Z(N,M).  On input, the eigenvectors to be
    !    back transformed in its first M columns.  On output, the transformed
    !    eigenvectors.
    !
    Subroutine rebakb ( n, b, dl, m, z )
      integer,                       intent(in)     :: n
      real(kind=dp), dimension(n,n), intent(in)     :: b
      real(kind=dp), dimension(n),   intent(in)     :: dl
      integer,                       intent(in)     :: m
      real(kind=dp), dimension(n,m), intent(in out) :: z

      integer ::i,j

      do j = 1, m
        do i = n, 1, -1
          z(i,j) = dl(i) * z(i,j) + dot_product ( b(i,1:i-1), z(1:i-1,j) )
        end do
      end do

    End Subroutine rebakb



    !    Subroutine reduc ( n, a, b, dl, ierr )
    !
    !! REDUC reduces the eigenvalue problem A*x=lambda*B*x to A*x=lambda*x.
    !
    !  Discussion:
    !
    !    REDUC reduces the generalized symmetric eigenproblem
    !    a x=(lambda) b x, where B is positive definite, to the standard
    !    symmetric eigenproblem using the Cholesky factorization of B.
    !
    !  Modified:
    !
    !    04 March 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrices A and B.  If the
    !    Cholesky factor L of B is already available, N should be prefixed with a
    !    minus sign.
    !
    !    Input/output, real(kind=dp) A(N,N).  On input, A contains a real
    !    symmetric matrix.  Only the full upper triangle of the matrix need be
    !    supplied.  On output, A contains in its full lower triangle the full lower
    !    triangle of the symmetric matrix derived from the reduction to the
    !    standard form.  The strict upper triangle of a is unaltered.
    !
    !    Input/output, real(kind=dp) B(N,N).  On input, the real symmetric
    !    input matrix.  Only the full upper triangle of the matrix need be supplied.
    !    If N is negative, the strict lower triangle of B contains, instead, the
    !    strict lower triangle of its Cholesky factor L.  In any case, on output,
    !    B contains in its strict lower triangle the strict lower triangle of
    !    its Cholesky factor L.  The full upper triangle of B is unaltered.
    !
    !    Input/output, real(kind=dp) DL(N).  If N is negative, then the DL
    !    contains the diagonal elements of L on input.  In any case, DL will contain
    !    the diagonal elements of L on output,
    !
    !    Output, integer ::IERR, error flag.
    !    0, for normal return,
    !    7*N+1, if B is not positive definite.
    !
    Subroutine reduc ( n, a, b, dl, ierr )
      integer,                       intent(in)     :: n
      real(kind=dp), dimension(n,n), intent(in out) :: a,b
      real(kind=dp), dimension(n),   intent(in out) :: dl
      integer,                       intent(out)    :: ierr

      integer       :: i,j,k,nn
      real(kind=dp) :: x,y

      ierr = 0
      nn = abs ( n )
    !
    !  Form L in the arrays B and DL.
    !
      do i = 1, n

         do j = i, n

            x = b(i,j)

            do k = 1, i - 1
              x = x - b(i,k) * b(j,k)
            end do

            if ( j == i ) then

              if ( x <= 0.0_dp ) then
                ERR_EisPack=.true.
                ERR_EisPack_Mess=" REDUC - Fatal error!: The matrix is not positive definite.."
                ierr = 7 * n + 1
                return
              end if

              y = sqrt ( x )
              dl(i) = y
            else
              b(j,i) = x / y
            end if

        end do

      end do
    !
    !  Form the transpose of the upper triangle of INV(L)*A
    !  in the lower triangle of the array A.
    !
      do i = 1, nn

        y = dl(i)

        do j = i, nn

          x = a(i,j)

          do k = 1, i - 1
            x = x - b(i,k) * a(j,k)
          end do

          a(j,i) = x / y

        end do

      end do
    !
    !  Pre-multiply by INV(L) and overwrite.
    !
      do j = 1, nn

        do i = j, nn

          x = a(i,j)

          do k = j, i - 1
            x = x - a(k,j) * b(i,k)
          end do

          do k = 1, j - 1
            x = x - a(j,k) * b(i,k)
          end do

          a(i,j) = x / dl(i)

        end do

      end do

    End Subroutine reduc



    !    Subroutine reduc2 ( n, a, b, dl, ierr )
    !
    !! REDUC2 reduces the eigenvalue problem A*B*x=lamdba*x to A*x=lambda*x.
    !
    !  Discussion:
    !
    !    REDUC2 reduces the generalized symmetric eigenproblems
    !    a b x=(lambda) x or b a y=(lambda) y, where B is positive definite,
    !    to the standard symmetric eigenproblem using the Cholesky
    !    factorization of B.
    !
    !  Modified:
    !
    !    04 March 2018
    !
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrices A and B.  If the
    !    Cholesky factor L of B is already available, N should be prefixed with a
    !    minus sign.
    !
    !    Input/output, real(kind=dp) A(N,N).  On input, A contains a real
    !    symmetric matrix.  Only the full upper triangle of the matrix need be
    !    supplied.  On output, A contains in its full lower triangle the full lower
    !    triangle of the symmetric matrix derived from the reduction to the
    !    standard form.  The strict upper triangle of a is unaltered.
    !
    !    Input/output, real(kind=dp) B(N,N).  On input, the real symmetric
    !    input matrix.  Only the full upper triangle of the matrix need be supplied.
    !    If N is negative, the strict lower triangle of B contains, instead, the
    !    strict lower triangle of its Cholesky factor L.  In any case, on output,
    !    B contains in its strict lower triangle the strict lower triangle of
    !    its Cholesky factor L.  The full upper triangle of B is unaltered.
    !
    !    Input/output, real(kind=dp) DL(N).  If N is negative, then the DL
    !    contains the diagonal elements of L on input.  In any case, DL will contain
    !    the diagonal elements of L on output,
    !
    !    Output, integer ::IERR, error flag.
    !    0, for normal return,
    !    7*N+1, if B is not positive definite.
    !
    Subroutine reduc2 ( n, a, b, dl, ierr )
      integer,                       intent(in)     :: n
      real(kind=dp), dimension(n,n), intent(in out) :: a,b
      real(kind=dp), dimension(n),   intent(in out) :: dl
      integer,                       intent(out)    :: ierr

      integer       :: i,j,k,nn
      real(kind=dp) :: x,y

      ierr = 0
      nn = abs ( n )
    !
    !  Form L in the arrays B and DL.
    !
      do i = 1, n

         do j = i, n

            x = b(i,j)
            do k = 1, i - 1
              x = x - b(i,k) * b(j,k)
            end do

            if ( j == i ) then

              if ( x <= 0.0_dp ) then
                ERR_EisPack=.true.
                ERR_EisPack_Mess=" REDUC2 - Fatal error!: The matrix is not positive definite."
                ierr = 7 * n + 1
                return
              end if

              y = sqrt ( x )
              dl(i) = y

            else

              b(j,i) = x / y

            end if

        end do

      end do
    !
    !  Form the lower triangle of A*L in the lower triangle of A.
    !
      do i = 1, nn

        do j = 1, i

          x = a(j,i) * dl(j)
          do k = j + 1, i
            x = x + a(k,i) * b(k,j)
          end do
          do k = i + 1, nn
            x = x + a(i,k) * b(k,j)
          end do
          a(i,j) = x

        end do

      end do
    !
    !  Pre-multiply by L' and overwrite.
    !
      do i = 1, nn

        y = dl(i)

        do j = 1, i

          x = y * a(i,j)

          do k = i + 1, nn
            x = x + a(k,j) * b(k,i)
          end do

          a(i,j) = x

        end do

      end do

    End Subroutine reduc2



    !     Subroutine rg_elm ( n, a, wr, wi, matz, z, ierr )
    !
    !! RG_ELM computes eigenvalues and eigenvectors of a real general matrix.
    !
    !  Discussion:
    !
    !    RG_ELM calls EISPACK routines to find the eigenvalues and eigenvectors
    !    of a real general matrix, using elementary transformations,
    !
    !
    !  Modified:
    !
    !    09 February 2018
    !

    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input/output, real(kind=dp) A(N,N), the real general matrix.  On
    !    output, A has been overwritten.
    !
    !    Input, logical MATZ, is false if only eigenvalues are desired,
    !    and true if both eigenvalues and eigenvectors are desired.
    !
    !    Output, real(kind=dp) WR(N), WI(N), the real and imaginary parts,
    !    respectively, of the eigenvalues.  Complex conjugate pairs of eigenvalues
    !    appear consecutively with the eigenvalue having the positive imaginary
    !    part first.
    !
    !    Output, real(kind=dp) Z(N,N), contains the real and imaginary parts of
    !    the eigenvectors if MATZ is true.  If the J-th eigenvalue is real, the
    !    J-th column of Z contains its eigenvector.  If the J-th eigenvalue is
    !    complex with positive imaginary part, the J-th and (J+1)-th columns of
    !    Z contain the real and imaginary parts of its eigenvector.  The
    !    conjugate of this vector is the eigenvector for the conjugate eigenvalue.
    !
    !    Output, integer ::IERR, an error completion code described in
    !    the documentation for HQR and HQR2.  The normal completion code is zero.
    !
    Subroutine rg_elm ( n, a, wr, wi, matz, z, ierr )
      integer,                       intent(in)     :: n
      real(kind=dp), dimension(n,n), intent(in out) :: a
      real(kind=dp), dimension(n),   intent(out)    :: wr,wi
      logical,                       intent(in)     :: matz
      real(kind=dp), dimension(n,n), intent(out)    :: z
      integer,                       intent(out)    :: ierr

      real(kind=dp), dimension(n) :: fv1
      integer ::is1,is2
      integer, dimension(n) ::iv1

      call balanc ( n, a, is1, is2, fv1 )

      call elmhes ( n, is1, is2, a, iv1 )

      if ( .not. matz ) then

        call hqr ( n, is1, is2, a, wr, wi, ierr )

        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RG_ELM - Fatal error!: Error return from HQR"
          return
        end if

      else

        call eltran ( n, is1, is2, a, iv1, z )

        call hqr2 ( n, is1, is2, a, wr, wi, z, ierr )

        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RG_ELM - Fatal error!: Error return from HQR2"
          return
        end if

        call balbak ( n, is1, is2, fv1, n, z )

      end if

    End Subroutine rg_elm

    !!----
    !!----  Subroutine rg_ort( n,a,wr,wi,matz,z,ierr)
    !!----  integer,                         intent(in)    :: n
    !!----  real(kind = dp), dimension(n,n), intent(in out):: a
    !!----  real(kind = dp), dimension(n),   intent(out)   :: wi, wr
    !!----  logical,                         intent(in)    :: matz
    !!----  real(kind = dp), dimension(n,n), intent(out)   :: z
    !!----  integer,                         intent(out)   :: ierr
    !!----
    !!---- RG_ORT computes eigenvalues and eigenvectors of a real general matrix.
    !!----
    !!----  Discussion:
    !!----
    !!----    RG_ORT calls EISPACK routines to find the eigenvalues and eigenvectors
    !!----    of a real general matrix, using orthogonal transformations.
    !!----
    !!----  Modified:
    !!----
    !!----    11 February 2018, 29 November 2019 (JRC)
    !!----
    !!----  Author:
    !!----
    !!----    John Burkardt
    !!----
    !!----  Arguments:
    !!----
    !!----    Input, integer ::N, the order of the matrix.
    !!----
    !!----    Input/output, real(kind=dp) A(N,N), the real general matrix.  On
    !!----    output, A has been overwritten.
    !!----
    !!----    Input, logical MATZ, is false if only eigenvalues are desired,
    !!----    and true if both eigenvalues and eigenvectors are desired.
    !!----
    !!----    Output, real(kind=dp) WR(N), WI(N), the real and imaginary parts,
    !!----    respectively, of the eigenvalues.  Complex conjugate pairs of eigenvalues
    !!----    appear consecutively with the eigenvalue having the positive imaginary
    !!----    part first.
    !!----
    !!----    Output, real(kind=dp) Z(N,N), contains the real and imaginary parts of
    !!----    the eigenvectors if MATZ is true.  If the J-th eigenvalue is real, the
    !!----    J-th column of Z contains its eigenvector.  If the J-th eigenvalue is
    !!----    complex with positive imaginary part, the J-th and (J+1)-th columns of
    !!----    Z contain the real and imaginary parts of its eigenvector.  The
    !!----    conjugate of this vector is the eigenvector for the conjugate eigenvalue.
    !!----
    !!----    Output, integer ::IERR, an error completion code described in
    !!----    the documentation for HQR and HQR2.  The normal completion code is zero.
    !!----
    Subroutine rg_ort( n,a,wr,wi,matz,z,ierr)
      !
      integer,                         intent(in)    :: n
      real(kind = dp), dimension(n,n), intent(in out):: a
      real(kind = dp), dimension(n),   intent(out)   :: wi, wr
      logical,                         intent(in)    :: matz
      real(kind = dp), dimension(n,n), intent(out)   :: z
      integer,                         intent(out)   :: ierr

      real(kind = dp), dimension(n) :: fv1,ort
      integer:: is1,is2


      call balanc( n, a, is1, is2, fv1 )
      call orthes( n, is1, is2, a, ort )
      if ( .not. matz ) then
        call hqr( n, is1, is2, a, wr, wi, ierr )
        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RG_ORT: Error return from HQR"
          return
        end if
      else
        call ortran( n, is1, is2, a, ort, z )
        call hqr2( n, is1, is2, a, wr, wi, z, ierr )
        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RG_ORT: Error return from HQR2"
          return
        end if
        call balbak ( n, is1, is2, fv1, n, z )
      end if

    End Subroutine rg_ort



    !    Subroutine rgg ( n, a, b, alfr, alfi, beta, matz, z, ierr )
    !
    !! RGG: eigenvalues/vectors for the generalized problem A*x = lambda*B*x.
    !
    !  Discussion:
    !
    !    RGG calls the recommended sequence of EISPACK routines
    !    to find the eigenvalues and eigenvectors (if desired)
    !    for the real general generalized eigenproblem
    !
    !      A * x = lambda * B * x.
    !
    !  Modified:
    !
    !    09 February 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrices A and B.
    !
    !    Input/output, real(kind=dp) A(N,N), B(N,N), the two real general
    !    matrices.  On output, A and B have been overwritten.
    !
    !    Input, logical MATZ, is false if only eigenvalues are desired,
    !    and true if both eigenvalues and eigenvectors are desired.
    !
    !    Output, real(kind=dp) ALFR(N), ALFI(N), the real and imaginary parts,
    !    respectively, of the numerators of the eigenvalues.
    !
    !    Output, real(kind=dp) BETA(N), the denominators of the eigenvalues,
    !    which are thus given by the ratios (ALFR + I * ALFI ) / BETA.
    !    Complex conjugate pairs of eigenvalues appear consecutively
    !    with the eigenvalue having the positive imaginary part first.
    !
    !    Output, real(kind=dp) Z(N,N), contains the real and imaginary parts of
    !    the eigenvectors if MATZ is true.  If the J-th eigenvalue is real, the
    !    J-th column of Z contains its eigenvector.  If the J-th eigenvalue is
    !    complex with positive imaginary part, the J-th and (J+1)-th columns of
    !    Z contain the real and imaginary parts of its eigenvector.  The
    !    conjugate of this vector is the eigenvector for the conjugate eigenvalue.
    !
    !    Output, integer ::IERR, is set equal to an error completion code
    !    described in the documentation for QZIT.  The normal completion
    !    code is zero.
    !
    Subroutine rgg ( n, a, b, alfr, alfi, beta, matz, z, ierr )
      integer,                         intent(in)    :: n
      real(kind = dp), dimension(n,n), intent(in out):: a,b
      real(kind = dp), dimension(n),   intent(out)   :: alfr, alfi
      real(kind = dp), dimension(n),   intent(out)   :: beta
      logical,                         intent(in)    :: matz
      real(kind = dp), dimension(n,n), intent(out)   :: z
      integer,                         intent(out)   :: ierr

      real(kind=dp) :: eps1

      eps1 = 0.0_dp

      call qzhes ( n, a, b, matz, z )

      call qzit ( n, a, b, eps1, matz, z, ierr )

      if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RGG - Fatal error!: Error return from QZIT"
          return
      end if

      call qzval ( n, a, b, alfr, alfi, beta, matz, z )

      if ( matz ) then
        call qzvec ( n, a, b, alfr, alfi, beta, z )
      end if

    End Subroutine rgg



    !    Subroutine rs ( n, a, w, matz, z, ierr )
    !
    !! RS computes eigenvalues and eigenvectors of real symmetric matrix.
    !
    !  Discussion:
    !
    !    RS calls the recommended sequence of EISPACK routines
    !    to find the eigenvalues and eigenvectors (if desired)
    !    of a real symmetric matrix.
    !
    !
    !  Modified:
    !
    !    09 February 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input/output (?), real(kind=dp) A(N,N), the real symmetric matrix.
    !
    !    Input, logical MATZ, is false if only eigenvalues are desired,
    !    and true if both eigenvalues and eigenvectors are desired.
    !
    !    Output, real(kind=dp) W(N), the eigenvalues in ascending order.
    !
    !    Output, real(kind=dp) Z(N,N), contains the eigenvectors, if MATZ
    !    is true.
    !
    !    Output, integer ::IERR, is set equal to an error
    !    completion code described in the documentation for TQLRAT and TQL2.
    !    The normal completion code is zero.
    !
    Subroutine rs ( n, a, w, matz, z, ierr )
      integer,                         intent(in)    :: n
      real(kind = dp), dimension(n,n), intent(in out):: a
      real(kind = dp), dimension(n),   intent(out)   :: w
      logical,                         intent(in)    :: matz
      real(kind = dp), dimension(n,n), intent(out)   :: z
      integer,                         intent(out)   :: ierr

      real(kind=dp), dimension(n) :: fv1,fv2

      if ( .not. matz ) then

        call tred1 ( n, a, w, fv1, fv2 )

        call tqlrat ( n, w, fv2, ierr )

        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RS - Fatal error!: Error return from TQLRAT"
          return
        end if

      else

        call tred2 ( n, a, w, fv1, z )

        call tql2 ( n, w, fv1, z, ierr )

        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RS - Fatal error!: Error return from TQL2"
          return
       end if

      end if

    End Subroutine rs



    !    Subroutine rsb ( n, mb, a, w, matz, z, ierr )
    !
    !! RSB computes eigenvalues and eigenvectors of a real symmetric band matrix.
    !
    !  Discussion:
    !
    !    RSB calls the recommended sequence of EISPACK routines
    !    to find the eigenvalues and eigenvectors (if desired)
    !    of a real symmetric band matrix.
    !
    !  Modified:
    !
    !    09 February 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input, integer ::MB, the half band width of the matrix,
    !    defined as the number of adjacent diagonals, including the principal
    !    diagonal, required to specify the non-zero portion of the lower triangle
    !    of the matrix.
    !
    !    Input/Output (bandr modifies A), real(kind=dp) A(N,MB), contains the lower triangle
    !    of the real symmetric band matrix.  Its lowest subdiagonal is stored in the last N+1-MB
    !    positions of the first column, its next subdiagonal in the last
    !    N+2-MB positions of the second column, further subdiagonals similarly,
    !    and finally its principal diagonal in the N positions of the last
    !    column.  Contents of storages not part of the matrix are arbitrary.
    !
    !    Input, logical MATZ, is false if only eigenvalues are desired,
    !    and true if both eigenvalues and eigenvectors are desired.
    !
    !    Output, real(kind=dp) W(N), the eigenvalues in ascending order.
    !
    !    Output, real(kind=dp) Z(N,N), contains the eigenvectors, if MATZ
    !    is true.
    !
    !    Output, integer ::IERR, is set to an error
    !    completion code described in the documentation for TQLRAT and TQL2.
    !    The normal completion code is zero.
    !
    Subroutine rsb ( n, mb, a, w, matz, z, ierr )
      integer,                         intent(in)    :: n,mb
      real(kind = dp), dimension(n,mb),intent(in out):: a
      real(kind = dp), dimension(n),   intent(out)   :: w
      logical,                         intent(in)    :: matz
      real(kind = dp), dimension(n,n), intent(out)   :: z
      integer,                         intent(out)   :: ierr

      real(kind=dp), dimension(n) :: fv1,fv2

      if ( mb <= 0 ) then
        ierr = 12 * n
        ERR_EisPack=.true.
        ERR_EisPack_Mess=" RSB - Fatal error!: MB <= 0"
        return
      end if

      if ( n < mb ) then
        ierr = 12 * n
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RSB - Fatal error!: N < MB"
          return
      end if

      call bandr ( n, mb, a, w, fv1, fv2, matz, z )

      if ( .not. matz ) then

        call tqlrat ( n, w, fv2, ierr )

        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RSB - Fatal error!: Error return from TQLRAT"
          return
        end if

      else

        call tql2 ( n, w, fv1, z, ierr )

        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RSB - Fatal error!: Error return from TQL2"
          return
        end if

      end if

    End Subroutine rsb



    !    Subroutine rsg ( n, a, b, w, matz, z, ierr )
    !
    !! RSG computes eigenvalues/vectors, A*x=lambda*B*x, A symmetric, B pos-def.
    !
    !  Discussion:
    !
    !    RSG calls the recommended sequence of EISPACK routines
    !    to find the eigenvalues and eigenvectors (if desired)
    !    for the real symmetric generalized eigenproblem  ax = (lambda)bx.
    !
    !  Modified:
    !
    !    09 February 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrices A and B.
    !
    !    Input, real(kind=dp) A(N,N), contains a real symmetric matrix.
    !
    !    Input, real(kind=dp) B(N,N), contains a positive definite real
    !    symmetric matrix.
    !
    !    Input, logical MATZ, is false if only eigenvalues are desired,
    !    and true if both eigenvalues and eigenvectors are desired.
    !
    !    Output, real(kind=dp) W(N), the eigenvalues in ascending order.
    !
    !    Output, real(kind=dp) Z(N,N), contains the eigenvectors, if MATZ
    !    is true.
    !
    !    Output, integer ::IERR, is set to an error
    !    completion code described in the documentation for TQLRAT and TQL2.
    !    The normal completion code is zero.
    !
    Subroutine rsg ( n, a, b, w, matz, z, ierr )
      integer,                         intent(in)    :: n
      real(kind = dp), dimension(n,n),intent(in out) :: a,b
      real(kind = dp), dimension(n),   intent(out)   :: w
      logical,                         intent(in)    :: matz
      real(kind = dp), dimension(n,n), intent(out)   :: z
      integer,                         intent(out)   :: ierr

      real(kind=dp), dimension(n) :: fv1,fv2

      call reduc ( n, a, b, fv2, ierr )

      if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RSG - Fatal error!: Error return from REDUC: "//trim(ERR_EisPack_Mess)
          return
      end if

      if ( .not. matz ) then

        call tred1 ( n, a, w, fv1, fv2 )

        call tqlrat ( n, w, fv2, ierr )

        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RSG - Warning!: Error return from TQLRAT"
          return
        end if

      else

        call tred2 ( n, a, w, fv1, z )

        call tql2 ( n, w, fv1, z, ierr )

        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RSG - Fatal error!: Error return from TQL2"
          return
        end if

        call rebak ( n, b, fv2, n, z )

      end if

    End Subroutine rsg



    !    Subroutine rsgab ( n, a, b, w, matz, z, ierr )
    !
    !! RSGAB computes eigenvalues/vectors, A*B*x=lambda*x, A symmetric, B pos-def.
    !
    !  Discussion:
    !
    !    RSGAB calls the recommended sequence of EISPACK routines
    !    to find the eigenvalues and eigenvectors (if desired)
    !    for the real symmetric generalized eigenproblem
    !      A * B * x = (lambda) *  x.
    !
    !  Modified:
    !
    !    09 February 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrices A and B.
    !
    !    Input, real(kind=dp) A(N,N), contains a real symmetric matrix.
    !
    !    Input, real(kind=dp) B(N,N), contains a positive definite real
    !    symmetric matrix.
    !
    !    Input, logical MATZ, is false if only eigenvalues are desired,
    !    and true if both eigenvalues and eigenvectors are desired.
    !
    !    Output, real(kind=dp) W(N), the eigenvalues in ascending order.
    !
    !    Output, real(kind=dp) Z(N,N), contains the eigenvectors, if MATZ
    !    is true.
    !
    !    Output, integer ::IERR, is set to an error
    !    completion code described in the documentation for TQLRAT and TQL2.
    !    The normal completion code is zero.
    !
    Subroutine rsgab ( n, a, b, w, matz, z, ierr )
      integer,                         intent(in)    :: n
      real(kind = dp), dimension(n,n), intent(in out):: a,b
      real(kind = dp), dimension(n),   intent(out)   :: w
      logical,                         intent(in)    :: matz
      real(kind = dp), dimension(n,n), intent(out)   :: z
      integer,                         intent(out)   :: ierr

      real(kind=dp), dimension(n) :: fv1,fv2

      call reduc2 ( n, a, b, fv2, ierr )

      if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RSGAB - Fatal error!: Error return from REDUC2 -> "//trim(ERR_EisPack_Mess)
          return
      end if

      if ( .not. matz ) then

        call tred1 ( n, a, w, fv1, fv2 )

        call tqlrat ( n, w, fv2, ierr )

        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RSGAB - Fatal error!: Error return from TQLRAT"
          return
        end if

      else

        call tred2 ( n, a, w, fv1, z )

        call tql2 ( n, w, fv1, z, ierr )

        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RSGAB - Fatal error!: Error return from TQL2"
          return
        end if

        call rebak ( n, b, fv2, n, z )

      end if

    End Subroutine rsgab



    !    Subroutine rsgba ( n, a, b, w, matz, z, ierr )
    !
    !! RSGBA computes eigenvalues/vectors, B*A*x=lambda*x, A symmetric, B pos-def.
    !
    !  Discussion:
    !
    !    RSGBA calls the recommended sequence of EISPACK routines
    !    to find the eigenvalues and eigenvectors (if desired)
    !    for the real symmetric generalized eigenproblem:
    !
    !      B * A * x = lambda * x
    !
    !  Modified:
    !
    !    09 February 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrices A and B.
    !
    !    Input, real(kind=dp) A(N,N), a real symmetric matrix.
    !
    !    Input, real(kind=dp) B(N,N), a positive definite symmetric matrix.
    !
    !    Input, logical MATZ, is false if only eigenvalues are desired,
    !    and true if both eigenvalues and eigenvectors are desired.
    !
    !    Output, real(kind=dp) W(N), the eigenvalues in ascending order.
    !
    !    Output, real(kind=dp) Z(N,N), contains the eigenvectors, if MATZ
    !    is true.
    !
    !    Output, integer ::IERR, is set to an error
    !    completion code described in the documentation for TQLRAT and TQL2.
    !    The normal completion code is zero.
    !
    Subroutine rsgba ( n, a, b, w, matz, z, ierr )
      integer,                         intent(in)    :: n
      real(kind = dp), dimension(n,n), intent(in out):: a,b
      real(kind = dp), dimension(n),   intent(out)   :: w
      logical,                         intent(in)    :: matz
      real(kind = dp), dimension(n,n), intent(out)   :: z
      integer,                         intent(out)   :: ierr

      real(kind=dp), dimension(n) :: fv1,fv2

      call reduc2 ( n, a, b, fv2, ierr )

      if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RSGBA - Fatal error!: Error return from REDUC2 -> "//trim(ERR_EisPack_Mess)
          return
      end if

      if ( .not. matz ) then

        call tred1 ( n, a, w, fv1, fv2 )

        call tqlrat ( n, w, fv2, ierr )

        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RSGBA - Fatal error!: Error return from TQLRAT"
          return
        end if

      else

        call tred2 ( n, a, w, fv1, z )

        call tql2 ( n, w, fv1, z, ierr )

        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RSGBA - Fatal error!: Error return from TQL2"
          return
        end if

        call rebakb ( n, b, fv2, n, z )

      end if

    End Subroutine rsgba



    !    Subroutine rsm ( n, a, w, m, z, ierr )
    !
    !! RSM computes eigenvalues, some eigenvectors, real symmetric matrix.
    !
    !  Discussion:
    !
    !    RSM calls the recommended sequence of EISPACK routines
    !    to find all of the eigenvalues and some of the eigenvectors
    !    of a real symmetric matrix.
    !
    !  Modified:
    !
    !    09 February 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input, real(kind=dp) A(N,N), the symmetric matrix.
    !
    !    Input, integer ::M, specifies the number of eigenvectors to
    !    compute.
    !
    !    Output, real(kind=dp) W(N), the eigenvalues in ascending order.
    !
    !    Output, real(kind=dp) Z(N,M), contains the orthonormal eigenvectors
    !    associated with the first M eigenvalues.
    !
    !    Output, integer ::IERR, is set to an error
    !    completion code described in the documentation for TQLRAT, IMTQLV and
    !    TINVIT.  The normal completion code is zero.
    !
    Subroutine rsm ( n, a, w, m, z, ierr )
      integer,                         intent(in)    :: n,m
      real(kind = dp), dimension(n,n), intent(in out):: a
      real(kind = dp), dimension(n),   intent(out)   :: w
      real(kind = dp), dimension(n,m), intent(out)   :: z
      integer,                         intent(out)   :: ierr

      real(kind=dp), dimension(n) :: fwork1,fwork2,fwork3
      integer,       dimension(n) :: iwork

      if ( m <= 0 ) then

        call tred1 ( n, a, w, fwork1, fwork2 )

        call tqlrat ( n, w, fwork2, ierr )

        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RS - Fatal error!: Error return from TQLRAT"
          return
        end if

      else

        call tred1 ( n, a, fwork1, fwork2, fwork3 )

        call imtqlv ( n, fwork1, fwork2, fwork3, w, iwork, ierr )

        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RSM - Fatal error!: Error return from IMTQLV"
          return
        end if

        call tinvit ( n, fwork1, fwork2, fwork3, m, w, iwork, z, ierr )

        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RSM - Fatal error!: Error return from TINVIT"
          return
        end if

        call trbak1 ( n, a, fwork2, m, z )

      end if

    End Subroutine rsm



    !    Subroutine rsp ( n, nv, a, w, matz, z, ierr )
    !
    !! RSP computes eigenvalues and eigenvectors of real symmetric packed matrix.
    !
    !  Discussion:
    !
    !    RSP calls the recommended sequence of EISPACK routines
    !    to find the eigenvalues and eigenvectors (if desired)
    !    of a real symmetric packed matrix.
    !
    !  Modified:
    !
    !    04 March 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input, integer ::NV, the dimension of the array A, which
    !    must be at least (N*(N+1))/2.
    !
    !    Input, real(kind=dp) A(NV), contains the lower triangle of the
    !    real symmetric packed matrix stored row-wise.
    !
    !    Input, logical MATZ, is false if only eigenvalues are desired,
    !    and true if both eigenvalues and eigenvectors are desired.
    !
    !    Output, real(kind=dp) W(N), the eigenvalues in ascending order.
    !
    !    Output, real(kind=dp) Z(N,N), contains the eigenvectors, if MATZ is
    !    true.
    !
    !    Output, integer ::IERR, is set to an error
    !    completion code described in the documentation for TQLRAT and TQL2.
    !    The normal completion code is zero.
    !
    Subroutine rsp ( n, nv, a, w, matz, z, ierr )
      integer,                         intent(in)    :: n,nv
      real(kind = dp), dimension(nv),  intent(in out):: a
      real(kind = dp), dimension(n),   intent(out)   :: w
      logical,                         intent(in)    :: matz
      real(kind = dp), dimension(n,n), intent(out)   :: z
      integer,                         intent(out)   :: ierr


      real(kind=dp), dimension(n) :: fv1, fv2

      if ( nv < ( n * ( n + 1 ) ) / 2 ) then
        ierr = 20 * n
        ERR_EisPack=.true.
        ERR_EisPack_Mess=" RSP - Fatal error!: NV too small"
        return
      end if

      call tred3 ( n, nv, a, w, fv1, fv2 )

      if ( .not. matz ) then

        call tqlrat ( n, w, fv2, ierr )

        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RSP - Fatal error!: Error return from TQLRAT"
          return
        end if

      else

        call r8mat_identity ( n, z )

        call tql2 ( n, w, fv1, z, ierr )

        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RSP - Fatal error!: Error return from TQL2"
          return
        end if

        call trbak3 ( n, nv, a, n, z )

      end if

    End Subroutine rsp



    !    Subroutine rspp ( n, nv, a, w, matz, z, m, typ, ierr )
    !
    !! RSPP computes some eigenvalues/vectors, real symmetric packed matrix.
    !
    !  Discussion:
    !
    !    RSPP calls the appropriate routines for the following problem:
    !
    !    Given a symmetric matrix A, which is stored in a packed mode, find
    !    the M smallest or largest eigenvalues, and corresponding eigenvectors.
    !
    !    The routine RSP returns all eigenvalues and eigenvectors.
    !
    !  Modified:
    !
    !    10 February 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of A, the number of rows and
    !    columns in the original matrix.
    !
    !    Input, integer ::NV, is the of the array A as specified in the
    !    calling program.  NV must not be less than N*(N+1)/2.
    !
    !    Input, real(kind=dp) A((N*(N+1))/2), on input the lower triangle of the
    !    real symmetric matrix, stored row-wise in the vector,
    !    in the order A(1,1), / A(2,1), A(2,2), / A(3,1), A(3,2), A(3,3)/
    !    and so on.
    !
    !    Output, real(kind=dp) W(M), the eigenvalues requested.
    !
    !    Input, logical MATZ, is false if only eigenvalues are desired,
    !    and true if both eigenvalues and eigenvectors are desired.
    !
    !    Output, real(kind=dp) Z(N,M), the eigenvectors.
    !
    !    Input, integer ::M, the number of eigenvalues to be found.
    !
    !    Input, logical typ, set to .true. if the smallest eigenvalues
    !    are to be found, or .false. if the largest ones are sought.
    !
    !    Output, integer ::IERR, error flag from RATQR.  IERR=0 on
    !    normal return.  IERR nonzero, in this case, means that the algorithm broke
    !    down while computing an eigenvalue.
    !
    Subroutine rspp ( n, nv, a, w, matz, z, m, typ, ierr )
      integer,                         intent(in)    :: n,nv,m
      real(kind = dp), dimension(nv),  intent(in out):: a
      real(kind = dp), dimension(m),   intent(out)   :: w
      logical,                         intent(in)    :: matz
      real(kind = dp), dimension(n,m), intent(out)   :: z
      logical,                         intent(in)    :: typ
      integer,                         intent(out)   :: ierr

      real(kind=dp) :: eps1
      integer       :: idef
      integer,      dimension(n) :: iwork
      real(kind=dp),dimension(n) :: work1,work2,work3,bd
    !
    !  IDEF =
    !    -1 if the matrix is known to be negative definite,
    !    +1 if the matrix is known to be positive definite, or
    !    0 otherwise.
    !
      idef = 0
    !
    !  Reduce to symmetric tridiagonal form.
    !
      call tred3 ( n, nv, a, work1, work2, work3 )
    !
    !  Find the eigenvalues.
    !
      eps1 = 0.0_dp

      call ratqr ( n, eps1, work1, work2, work3, m, w, iwork, &
        bd, typ, idef, ierr )

      if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RSPP - Fatal error!: Error return from RATQR"
          return
      end if
    !
    !  Find eigenvectors for the first M eigenvalues.
    !
      if ( matz ) then

        call tinvit ( n, work1, work2, work3, m, w, iwork, z, ierr )

        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RSPP - Fatal error!: Error return from TINVIT"
          return
        end if
    !
    !  Reverse the transformation.
    !
        call trbak3 ( n, nv, a, m, z )

      end if

    End Subroutine rspp



    !    Subroutine rst ( n, w, e, matz, z, ierr )
    !
    !! RST computes eigenvalues/vectors, real symmetric tridiagonal matrix.
    !
    !  Discussion:
    !
    !    RST calls the recommended sequence of EISPACK routines
    !    to find the eigenvalues and eigenvectors (if desired)
    !    of a real symmetric tridiagonal matrix.
    !
    !  Modified:
    !
    !    04 March 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input/output, real(kind=dp) W(N).  On input, the diagonal elements
    !    of the real symmetric tridiagonal matrix.  On output, the eigenvalues in
    !    ascending order.
    !
    !    Input/output (may be modified within imtql1,2), real(kind=dp) E(N), the subdiagonal
    !    elements of the matrix in E(2:N).  E(1) is arbitrary.
    !
    !    Input, logical MATZ, is false if only eigenvalues are desired,
    !    and true if both eigenvalues and eigenvectors are desired.
    !
    !    Output, real(kind=dp) Z(N,N), contains the eigenvectors, if MATZ
    !    is true.
    !
    !    Output, integer ::IERR, is set to an error
    !    completion code described in the documentation for IMTQL1 and IMTQL2.
    !    The normal completion code is zero.
    !
    Subroutine rst ( n, w, e, matz, z, ierr )
      integer,                         intent(in)    :: n
      real(kind = dp), dimension(n),   intent(in out):: w
      real(kind = dp), dimension(n),   intent(in out):: e
      logical,                         intent(in)    :: matz
      real(kind = dp), dimension(n,n), intent(out)   :: z
      integer,                         intent(out)   :: ierr

      if ( .not. matz ) then

        call imtql1 ( n, w, e, ierr )

        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RST - Fatal error!: Error return from IMTQL1"
          return
        end if

      else

        call r8mat_identity ( n, z )

        call imtql2 ( n, w, e, z, ierr )

        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RST - Fatal error!: Error return from IMTQL2"
          return
        end if

      end if

    End Subroutine rst



    !    Subroutine rt ( n, a, w, matz, z, ierr )
    !
    !! RT computes eigenvalues/vectors, real sign-symmetric tridiagonal matrix.
    !
    !  Discussion:
    !
    !    RT calls the recommended sequence of EISPACK routines
    !    to find the eigenvalues and eigenvectors (if desired)
    !    of a special real tridiagonal matrix.
    !
    !  Modified:
    !
    !    04 March 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input, real(kind=dp) A(N,3), contains the special real tridiagonal
    !    matrix in its first three columns.  The subdiagonal elements are stored
    !    in the last N-1 positions of the first column, the diagonal elements
    !    in the second column, and the superdiagonal elements in the first N-1
    !    positions of the third column.  Elements A(1,1) and A(N,3) are arbitrary.
    !
    !    Input, logical MATZ, is false if only eigenvalues are desired,
    !    and true if both eigenvalues and eigenvectors are desired.
    !
    !    Output, real(kind=dp) W(N), the eigenvalues in ascending order.
    !
    !    Output, real(kind=dp) Z(N,N), contains the eigenvectors, if MATZ
    !    is true.
    !
    !    Output, integer ::IERR, is set to an error
    !    completion code described in the documentation for IMTQL1 and IMTQL2.
    !    The normal completion code is zero.
    !
    Subroutine rt ( n, a, w, matz, z, ierr )
      integer,                         intent(in)    :: n
      real(kind = dp), dimension(n,3), intent(in)    :: a
      real(kind = dp), dimension(n),   intent(out)   :: w
      logical,                         intent(in)    :: matz
      real(kind = dp), dimension(n,n), intent(out)   :: z
      integer,                         intent(out)   :: ierr

      real(kind=dp), dimension(n) :: fv1

      if ( .not. matz ) then

        call figi ( n, a, w, fv1, fv1, ierr )

        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RT - Fatal error!: Error return from FIGI"
          return
        end if

        call imtql1 ( n, w, fv1, ierr )

        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RT - Fatal error!: Error return from IMTQL1"
          return
        end if

      else

        call figi2 ( n, a, w, fv1, z, ierr )

        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RT - Fatal error!: Error return from FIGI2"
          return
        end if

        call imtql2 ( n, w, fv1, z, ierr )

        if ( ierr /= 0 ) then
          ERR_EisPack=.true.
          ERR_EisPack_Mess=" RT - Fatal error!: Error return from IMTQL2"
          return
        end if

      end if

    End Subroutine rt



    !    function sturm_sequence ( d, e, e2, n, p, q, x1 ) result(sturm_sequence_val)
    !
    !! STURM_SEQUENCE counts eigenvalues of a symmetric tridiagonal submatrix.
    !
    !  Discussion:
    !
    !    Let A be a symmetric tridiagonal matrix, and consider the submatrix
    !    defined by indices P through Q.
    !
    !    STURM_SEQUENCE will determine the number of eigenvalues associated
    !    with this submatrix which are no larger than a given upper bound X1.
    !
    !  Modified:
    !
    !    19 February 2018
    !
    !  Arguments:
    !
    !    Input, real(kind=dp) D(N), the diagonal elements of the input matrix.
    !
    !    Input, real(kind=dp) E(N), the subdiagonal elements of the matrix.
    !    E(1) is arbitrary.
    !
    !    Input, real(kind=dp) E2(N), the squares of the corresponding
    !    elements of E, with zeros corresponding to negligible elements of E.
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input, integer ::P, Q, the lower and upper limits on the
    !    submatrix.
    !
    !    Input, real(kind=dp) X1, an upper bound.
    !
    !    Output, integer ::STURM_SEQUENCE, the number of eigenvalues
    !    less than or equal to X1.
    !
    function sturm_sequence ( d, e, e2, n, p, q, x1 ) result(sturm_sequence_val)
      real(kind = dp), dimension(n),   intent(in)    :: d,e,e2
      integer,                         intent(in)    :: n,p,q
      real(kind = dp),                 intent(in)    :: x1
      integer                                        :: sturm_sequence_val

      integer ::i,s

      real(kind=dp) :: u,v

      s = p - 1
      u = 1.0_dp

      do i = p, q

        if ( u == 0.0_dp ) then
          v = abs ( e(i) ) / epsilon ( v )
          if ( e2(i) == 0.0_dp ) then
            v = 0.0_dp
          end if
        else
          v = e2(i) / u
        end if

        u = d(i) - x1 - v
        if ( u < 0.0_dp ) then
          s = s + 1
        end if

      end do

      sturm_sequence_val = s

      return
    End function sturm_sequence



    ! Subroutine svd ( m, n, a, w, matu, u, matv, v, ierr )
    !
    !! SVD computes the singular value decomposition for a real matrix.
    !
    !  Discussion:
    !
    !    SVD determines the singular value decomposition
    !
    !      A = U * S * V'
    !
    !    of a real M by N rectangular matrix.  Householder bidiagonalization
    !    and a variant of the QR algorithm are used.
    !
    !  Modified:
    !
    !    01 February 2018
    !
    !  Reference:
    !
    !    Golub, Christian Reinsch,
    !    Numerische Mathematik,
    !    Volume 14, 1970, pages 403-420.
    !
    !  Arguments:
    !
    !    Input, integer ::M, the number of rows of A and U.
    !
    !    Input, integer ::N, the number of columns of A and U, and
    !    the order of V.
    !
    !    Input, real(kind=dp) A(M,N), the M by N matrix to be decomposed.
    !
    !    Output, real(kind=dp) W(N), the singular values of A.  These are the
    !    diagonal elements of S.  They are unordered.  If an error exit is
    !    made, the singular values should be correct for indices
    !    IERR+1, IERR+2,..., N.
    !
    !    Input, logical MATU, should be set to TRUE if the U matrix in the
    !    decomposition is desired, and to FALSE otherwise.
    !
    !    Output, real(kind=dp) U(M,N), contains the matrix U, with orthogonal
    !    columns, of the decomposition, if MATU has been set to TRUE.  Otherwise
    !    U is used as a temporary array.  U may coincide with A.
    !    If an error exit is made, the columns of U corresponding
    !    to indices of correct singular values should be correct.
    !
    !    Input, logical MATV, should be set to TRUE if the V matrix in the
    !    decomposition is desired, and to FALSE otherwise.
    !
    !    Output, real(kind=dp) V(N,N), the orthogonal matrix V of the
    !    decomposition if MATV has been set to TRUE.  Otherwise V is not referenced.
    !    V may also coincide with A if U is not needed.  If an error
    !    exit is made, the columns of V corresponding to indices of
    !    correct singular values should be correct.
    !
    !    Output, integer ::IERR, error flag.
    !    0, for normal return,
    !    K, if the K-th singular value has not been determined after 30 iterations.
    !
    Subroutine svd ( m, n, a, w, matu, u, matv, v, ierr )
      integer,                      intent(in)  :: m, n
      real(kind=dp), dimension(m,n),intent(in)  :: a
      real(kind=dp), dimension(n),  intent(out) :: w
      logical,                      intent(in)  :: matu
      real(kind=dp), dimension(m,n),intent(out) :: u
      logical,                      intent(in)  :: matv
      real(kind=dp), dimension(n,n),intent(out) :: v
      integer,                      intent(out) :: ierr
      !
      integer :: i,its,i1,j,k,l,ll,l1,mn
      logical :: skip
      real(kind=dp), dimension(n) :: rv1
      real(kind=dp) :: s,scal,tst1,tst2,x,y,z,c,f,g,h

      ierr = 0
      u(1:m,1:n) = a(1:m,1:n)
    !
    !  Householder reduction to bidiagonal form.
    !
      g = 0.0_dp
      scal = 0.0_dp
      x = 0.0_dp

      do i = 1, n

        l = i + 1
        rv1(i) = scal * g
        g = 0.0_dp
        s = 0.0_dp
        scal = 0.0_dp

        if ( i <= m ) then

          scal = sum ( abs ( u(i:m,i) ) )

          if ( scal /= 0.0_dp ) then

            u(i:m,i) = u(i:m,i) / scal
            s = sum ( u(i:m,i)**2 )
            f = u(i,i)
            g = - sign ( sqrt ( s ), f )
            h = f * g - s
            u(i,i) = f - g

            do j = i + 1, n
              s = dot_product ( u(i:m,i), u(i:m,j) )
              u(i:m,j) = u(i:m,j) + s * u(i:m,i) / h
            end do

            u(i:m,i) = scal * u(i:m,i)

          end if

        end if

        w(i) = scal * g
        g = 0.0_dp
        s = 0.0_dp
        scal = 0.0_dp

        if ( i <= m .and. i /= n ) then

          scal = sum ( abs ( u(i,i+1:n) ) )

          if ( scal /= 0.0_dp ) then

            u(i,i+1:n) = u(i,i+1:n) / scal
            s = sum ( u(i,i+1:n)**2 )
            f = u(i,i+1)
            g = - sign ( sqrt ( s ), f )
            h = f * g - s
            u(i,i+1) = f - g
            rv1(i+1:n) = u(i,i+1:n) / h

            do j = i + 1, m
              s = dot_product ( u(j,i+1:n), u(i,i+1:n) )
              u(j,i+1:n) = u(j,i+1:n) + s * rv1(i+1:n)
            end do
            u(i,i+1:n) = u(i,i+1:n) * scal

          end if

        end if

        x = max ( x, abs ( w(i) ) + abs ( rv1(i) ) )

      end do
    !
    !  Accumulation of right-hand transformations.
    !
      if ( matv ) then

        do i = n, 1, -1

          if ( i < n ) then

            if ( g /= 0.0_dp ) then

              v(i+1:n,i) = ( u(i,i+1:n) / u(i,i+1) ) / g

              do j = i + 1, n
                s = dot_product ( u(i,i+1:n), v(i+1:n,j) )
                v(i+1:n,j) = v(i+1:n,j) + s * v(i+1:n,i)
              end do

            end if

            v(i,i+1:n) = 0.0_dp
            v(i+1:n,i) = 0.0_dp

          end if

          v(i,i) = 1.0_dp
          g = rv1(i)

        end do

      end if
    !
    !  Accumulation of left-hand transformations.
    !
      if ( matu ) then

        mn = min ( m, n )

        do i = min ( m, n ), 1, -1

          l = i + 1
          g = w(i)

          if ( i /= n ) then
            u(i,i+1:n) = 0.0_dp
          end if

          if ( g /= 0.0_dp ) then

            if ( i /= mn ) then

              do j = i + 1, n
                s = dot_product ( u(i+1:m,i), u(i+1:m,j) )
                f = ( s / u(i,i) ) / g
                u(i:m,j) = u(i:m,j) + f * u(i:m,i)
              end do

            end if

            u(i:m,i) = u(i:m,i) / g

          else

            u(i:m,i) = 0.0_dp

          end if

          u(i,i) = u(i,i) + 1.0_dp

        end do

      end if
    !
    !  Diagonalization of the bidiagonal form.
    !
      tst1 = x

      do k = n, 1, -1

        its = 0
    !
    !  Test for splitting.
    !
        do

          skip = .false.

          do ll = 1, k

            l1 = k - ll
            l = l1 + 1
            tst2 = tst1 + abs ( rv1(l) )

            if ( tst2 == tst1 ) then
              skip = .true.
              exit
            end if

            tst2 = tst1 + abs ( w(l-1) )

            if ( tst2 == tst1 ) then
              exit
            end if

          end do
    !
    !  Cancellation of rv1(l) if L greater than 1.
    !
          if ( .not. skip ) then

            c = 0.0_dp
            s = 1.0_dp

            do i = l, k

              f = s * rv1(i)
              rv1(i) = c * rv1(i)
              tst2 = tst1 + abs ( f )

              if ( tst2 == tst1 ) then
                exit
              end if

              g = w(i)
              h = pythag ( f, g )
              w(i) = h
              c = g / h
              s = - f / h

              if ( matu ) then

                do j = 1, m
                  y = u(j,l1)
                  z = u(j,i)
                  u(j,l1) =  y * c + z * s
                  u(j,i) = - y * s + z * c
                end do

              end if

            end do

          end if

          z = w(k)
    !
    !  Convergence.
    !
          if ( l == k ) then

            if ( z <= 0.0_dp ) then
              w(k) = - z
              if ( matv ) then
                v(1:n,k) = - v(1:n,k)
              end if
            end if

            exit

          end if
    !
    !  Shift from bottom 2 by 2 minor.
    !
          if ( 30 <= its ) then
            ierr = k
            return
          end if

          its = its + 1
          x = w(l)
          y = w(k-1)
          g = rv1(k-1)
          h = rv1(k)
          f = 0.5_dp * ( ( ( g + z ) / h ) * ( ( g - z ) / y ) + y / h - h / y )
          g = pythag ( f, 1.0_dp )
          f = x - ( z / x ) * z + ( h / x ) * ( y / ( f + sign ( g, f ) ) - h )
    !
    !  Next QR transformation.
    !
          c = 1.0_dp
          s = 1.0_dp

          do i1 = l, k - 1

            i = i1 + 1
            g = rv1(i)
            y = w(i)
            h = s * g
            g = c * g
            z = pythag ( f, h )
            rv1(i-1) = z
            c = f / z
            s = h / z
            f =   x * c + g * s
            g = - x * s + g * c
            h = y * s
            y = y * c

            if ( matv ) then

              do j = 1, n
                x = v(j,i1)
                z = v(j,i)
                v(j,i1) =  x * c + z * s
                v(j,i) = - x * s + z * c
              end do

            end if

            z = pythag ( f, h )
            w(i1) = z
    !
    !  Rotation can be arbitrary if Z is zero.
    !
            if ( z /= 0.0_dp ) then
              c = f / z
              s = h / z
            end if

            f =   c * g + s * y
            x = - s * g + c * y

            if ( matu ) then

              do j = 1, m
                y = u(j,i1)
                z = u(j,i)
                u(j,i1) =  y * c + z * s
                u(j,i) = - y * s + z * c
              end do

            end if

          end do

          rv1(l) = 0.0_dp
          rv1(k) = f
          w(k) = x

        end do

      end do

    End Subroutine svd



    !    Subroutine timestamp ( )
    !
    !! TIMESTAMP prints the current YMDHMS date as a time stamp.
    !
    !  Example:
    !
    !    31 May 2001   9:45:54.872 AM
    !
    !
    !  Modified:
    !
    !    18 May 2013
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Arguments:
    !
    !    None
    !
    Subroutine timestamp ( )
      character ( len = 8 ) ampm
      integer ::d
      integer ::h
      integer ::m
      integer ::mm
      character ( len = 9 ), parameter, dimension(12) :: month = (/ &
        'January  ', 'February ', 'March    ', 'April    ', &
        'May      ', 'June     ', 'July     ', 'August   ', &
        'September', 'October  ', 'November ', 'December ' /)
      integer ::n
      integer ::s
      integer ::values(8)
      integer ::y

      call date_and_time ( values = values )

      y = values(1)
      m = values(2)
      d = values(3)
      h = values(5)
      n = values(6)
      s = values(7)
      mm = values(8)

      if ( h < 12 ) then
        ampm = 'AM'
      else if ( h == 12 ) then
        if ( n == 0 .and. s == 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h < 12 ) then
          ampm = 'PM'
        else if ( h == 12 ) then
          if ( n == 0 .and. s == 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write(*, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
        d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

    End Subroutine timestamp


    !    Subroutine tinvit ( n, d, e, e2, m, w, ind, z, ierr )
    !
    !! TINVIT computes eigenvectors from eigenvalues, real tridiagonal symmetric.
    !
    !  Discussion:
    !
    !    TINVIT finds eigenvectors of a tridiagonal symmetric matrix corresponding
    !    to specified eigenvalues using inverse iteration.
    !
    !  Modified:
    !
    !    03 February 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input, real(kind=dp) D(N), the diagonal elements of the matrix.
    !
    !    Input, real(kind=dp) E(N), contains the subdiagonal elements of
    !    the input matrix in E(2:N).  E(1) is arbitrary.
    !
    !    Input, real(kind=dp) E2(N), contains the squares of the corresponding
    !    elements of E, with zeros corresponding to negligible elements of E.
    !    E(I) is considered negligible if it is not larger than the product of
    !    the relative machine precision and the sum of the magnitudes of D(I)
    !    and D(I-1).  E2(1) must contain 0.0_dp if the eigenvalues are in
    !    ascending order, or 2.0_dp if the eigenvalues are in descending order.
    !    If BISECT, TRIDIB, or IMTQLV has been used to find the eigenvalues,
    !    their output E2 array is exactly what is expected here.
    !
    !    Input, integer ::M, the number of specified eigenvalues.
    !
    !    Input, real(kind=dp) W(M), the eigenvalues.
    !
    !    Input, integer ::IND(M), the submatrix indices associated with
    !    the corresponding eigenvalues in W: 1 for eigenvalues belonging to the
    !    first submatrix from the top, 2 for those belonging to the second
    !    submatrix, and so on.
    !
    !    Output, real(kind=dp) Z(N,M), the associated set of orthonormal
    !    eigenvectors.  Any vector which fails to converge is set to zero.
    !
    !    Output, integer ::IERR, error flag.
    !    0, for normal return,
    !    -R, if the eigenvector corresponding to the R-th eigenvalue fails to
    !      converge in 5 iterations.
    !
    Subroutine tinvit ( n, d, e, e2, m, w, ind, z, ierr )
      integer,                       intent(in)     :: n
      real(kind=dp), dimension(n),   intent(in)     :: d,e
      real(kind=dp), dimension(n),   intent(in)     :: e2
      integer,                       intent(in)     :: m
      real(kind=dp), dimension(m),   intent(in)     :: w
      integer,       dimension(m),   intent(in)     :: ind
      real(kind=dp), dimension(n,m), intent(out)    :: z
      integer,                       intent(out)    :: ierr

      integer :: group,i,its,j,jj,p,q,r,s,tag
      real(kind=dp) :: eps2,eps3,eps4,norm,order,u,uk,v,x0,x1,xu
      real(kind=dp), dimension(n) :: rv1 ,rv2,rv3,rv4,rv6

      ierr = 0

      if ( m == 0 ) return

      u = 0.0_dp
      x0 = 0.0_dp

      tag = 0
      order = 1.0_dp - e2(1)
      q = 0
    !
    !  Establish and process next submatrix.
    !
      do

        p = q + 1

        do q = p, n
          if ( q == n ) then
            exit
          end if
          if ( e2(q+1) == 0.0_dp ) then
            exit
          end if
        end do
    !
    !  Find vectors by inverse iteration.
    !
        tag = tag + 1
        s = 0

        do r = 1, m

          if ( ind(r) /= tag ) then
            cycle
          end if

          its = 1
          x1 = w(r)
    !
    !  Look for close or coincident roots.
    !
          if ( s /= 0 ) then

            if ( eps2 <= abs ( x1 - x0 ) ) then
              group = 0
            else
              group = group + 1
              if ( order * ( x1 - x0 ) <= 0.0_dp ) then
                x1 = x0 + order * eps3
              end if
            end if

          else
    !
    !  Check for isolated root.
    !
            xu = 1.0_dp

            if ( p == q ) then
              rv6(p) = 1.0_dp
              z(1:n,r) = 0.0_dp
              z(p:q,r) = rv6(p:q) * xu
              x0 = x1
              cycle
            end if

            norm = abs ( d(p) )

            do i = p + 1, q
              norm = max ( norm, abs ( d(i) ) + abs ( e(i) ) )
            end do
    !
    !  EPS2 is the criterion for grouping,
    !  EPS3 replaces zero pivots and equal roots are modified by EPS3,
    !  EPS4 is taken very small to avoid overflow.
    !
            eps2 = 0.001_dp * norm
            eps3 = abs ( norm ) * epsilon ( eps3 )
            uk = q - p + 1
            eps4 = uk * eps3
            uk = eps4 / sqrt ( uk )
            s = p
            group = 0

          end if
    !
    !  Elimination with interchanges and initialization of vector.
    !
          v = 0.0_dp

          do i = p, q

            rv6(i) = uk

            if ( i == p ) then

              u = d(i) - x1 - xu * v
              if ( i /= q ) then
                v = e(i+1)
              end if

            else if ( abs ( u ) <= abs ( e(i) ) ) then

              xu = u / e(i)
              rv4(i) = xu
              rv1(i-1) = e(i)
              rv2(i-1) = d(i) - x1
              rv3(i-1) = 0.0_dp
              if ( i /= q ) then
                rv3(i-1) = e(i+1)
              end if
              u = v - xu * rv2(i-1)
              v = - xu * rv3(i-1)

            else

              xu = e(i) / u
              rv4(i) = xu
              rv1(i-1) = u
              rv2(i-1) = v
              rv3(i-1) = 0.0_dp

              u = d(i) - x1 - xu * v
              if ( i /= q ) then
                v = e(i+1)
              end if

            end if

          end do

          if ( u == 0.0_dp ) then
            u = eps3
          end if

          rv1(q) = u
          rv2(q) = 0.0_dp
          rv3(q) = 0.0_dp
    !
    !  Back substitution.
    !
          do

            do i = q, p, -1
              rv6(i) = ( rv6(i) - u * rv2(i) - v * rv3(i) ) / rv1(i)
              v = u
              u = rv6(i)
            end do
    !
    !  Orthogonalize with respect to previous members of group.
    !
            j = r

            do jj = 1, group

              do

                j = j - 1

                if ( ind(j) == tag ) then
                  exit
                end if

              end do

              xu = dot_product ( rv6(p:q), z(p:q,j) )

              rv6(p:q) = rv6(p:q) - xu * z(p:q,j)

            end do

            norm = sum ( abs ( rv6(p:q) ) )
    !
    !  Normalize so that sum of squares is 1.
    !
            if ( 1.0_dp <= norm ) then

              u = 0.0_dp
              do i = p, q
                u = pythag ( u, rv6(i) )
              end do

              xu = 1.0_dp / u

              z(1:n,r) = 0.0_dp
              z(p:q,r) = rv6(p:q) * xu
              x0 = x1
              exit
    !
    !  Set error: non-converged eigenvector.
    !
            else if ( 5 <= its ) then

              ierr = - r
              xu = 0.0_dp
              z(1:n,r) = 0.0_dp
              z(p:q,r) = rv6(p:q) * xu
              x0 = x1
              exit

            else

              if ( norm == 0.0_dp ) then
                rv6(s) = eps4
                s = s + 1
                if ( q < s ) then
                  s = p
                end if
              else
                xu = eps4 / norm
                rv6(p:q) = rv6(p:q) * xu
              end if
    !
    !  If RV1(I-1) == E(I), a row interchange was performed earlier in the
    !  triangularization process.
    !
              do i = p + 1, q

                u = rv6(i)

                if ( rv1(i-1) == e(i) ) then
                  u = rv6(i-1)
                  rv6(i-1) = rv6(i)
                end if

                rv6(i) = u - rv4(i) * rv6(i-1)

              end do

              its = its + 1

            end if

          end do

        end do

        if ( n <= q ) then
          exit
        end if

      end do

    End Subroutine tinvit



    !    Subroutine tql1 ( n, d, e, ierr )
    !
    !! TQL1 computes all eigenvalues of a real symmetric tridiagonal matrix.
    !
    !  Discussion:
    !
    !    TQL1 finds the eigenvalues of a symmetric tridiagonal
    !    matrix by the QL method.
    !
    !  Modified:
    !
    !    31 January 2018
    !
    !
    !  References:
    !
    !    Bowdler, Martin, Reinsch, James Wilkinson,
    !    Numerische Mathematik,
    !    Volume 11, 1968, pages 293-306.
    !
    !  Arguments:
    !
    !    Input, integer ::N, is the order of the matrix.
    !
    !    Input/output, real(kind=dp) D(N).
    !    On input, the diagonal elements of the matrix.
    !    On output, the eigenvalues in ascending order.
    !    If an error exit is made, the eigenvalues are correct and
    !    ordered for indices 1, 2,... IERR-1, but may not be
    !    the smallest eigenvalues.
    !
    !    Input/output, real(kind=dp) E(N).  On input, E(2:N) contains the
    !    subdiagonal elements of the input matrix, and E(1) is arbitrary.
    !    On output, E has been destroyed.
    !
    !    Output, integer ::IERR, error flag.
    !    0, normal return,
    !    J, if the J-th eigenvalue has not been determined after 30 iterations.
    !
    Subroutine tql1 ( n, d, e, ierr )
      integer,                       intent(in)     :: n
      real(kind=dp), dimension(n),   intent(in out) :: d,e
      integer,                       intent(out)    :: ierr

      real(kind=dp) :: c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2
      integer :: i,its,l,l1,l2,m

      ierr = 0
      if ( n == 1 ) return

      do i = 2, n
        e(i-1) = e(i)
      end do

      f = 0.0_dp
      tst1 = 0.0_dp
      e(n) = 0.0_dp

      do l = 1, n

        its = 0
        h = abs ( d(l) ) + abs ( e(l) )
        tst1 = max ( tst1, h )
    !
    !  Look for a small sub-diagonal element.
    !
        do m = l, n

          tst2 = tst1 + abs ( e(m) )

          if ( tst2 == tst1 ) then
            exit
          end if

        end do

        if ( m /= l ) then

          do

            if ( 30 <= its ) then
              ierr = l
              return
            end if

            its = its + 1
    !
    !  Form the shift.
    !
            l1 = l + 1
            l2 = l1 + 1
            g = d(l)
            p = ( d(l1) - g ) / ( 2.0_dp * e(l) )
            r = pythag ( p, 1.0_dp )
            d(l) = e(l) / ( p + sign ( r, p ) )
            d(l1) = e(l) * ( p + sign ( r, p ) )
            dl1 = d(l1)
            h = g - d(l)

            d(l2:n) = d(l2:n) - h

            f = f + h
    !
    !  QL transformation.
    !
            p = d(m)
            c = 1.0_dp
            c2 = c
            el1 = e(l1)
            s = 0.0_dp

            do i = m - 1, l, -1
              c3 = c2
              c2 = c
              s2 = s
              g = c * e(i)
              h = c * p
              r = pythag ( p, e(i) )
              e(i+1) = s * r
              s = e(i) / r
              c = p / r
              p = c * d(i) - s * g
              d(i+1) = h + s * ( c * g + s * d(i) )
            end do

            p = - s * s2 * c3 * el1 * e(l) / dl1
            e(l) = s * p
            d(l) = c * p
            tst2 = tst1 + abs ( e(l) )

            if ( tst2 <= tst1 ) then
              exit
            end if

          end do

        end if

        p = d(l) + f
    !
    !  Order the eigenvalues.
    !
        do i = l, 1, -1
          if ( i == 1 ) then
            d(i) = p
          else if ( d(i-1) <= p ) then
            d(i) = p
            exit
          else
            d(i) = d(i-1)
          end if
        end do

      end do

    End Subroutine tql1



    !    Subroutine tql2 ( n, d, e, z, ierr )
    !
    !! TQL2 computes all eigenvalues/vectors, real symmetric tridiagonal matrix.
    !
    !  Discussion:
    !
    !    TQL2 finds the eigenvalues and eigenvectors of a symmetric
    !    tridiagonal matrix by the QL method.  The eigenvectors of a full
    !    symmetric matrix can also be found if TRED2 has been used to reduce this
    !    full matrix to tridiagonal form.
    !
    !  Modified:
    !
    !    04 March 2018
    !
    !
    !  Reference:
    !
    !    Bowdler, Martin, Reinsch, James Wilkinson,
    !    TQL2,
    !    Numerische Mathematik,
    !    Volume 11, pages 293-306, 1968.
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input/output, real(kind=dp) D(N).  On input, the diagonal elements of
    !    the matrix.  On output, the eigenvalues in ascending order.  If an error
    !    exit is made, the eigenvalues are correct but unordered for indices
    !    1,2,...,IERR-1.
    !
    !    Input/output, real(kind=dp) E(N).  On input, E(2:N) contains the
    !    subdiagonal elements of the input matrix, and E(1) is arbitrary.
    !    On output, E has been destroyed.
    !
    !    Input/output, real(kind=dp) Z(N,N).  On input, the transformation matrix
    !    produced in the reduction by TRED2, if performed.  If the eigenvectors of
    !    the tridiagonal matrix are desired, Z must contain the identity matrix.
    !    On output, Z contains the orthonormal eigenvectors of the symmetric
    !    tridiagonal (or full) matrix.  If an error exit is made, Z contains
    !    the eigenvectors associated with the stored eigenvalues.
    !
    !    Output, integer ::IERR, error flag.
    !    0, normal return,
    !    J, if the J-th eigenvalue has not been determined after
    !    30 iterations.
    !
    Subroutine tql2 ( n, d, e, z, ierr )
      integer,                       intent(in)     :: n
      real(kind=dp), dimension(n),   intent(in out) :: d,e
      real(kind=dp), dimension(n,n), intent(in out) :: z
      integer,                       intent(out)    :: ierr

      integer :: i,its,j,k,l,l1,l2,m,mml
      real(kind=dp) :: c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,t,tst1,tst2

      ierr = 0

      if ( n == 1 ) return

      do i = 2, n
        e(i-1) = e(i)
      end do

      f = 0.0_dp
      tst1 = 0.0_dp
      e(n) = 0.0_dp

      do l = 1, n

        its = 0
        h = abs ( d(l) ) + abs ( e(l) )
        tst1 = max ( tst1, h )
    !
    !  Look for a small sub-diagonal element.
    !
        do m = l, n
          tst2 = tst1 + abs ( e(m) )
          if ( tst2 == tst1 ) then
            exit
          end if
        end do

        if ( m /= l ) then

          do

            if ( 30 <= its ) then
              ierr = l
              return
            end if

            its = its + 1
    !
    !  Form shift.
    !
            l1 = l + 1
            l2 = l1 + 1
            g = d(l)
            p = ( d(l1) - g ) / ( 2.0_dp * e(l) )
            r = pythag ( p, 1.0_dp )
            d(l) = e(l) / ( p + sign ( r, p ) )
            d(l1) = e(l) * ( p + sign ( r, p ) )
            dl1 = d(l1)
            h = g - d(l)
            d(l2:n) = d(l2:n) - h
            f = f + h
    !
    !  QL transformation.
    !
            p = d(m)
            c = 1.0_dp
            c2 = c
            el1 = e(l1)
            s = 0.0_dp
            mml = m - l

            do i = m - 1, l, -1

              c3 = c2
              c2 = c
              s2 = s
              g = c * e(i)
              h = c * p
              r = pythag ( p, e(i) )
              e(i+1) = s * r
              s = e(i) / r
              c = p / r
              p = c * d(i) - s * g
              d(i+1) = h + s * ( c * g + s * d(i) )
    !
    !  Form vector.
    !
              do k = 1, n
                h = z(k,i+1)
                z(k,i+1) = s * z(k,i) + c * h
                z(k,i) = c * z(k,i) - s * h
              end do

            end do

            p = - s * s2 * c3 * el1 * e(l) / dl1
            e(l) = s * p
            d(l) = c * p
            tst2 = tst1 + abs ( e(l) )

            if ( tst2 <= tst1 ) then
              exit
            end if

          end do

        end if

        d(l) = d(l) + f

      end do
    !
    !  Order eigenvalues and eigenvectors.
    !
      do i = 1, n - 1

        k = i
        p = d(i)

        do j = i + 1, n

          if ( d(j) < p ) then
            k = j
            p = d(j)
          end if

        end do

        if ( k /= i ) then

          d(k) = d(i)
          d(i) = p

          do j = 1, n
            t      = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = t
          end do

        end if

      end do

    End Subroutine tql2



    !    Subroutine tqlrat ( n, d, e2, ierr )
    !
    !! TQLRAT computes all eigenvalues of a real symmetric tridiagonal matrix.
    !
    !  Discussion:
    !
    !    TQLRAT finds the eigenvalues of a symmetric
    !    tridiagonal matrix by the rational QL method.
    !
    !  Modified:
    !
    !    04 March 2018
    !
    !
    !  Reference:
    !
    !    C Reinsch,
    !    Algorithm 464, TQLRAT,
    !    Communications of the ACM,
    !    Volume 16, page 689, 1973.
    !
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input/output, real(kind=dp) D(N).  On input, D contains the diagonal
    !    elements of the matrix.  On output, D contains the eigenvalues in ascending
    !    order.  If an error exit was made, then the eigenvalues are correct
    !    in positions 1 through IERR-1, but may not be the smallest eigenvalues.
    !
    !    Input/output, real(kind=dp) E2(N), contains in positions 2 through N
    !    the squares of the subdiagonal elements of the matrix.  E2(1) is
    !    arbitrary.  On output, E2 has been overwritten by workspace
    !    information.
    !
    !    Output, integer ::IERR, error flag.
    !    0, for no error,
    !    J, if the J-th eigenvalue could not be determined after 30 iterations.
    !
    Subroutine tqlrat ( n, d, e2, ierr )
      integer,                       intent(in)     :: n
      real(kind=dp), dimension(n),   intent(in out) :: d,e2
      integer,                       intent(out)    :: ierr

      real(kind=dp) :: b,c,f,g,h,p,r,s,t
      integer :: i,its,l,l1,m,mml

      ierr = 0

      if ( n == 1 ) return

      do i = 2, n
        e2(i-1) = e2(i)
      end do

      f = 0.0_dp
      t = 0.0_dp
      e2(n) = 0.0_dp

      do l = 1, n

        its = 0
        h = abs ( d(l) ) + sqrt ( e2(l) )

        if ( t <= h ) then

          t = h
          b = abs ( t ) * epsilon ( b )
          c = b * b

        end if
    !
    !  Look for small squared sub-diagonal element.
    !
        do m = l, n
          if ( e2(m) <= c ) then
            exit
          end if
        end do

        if ( m /= l ) then

          do

            if ( 30 <= its ) then
              ierr = l
              return
            end if

            its = its + 1
    !
    !  Form shift.
    !
            l1 = l + 1
            s = sqrt ( e2(l) )
            g = d(l)
            p = ( d(l1) - g ) / ( 2.0_dp * s )
            r = pythag ( p, 1.0_dp )
            d(l) = s / ( p + sign ( r, p ) )
            h = g - d(l)
            d(l1:n) = d(l1:n) - h
            f = f + h
    !
    !  Rational QL transformation.
    !
            g = d(m)
            if ( g == 0.0_dp ) then
              g = b
            end if

            h = g
            s = 0.0_dp
            mml = m - l

            do i = m - 1, l, -1
              p = g * h
              r = p + e2(i)
              e2(i+1) = s * r
              s = e2(i) / r
              d(i+1) = h + s * ( h + d(i) )
              g = d(i) - e2(i) / g
              if ( g == 0.0_dp ) then
                g = b
              end if
              h = g * p / r
            end do

            e2(l) = s * g
            d(l) = h
    !
    !  Guard against underflow in convergence test.
    !
            if ( h == 0.0_dp ) then
              exit
            end if

            if ( abs ( e2(l) ) <= abs ( c / h ) ) then
              exit
            end if

            e2(l) = h * e2(l)

            if ( e2(l) == 0.0_dp ) then
              exit
            end if

          end do

        end if

        p = d(l) + f
    !
    !  Order the eigenvalues.
    !
        do i = l, 1, -1
          if ( i == 1 ) then
            d(i) = p
            exit
          else if ( d(i-1) <= p ) then
            d(i) = p
            exit
          end if
          d(i) = d(i-1)
        end do

      end do

    End Subroutine tqlrat



    !    Subroutine trbak1 ( n, a, e, m, z )
    !
    !! TRBAK1 determines eigenvectors by undoing the TRED1 transformation.
    !
    !  Discussion:
    !
    !    TRBAK1 forms the eigenvectors of a real symmetric
    !    matrix by back transforming those of the corresponding
    !    symmetric tridiagonal matrix determined by TRED1.
    !
    !
    !  Modified:
    !
    !    04 March 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input, real(kind=dp) A(N,N), contains information about the orthogonal
    !    transformations used in the reduction by TRED1 in its strict lower
    !    triangle.
    !
    !    Input, real(kind=dp) E(N), the subdiagonal elements of the tridiagonal
    !    matrix in E(2:N).  E(1) is arbitrary.
    !
    !    Input, integer ::M, the number of eigenvectors to be back
    !    transformed.
    !
    !    Input/output, real(kind=dp) Z(N,M).  On input, the eigenvectors to be
    !    back transformed.  On output, the transformed eigenvectors.
    !
    Subroutine trbak1 ( n, a, e, m, z )
      integer,                       intent(in)     :: n,m
      real(kind=dp), dimension(n,n), intent(in)     :: a
      real(kind=dp), dimension(n),   intent(in)     :: e
      real(kind=dp), dimension(n,m), intent(in out) :: z

      integer :: i,j
      real(kind=dp) :: s

      if ( m <= 0 ) return

      if ( n <= 1 )  return

      do i = 2, n

        if ( e(i) /= 0.0_dp ) then

          do j = 1, m
            s = dot_product ( a(i,1:i-1), z(1:i-1,j) )
            s = ( s / a(i,i-1) ) / e(i)
            z(1:i-1,j) = z(1:i-1,j) + s * a(i,1:i-1)
          end do

        end if

      end do

    End Subroutine trbak1



    !    Subroutine trbak3 ( n, nv, a, m, z )
    !
    !! TRBAK3 determines eigenvectors by undoing the TRED3 transformation.
    !
    !  Discussion:
    !
    !    TRBAK3 forms the eigenvectors of a real symmetric
    !    matrix by back transforming those of the corresponding
    !    symmetric tridiagonal matrix determined by TRED3.
    !
    !
    !  Modified:
    !
    !    28 January 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input, integer ::NV, the dimension of the array paramater A,
    !    which must be at least N*(N+1)/2.
    !
    !    Input, real(kind=dp) A(NV), information about the orthogonal
    !    transformations used in the reduction by TRED3.
    !
    !    Input, integer ::M, the number of eigenvectors to be back
    !    transformed.
    !
    !    Input/output, real(kind=dp) Z(N,M).  On input, the eigenvectors to be
    !    back transformed.  On output, the transformed eigenvectors.
    !
    Subroutine trbak3 ( n, nv, a, m, z )
      integer,                       intent(in)     :: n,nv,m
      real(kind=dp), dimension(nv),  intent(in)     :: a
      real(kind=dp), dimension(n,m), intent(in out) :: z

      real(kind=dp) :: h,s
      integer       :: i,ik,iz,j,k

      if ( m == 0 ) return

      do i = 2, n

        iz = ( i * ( i - 1 ) ) / 2
        ik = iz + i
        h = a(ik)

        if ( h /= 0.0_dp ) then

          do j = 1, m

            s = 0.0_dp
            ik = iz

            do k = 1, i - 1
              ik = ik + 1
              s = s + a(ik) * z(k,j)
            end do

            s = ( s / h ) / h
            ik = iz

            do k = 1, i - 1
              ik = ik + 1
              z(k,j) = z(k,j) - s * a(ik)
            end do

          end do

        end if

      end do

    End Subroutine trbak3



    !    Subroutine tred1 ( n, a, d, e, e2 )
    !
    !! TRED1 transforms a real symmetric matrix to symmetric tridiagonal form.
    !
    !  Discussion:
    !
    !    TRED1 reduces a real symmetric matrix to a symmetric
    !    tridiagonal matrix using orthogonal similarity transformations.
    !
    !  Modified:
    !
    !    04 March 2018
    !
    !  Reference:
    !
    !    Martin, Reinsch, James Wilkinson,
    !    TRED1,
    !    Numerische Mathematik,
    !    Volume 11, pages 181-195, 1968.
    !
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix A.
    !
    !    Input/output, real(kind=dp) A(N,N), on input, contains the real
    !    symmetric matrix.  Only the lower triangle of the matrix need be supplied.
    !    On output, A contains information about the orthogonal transformations
    !    used in the reduction in its strict lower triangle.
    !    The full upper triangle of A is unaltered.
    !
    !    Output, real(kind=dp) D(N), contains the diagonal elements of the
    !    tridiagonal matrix.
    !
    !    Output, real(kind=dp) E(N), contains the subdiagonal elements of the
    !    tridiagonal matrix in its last N-1 positions.  E(1) is set to zero.
    !
    !    Output, real(kind=dp) E2(N), contains the squares of the corresponding
    !    elements of E.  E2 may coincide with E if the squares are not needed.
    !
    Subroutine tred1 ( n, a, d, e, e2 )
      integer,                       intent(in)     :: n
      real(kind=dp), dimension(n,n), intent(in out) :: a
      real(kind=dp), dimension(n),   intent(out)    :: d,e,e2

      real(kind=dp) :: f,g,h,scal
      integer :: i,j,k,l

      d(1:n) = a(n,1:n)

      do i = 1, n
        a(n,i) = a(i,i)
      end do

      do i = n, 1, -1

        l = i - 1
        h = 0.0_dp
    !
    !  scal row.
    !
        scal = sum ( abs ( d(1:l) ) )

        if ( scal == 0.0_dp ) then

          do j = 1, l
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = 0.0_dp
          end do

          e(i) = 0.0_dp
          e2(i) = 0.0_dp

          cycle

        end if

        d(1:l) = d(1:l) / scal

        do k = 1, l
          h = h + d(k) * d(k)
        end do

        e2(i) = h * scal * scal
        f = d(l)
        g = - sign ( sqrt ( h ), f )
        e(i) = scal * g
        h = h - f * g
        d(l) = f - g

        if ( 1 <= l ) then
    !
    !  Form A * U.
    !
          e(1:l) = 0.0_dp

          do j = 1, l

            f = d(j)
            g = e(j) + a(j,j) * f
            do k = j + 1, l
              g = g + a(k,j) * d(k)
              e(k) = e(k) + a(k,j) * f
            end do

            e(j) = g

          end do
    !
    !  Form P.
    !
          f = 0.0_dp

          do j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
          end do

          h = f / ( h + h )
    !
    !  Form Q.
    !
          e(1:l) = e(1:l) - h * d(1:l)
    !
    !  Form reduced A.
    !
          do j = 1, l

            f = d(j)
            g = e(j)

            a(j:l,j) = a(j:l,j) - f * e(j:l) - g * d(j:l)

          end do

        end if

        do j = 1, l
          f = d(j)
          d(j) = a(l,j)
          a(l,j) = a(i,j)
          a(i,j) = f * scal
        end do

      end do

    End Subroutine tred1



    !    Subroutine tred2 ( n, a, d, e, z )
    !
    !! TRED2 transforms a real symmetric matrix to symmetric tridiagonal form.
    !
    !  Discussion:
    !
    !    TRED2 reduces a real symmetric matrix to a
    !    symmetric tridiagonal matrix using and accumulating
    !    orthogonal similarity transformations.
    !
    !    A and Z may coincide, in which case a single storage area is used
    !    for the input of A and the output of Z.
    !
    !
    !  Modified:
    !
    !    04 March 2018
    !
    !
    !  Reference:
    !
    !    Martin, Reinsch, James Wilkinson,
    !    TRED2,
    !    Numerische Mathematik,
    !    Volume 11, pages 181-195, 1968.
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input, real(kind=dp) A(N,N), the real symmetric input matrix.  Only the
    !    lower triangle of the matrix need be supplied.
    !
    !    Output, real(kind=dp) D(N), the diagonal elements of the tridiagonal
    !    matrix.
    !
    !    Output, real(kind=dp) E(N), contains the subdiagonal elements of the
    !    tridiagonal matrix in E(2:N).  E(1) is set to zero.
    !
    !    Output, real(kind=dp) Z(N,N), the orthogonal transformation matrix
    !    produced in the reduction.
    !
    Subroutine tred2 ( n, a, d, e, z )
      integer,                       intent(in)     :: n
      real(kind=dp), dimension(n,n), intent(in)     :: a
      real(kind=dp), dimension(n),   intent(out)    :: d,e
      real(kind=dp), dimension(n,n), intent(out)    :: z

      real(kind=dp) :: f,g,h,hh,scal
      integer :: i,j,k,l

      do i = 1, n
        z(i:n,i) = a(i:n,i)
      end do

      d(1:n) = a(n,1:n)

      do i = n, 2, -1

        l = i - 1
        h = 0.0_dp
    !
    !  scal row.
    !
        scal = sum ( abs ( d(1:l) ) )

        if ( scal == 0.0_dp ) then

          e(i) = d(l)

          do j = 1, l
            d(j) = z(l,j)
            z(i,j) = 0.0_dp
            z(j,i) = 0.0_dp
          end do

          d(i) = 0.0_dp

          cycle

        end if

        d(1:l) = d(1:l) / scal

        h = h + dot_product ( d(1:l), d(1:l) )

        f = d(l)
        g = - sign ( sqrt ( h ), f )
        e(i) = scal * g
        h = h - f * g
        d(l) = f - g
    !
    !  Form A*U.
    !
        e(1:l) = 0.0_dp

        do j = 1, l

          f = d(j)
          z(j,i) = f
          g = e(j) + z(j,j) * f

          do k = j + 1, l
            g = g + z(k,j) * d(k)
            e(k) = e(k) + z(k,j) * f
          end do

          e(j) = g

        end do
    !
    !  Form P.
    !
        e(1:l) = e(1:l) / h

        f = dot_product ( e(1:l), d(1:l) )

        hh = 0.5_dp * f / h
    !
    !  Form Q.
    !
        e(1:l) = e(1:l) - hh * d(1:l)
    !
    !  Form reduced A.
    !
        do j = 1, l

          f = d(j)
          g = e(j)

          z(j:l,j) = z(j:l,j) - f * e(j:l) - g * d(j:l)

          d(j) = z(l,j)
          z(i,j) = 0.0_dp

        end do

        d(i) = h

      end do
    !
    !  Accumulation of transformation matrices.
    !
      do i = 2, n

    !   l = i - 1
        z(n,i-1) = z(i-1,i-1)
        z(i-1,i-1) = 1.0_dp
        h = d(i)

        if ( h /= 0.0_dp ) then

          d(1:i-1) = z(1:i-1,i) / h

          do j = 1, i - 1

            g = dot_product ( z(1:i-1,i), z(1:i-1,j) )

            do k = 1, i - 1
              z(k,j) = z(k,j) - g * d(k)
            end do

          end do

        end if

        z(1:i-1,i) = 0.0_dp

      end do

      d(1:n) = z(n,1:n)

      z(n,1:n-1) = 0.0_dp
      z(n,n) = 1.0_dp

      e(1) = 0.0_dp

      return
    End  Subroutine tred2



    !    Subroutine tred3 ( n, nv, a, d, e, e2 )
    !
    !! TRED3: transform real symmetric packed matrix to symmetric tridiagonal form.
    !
    !  Discussion:
    !
    !    TRED3 reduces a real symmetric matrix, stored as
    !    a one-dimensional array, to a symmetric tridiagonal matrix
    !    using orthogonal similarity transformations.
    !
    !
    !  Modified:
    !
    !    04 March 2018
    !
    !
    !  Reference:
    !
    !    Martin, Reinsch, James Wilkinson,
    !    TRED3,
    !    Numerische Mathematik,
    !    Volume 11, pages 181-195, 1968.
    !
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input, integer ::NV, the dimension of A, which must be at least
    !    (N*(N+1))/2.
    !
    !    Input/output, real(kind=dp) A(NV).  On input, the lower triangle of
    !    the real symmetric matrix, stored row-wise.  On output, information about
    !    the orthogonal transformations used in the reduction.
    !
    !    Output, real(kind=dp) D(N), the diagonal elements of the tridiagonal
    !    matrix.
    !
    !    Output, real(kind=dp) E(N), the subdiagonal elements of the tridiagonal
    !    matrix in E(2:N).  E(1) is set to zero.
    !
    !    Output, real(kind=dp) E2(N),  the squares of the corresponding
    !    elements of E.  E2 may coincide with E if the squares are not needed.
    !
    Subroutine tred3 ( n, nv, a, d, e, e2 )
      integer,                       intent(in)     :: n,nv
      real(kind=dp), dimension(nv),  intent(in out) :: a
      real(kind=dp), dimension(n),   intent(out)    :: d,e,e2

      real(kind=dp) :: f,g,h,hh,scal
      integer :: i,iz,j,jk,k

      do i = n, 1, -1

        iz = ( i * ( i - 1 ) ) / 2
        h = 0.0_dp
        scal = 0.0_dp
    !
    !  scal row.
    !
        do k = 1, i - 1
          iz = iz + 1
          d(k) = a(iz)
          scal = scal + abs ( d(k) )
        end do

        if ( scal == 0.0_dp ) then
          e(i) = 0.0_dp
          e2(i) = 0.0_dp
          d(i) = a(iz+1)
          a(iz+1) = scal * sqrt ( h )
          cycle
        end if

        do k = 1, i - 1
          d(k) = d(k) / scal
          h = h + d(k) * d(k)
        end do

        e2(i) = scal * scal * h
        f = d(i-1)
        g = - sign ( sqrt ( h ), f )
        e(i) = scal * g
        h = h - f * g
        d(i-1) = f - g
        a(iz) = scal * d(i-1)

        if ( i == 2 ) then
          d(i) = a(iz+1)
          a(iz+1) = scal * sqrt ( h )
          cycle
        end if

        jk = 1

        do j = 1, i - 1

          f = d(j)
          g = 0.0_dp

          do k = 1, j - 1
            g = g + a(jk) * d(k)
            e(k) = e(k) + a(jk) * f
            jk = jk + 1
          end do

          e(j) = g + a(jk) * f
          jk = jk + 1

        end do
    !
    !  Form P.
    !
        e(1:i-1) = e(1:i-1) / h
        f = dot_product ( e(1:i-1), d(1:i-1) )
        hh = f / ( h + h )
    !
    !  Form Q.
    !
        e(1:i-1) = e(1:i-1) - hh * d(1:i-1)
        jk = 1
    !
    !  Form reduced A.
    !
        do j = 1, i - 1
          f = d(j)
          g = e(j)
          do k = 1, j
            a(jk) = a(jk) - f * e(k) - g * d(k)
            jk = jk + 1
          end do
        end do

        d(i) = a(iz+1)
        a(iz+1) = scal * sqrt ( h )

      end do

      return
    End  Subroutine tred3



    !    Subroutine tridib ( n, eps1, d, e, e2, lb, ub, m11, m, w, ind, ierr )
    !
    !! TRIDIB computes some eigenvalues of a real symmetric tridiagonal matrix.
    !
    !  Discussion:
    !
    !    TRIDIB finds those eigenvalues of a tridiagonal symmetric matrix between
    !    specified boundary indices, using bisection.
    !
    !
    !  Modified:
    !
    !    12 March 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input/output, real(kind=dp) EPS1.  On input, an absolute error
    !    tolerance for the computed eigenvalues.  It should be chosen commensurate
    !    with relative perturbations in the matrix elements of the order of the
    !    relative machine precision.  If the input EPS1 is non-positive, it
    !    is reset for each submatrix to a default value, namely, minus the
    !    product of the relative machine precision and the 1-norm of the submatrix.
    !
    !    Input, real(kind=dp) D(N), the diagonal elements of the input matrix.
    !
    !    Input, real(kind=dp) E(N), the subdiagonal elements of the input matrix
    !    in E(2:N).  E(1) is arbitrary.
    !
    !    Input/output, real(kind=dp) E2(N).  On input, the squares of the
    !    corresponding elements of E.  E2(1) is arbitrary.  On output, elements of
    !    E2 corresponding to elements of E regarded as negligible, have been
    !    replaced by zero, causing the matrix to split into a direct sum of
    !    submatrices.  E2(1) is also set to zero.
    !
    !    Output, real(kind=dp) LB, UB, define an interval containing exactly
    !    the desired eigenvalues.
    !
    !    Input, integer ::M11, the lower boundary index for the desired
    !    eigenvalues.
    !
    !    Input, integer ::M, the number of eigenvalues desired.  The
    !    upper boundary index M22 is then obtained as M22 = M11 + M - 1.
    !
    !    Output, real(kind=dp) W(M), the eigenvalues between indices M11 and M22
    !    in ascending order.
    !
    !    Output, integer ::IND(M), the submatrix indices associated with
    !    the corresponding eigenvalues in W: 1 for eigenvalues belonging to the
    !    first submatrix from the top, 2 for those belonging to the second
    !    submatrix, and so on.
    !
    !    Output, integer ::IERR, error flag.
    !    0, for normal return,
    !    3*N+1, if multiple eigenvalues at index M11 make unique selection
    !      impossible,
    !    3*N+2, if multiple eigenvalues at index M22 make unique selection
    !      impossible.
    !
    Subroutine tridib ( n, eps1, d, e, e2, lb, ub, m11, m, w, ind, ierr )
      integer,                       intent(in)     :: n,m,m11
      real(kind=dp),                 intent(in out) :: eps1
      real(kind=dp), dimension(n),   intent(in)     :: d,e
      real(kind=dp), dimension(n),   intent(in out) :: e2
      real(kind=dp),                 intent(out)    :: lb,ub
      real(kind=dp), dimension(m),   intent(out)    :: w
      integer,                       intent(out)    :: ierr

      integer       :: i,isturm,j,k,l,m1,m2,m22,p,q,r,s,tag
      real(kind=dp) :: t1,t2,tst1,tst2,u,v,x0,x1,xu
      integer,       dimension(m) :: ind
      real(kind=dp), dimension(n) :: rv4,rv5

      ierr = 0
      tag = 0
      xu = d(1)
      x0 = d(1)
      s = 0
      u = 0.0_dp
    !
    !  Look for small sub-diagonal entries and determine an
    !  interval containing all the eigenvalues.
    !
      do i = 1, n

        x1 = u

        if ( i == n ) then
          u = 0.0_dp
        else
          u = abs ( e(i+1) )
        end if

        xu = min ( xu, d(i) - ( x1 + u ) )
        x0 = max ( x0, d(i) + ( x1 + u ) )

        if ( 1 < i ) then
          tst1 = abs ( d(i) ) + abs ( d(i-1) )
          tst2 = tst1 + abs ( e(i) )
          if ( tst2 <= tst1 ) then
            e2(i) = 0.0_dp
          end if
        else
          e2(i) = 0.0_dp
        end if

      end do

      x1 = n
      x1 = x1 * max ( abs ( xu ), abs ( x0 ) ) * epsilon ( x1 )
      xu = xu - x1
      t1 = xu
      x0 = x0 + x1
      t2 = x0
    !
    !  Determine an interval containing exactly the desired eigenvalues.
    !
      p = 1
      q = n
      m1 = m11 - 1
      m22 = m1 + m

      if ( m1 /= 0 .or. m /= n ) then

        if ( m1 == 0 .and. m /= n ) then
           x0 = t2
           isturm = 2
        else
          isturm = 1
        end if

        do

          v = x1
          x1 = xu + ( x0 - xu ) * 0.5_dp

          if ( x1 == v ) then
            ierr = 3 * n + isturm
            lb = t1
            ub = t2
            return
          end if

          s = sturm_sequence ( d, e, e2, n, p, q, x1 )

          if ( isturm == 1 ) then

            if ( s < m1 ) then
              xu = x1
              cycle
            else if ( m1 < s ) then
              x0 = x1
              cycle
            end if

            xu = x1
            t1 = x1
            m22 = m1 + m

            if ( m22 /= n ) then
              x0 = t2
              isturm = 2
              cycle
            end if

          else

            if ( s < m22 ) then
              xu = x1
              cycle
            else if ( m22 < s ) then
              x0 = x1
              cycle
            end if

            t2 = x1

          end if

          exit

        end do

      end if

      q = 0
      r = 0
    !
    !  Establish and process next submatrix, refining interval by the
    !  Gerschgorin bounds.
    !
      do

        if ( r == m ) then
          lb = t1
          ub = t2
          return
        end if

        tag = tag + 1
        p = q + 1
        xu = d(p)
        x0 = d(p)
        u = 0.0_dp

        do q = p, n

          x1 = u
          u = 0.0_dp
          v = 0.0_dp

          if ( q < n ) then
            u = abs ( e(q+1) )
            v = e2(q+1)
          end if

          xu = min ( d(q) - ( x1 + u ), xu )
          x0 = max ( d(q) + ( x1 + u ), x0 )

          if ( v == 0.0_dp ) then
            exit
          end if

        end do

        x1 = max ( abs ( xu ), abs ( x0 ) ) * epsilon ( x1 )

        if ( eps1 <= 0.0_dp ) then
          eps1 = - x1
        end if
    !
    !  Check for isolated root within interval.
    !
        if ( p == q ) then

          if ( d(p) < t1 .or. t2 <= d(p) ) then
            if ( n <= q ) then
              exit
            end if
            cycle
          end if
          m1 = p
          m2 = p
          rv5(p) = d(p)

        else

          x1 = x1 * ( q - p + 1 )
          lb = max ( t1, xu - x1 )
          ub = min ( t2, x0 + x1 )

          x1 = lb
          s = sturm_sequence ( d, e, e2, n, p, q, x1 )
          m1 = s + 1

          x1 = ub
          s = sturm_sequence ( d, e, e2, n, p, q, x1 )
          m2 = s

          if ( m2 < m1 ) then
            if ( n <= q ) then
              exit
            end if
            cycle
          end if
    !
    !  Find roots by bisection.
    !
          x0 = ub
          isturm = 5

          rv5(m1:m2) = ub
          rv4(m1:m2) = lb
    !
    !  Loop for the K-th eigenvalue.
    !
          k = m2

          do

            xu = lb

            do i = k, m1, -1
              if ( xu < rv4(i) ) then
                xu = rv4(i)
                exit
              end if

            end do

            x0 = min ( x0, rv5(k) )
    !
    !  Next bisection step.
    !
            do

              x1 = ( xu + x0 ) * 0.5_dp

              if ( ( x0 - xu ) <= abs ( eps1 ) ) then
                exit
              end if

              tst1 = 2.0_dp * ( abs ( xu ) + abs ( x0 ) )
              tst2 = tst1 + ( x0 - xu )

              if ( tst2 == tst1 ) then
                exit
              end if

              s = sturm_sequence ( d, e, e2, n, p, q, x1 )
    !
    !  Refine intervals.
    !
              if ( k <= s ) then
                x0 = x1
              else
                xu = x1
                if ( m1 <= s ) then
                  rv4(s+1) = x1
                  rv5(s) = min ( rv5(s), x1 )
                else
                  rv4(m1) = x1
                end if
              end if

            end do
    !
    !  K-th eigenvalue found.
    !
            rv5(k) = x1
            k = k - 1

            if ( k < m1 ) then
              exit
            end if

          end do

        end if
    !
    !  Order eigenvalues tagged with their submatrix associations.
    !
        s = r
        r = r + m2 - m1 + 1
        j = 1
        k = m1

        do l = 1, r

          if ( j <= s ) then

            if ( m2 < k ) then
              exit
            end if

            if ( w(l) <= rv5(k) ) then
              j = j + 1
              cycle
            end if

            do i = l + s - j, l, -1
              w(i+1) = w(i)
              ind(i+1) = ind(i)
            end do

          end if

          w(l) = rv5(k)
          ind(l) = tag
          k = k + 1

        end do

        if ( n <= q ) then
          exit
        end if

      end do

      lb = t1
      ub = t2

    End Subroutine tridib

    !  Subroutine Tsturm ( n, eps1, d, e, e2, t1, t2, mm, m, w, z, ierr )
    !
    !! TSTURM computes some eigenvalues/vectors, real symmetric tridiagonal matrix.
    !
    !  Discussion:
    !
    !    TSTURM finds those eigenvalues of a tridiagonal symmetric matrix which
    !    lie in a specified interval and their associated eigenvectors, using
    !    bisection and inverse iteration.
    !
    !  Modified:
    !
    !    04 March 2018
    !
    !  Arguments:
    !
    !    Input, integer ::N, the order of the matrix.
    !
    !    Input/output, real(kind=dp) EPS1.  On input, an absolute error
    !    tolerance for the computed eigenvalues.  It should be chosen commensurate
    !    with relative perturbations in the matrix elements of the order of the
    !    relative machine precision.  If the input EPS1 is non-positive, it
    !    is reset for each submatrix to a default value, namely, minus the
    !    product of the relative machine precision and the 1-norm of the submatrix.
    !
    !    Input, real(kind=dp) D(N), the diagonal elements of the input matrix.
    !
    !    Input, real(kind=dp) E(N), the subdiagonal elements of the input matrix
    !    in E(2:N).  E(1) is arbitrary.
    !
    !    Input/output, real(kind=dp) E2(N).  On input, the squares of the
    !    corresponding elements of E.  E2(1) is arbitrary.  On output, elements of
    !    E2 corresponding to elements of E regarded as negligible have been
    !    replaced by zero, causing the matrix to split into a direct sum of
    !    submatrices.  E2(1) is also set to zero.
    !
    !    Input, real(kind=dp) T1, T2, define the interval to be searched
    !    for eigenvalues.  If T1 is not less than T2, no eigenvalues will be found.
    !
    !    Input, integer ::MM, an upper bound for the number of
    !    eigenvalues in the interval.  If more than MM eigenvalues are determined
    !    to lie in the interval, an error return is made with no values or vectors
    !    found.
    !
    !    Output, integer ::M, the number of eigenvalues determined to lie
    !    in (T1,T2).
    !
    !    Output, real(kind=dp) W(M), the eigenvalues in ascending order if the
    !    matrix does not split.  If the matrix splits, the eigenvalues are in
    !    ascending order for each submatrix.  If a vector error exit is made, W
    !    contains those values already found.
    !
    !    Output, real(kind=dp) Z(N,MM), the associated set of orthonormal
    !    eigenvectors.  If an error exit is made, Z contains those vectors already
    !    found.
    !
    !    Output, integer ::IERR, error flag.
    !    0, normal return.
    !    3*N+1, if M exceeds MM.
    !    4*N+R, if the eigenvector corresponding to the R-th
    !      eigenvalue fails to converge in 5 iterations.
    !
    Subroutine Tsturm ( n, eps1, d, e, e2, t1, t2, mm, m, w, z, ierr )
      integer,                       intent(in)     :: n,mm
      real(kind=dp),                 intent(in out) :: eps1
      real(kind=dp), dimension(n),   intent(in)     :: d,e
      real(kind=dp), dimension(n),   intent(in out) :: e2
      real(kind=dp),                 intent(in)     :: t1,t2
      integer,                       intent(out)    :: m
      real(kind=dp), dimension(mm),  intent(out)    :: w
      real(kind=dp), dimension(n,mm),intent(out)    :: z
      integer,                       intent(out)    :: ierr

      integer       :: group,i,ip,its,j,k,m1,m2,p,q,r,s
      real(kind=dp) :: eps2,eps3,eps4,lb,norm,tst1,tst2,u,ub,uk,v,x0,x1,xu
      real(kind=dp), dimension(n) :: rv1 ,rv2,rv3,rv4,rv5,rv6

      ierr = 0
      s = 0
      lb = t1
      ub = t2
      !
      !  Look for small sub-diagonal entries.
      !
      e2(1) = 0.0_dp
      do i = 2, n

        tst1 = abs ( d(i) ) + abs ( d(i-1) )
        tst2 = tst1 + abs ( e(i) )

        if ( tst2 <= tst1 ) then
          e2(i) = 0.0_dp
        end if

      end do
      !
      !  Determine the number of eigenvalues in the interval.
      !
      p = 1
      q = n
      x1 = ub
      s = sturm_sequence ( d, e, e2, n, p, q, x1 )
      m = s

      x1 = lb
      s = sturm_sequence ( d, e, e2, n, p, q, x1 )
      m = m - s

      if ( mm < m ) then
        ierr = 3 * n + 1
        return
      end if

      q = 0
      r = 0
      !
      !  Establish and process next submatrix, refining interval by the
      !  Gerschgorin bounds.
      !
      do while ( q < n )

        if ( r == m ) then
          return
        end if

        p = q + 1
        xu = d(p)
        x0 = d(p)
        u = 0.0_dp

        do q = p, n

          x1 = u
          u = 0.0_dp
          v = 0.0_dp

          if ( q < n ) then
            u = abs ( e(q+1) )
            v = e2(q+1)
          end if

          xu = min ( d(q) - ( x1 + u ), xu )
          x0 = max ( d(q) + ( x1 + u ), x0 )

          if ( v == 0.0_dp ) then
            exit
          end if

        end do

        x1 = max ( abs ( xu ), abs ( x0 ) ) * epsilon ( x1 )

        if ( eps1 <= 0.0_dp ) then
          eps1 = - x1
        end if
        !
        !  Check for isolated root within interval.
        !
        if ( p == q ) then

          if ( d(p) < t1 .or. t2 <= d(p) ) cycle

          r = r + 1

          z(1:n,r) = 0.0_dp

          w(r) = d(p)
          z(p,r) = 1.0_dp
          cycle

        end if

        u = q - p + 1
        x1 = u * x1
        lb = max ( t1, xu - x1 )
        ub = min ( t2, x0 + x1 )
        x1 = lb
        s = sturm_sequence ( d, e, e2, n, p, q, x1 )
        m1 = s + 1

        x1 = ub
        s = sturm_sequence ( d, e, e2, n, p, q, x1 )
        m2 = s

        if ( m2 < m1 ) cycle
    !
    !  Find roots by bisection.
    !
        x0 = ub
        rv5(m1:m2) = ub
        rv4(m1:m2) = lb
    !
    !  Loop for K-th eigenvalue.
    !
        k = m2

        do while ( m1 <= k )

          xu = lb

          do i = k, m1, -1

            if ( xu < rv4(i) ) then
              xu = rv4(i)
              exit
            end if

          end do

          x0 = min ( x0, rv5(k) )
    !
    !  Next bisection step.
    !
          do

            x1 = ( xu + x0 ) * 0.5_dp

            if ( ( x0 - xu ) <= abs ( eps1 ) ) exit

            tst1 = 2.0_dp * ( abs ( xu ) + abs ( x0 ) )
            tst2 = tst1 + ( x0 - xu )

            if ( tst2 == tst1 ) exit

            s = sturm_sequence ( d, e, e2, n, p, q, x1 )
    !
    !  Refine intervals.
    !
            if ( k <= s ) then
              x0 = x1
              cycle
            end if

            xu = x1

            if ( m1 <= s ) then
              rv4(s+1) = x1
              if ( x1 < rv5(s) ) then
                rv5(s) = x1
              end if
            else
              rv4(m1) = x1
            end if

          end do
    !
    !  K-th eigenvalue found.
    !
          rv5(k) = x1
          k = k - 1

        end do
    !
    !  Find vectors by inverse iteration.
    !
        norm = abs ( d(p) )
        ip = p + 1

        do i = p + 1, q
          norm = max ( norm, abs ( d(i) ) + abs ( e(i) ) )
        end do
    !
    !  EPS2 is the criterion for grouping,
    !  EPS3 replaces zero pivots and equal roots are modified by eps3,
    !  EPS4 is taken very small to avoid overflow.
    !
        eps2 = 0.001_dp * norm
        eps3 = abs ( norm ) * epsilon ( eps3 )
        uk = q - p + 1
        eps4 = uk * eps3
        uk = eps4 / sqrt ( uk )
        group = 0
        s = p

        do k = m1, m2

          r = r + 1
          its = 1
          w(r) = rv5(k)
          x1 = rv5(k)
    !
    !  Look for close or coincident roots.
    !
          if ( k /= m1 ) then
            if ( eps2 <= x1 - x0 ) then
              group = -1
            end if
            group = group + 1
            if ( x1 <= x0 ) then
              x1 = x0 + eps3
            end if
          end if
    !
    !  Elimination with interchanges and initialization of vector.
    !
          v = 0.0_dp

          do i = p, q

            rv6(i) = uk

            if ( i /= p ) then

              if ( abs ( u ) <= abs ( e(i) ) ) then
                xu = u / e(i)
                rv4(i) = xu
                rv1(i-1) = e(i)
                rv2(i-1) = d(i) - x1
                rv3(i-1) = 0.0_dp
                if ( i /= q ) then
                  rv3(i-1) = e(i+1)
                end if
                u = v - xu * rv2(i-1)
                v = - xu * rv3(i-1)
                cycle
              end if

              xu = e(i) / u
              rv4(i) = xu
              rv1(i-1) = u
              rv2(i-1) = v
              rv3(i-1) = 0.0_dp

            end if

            u = d(i) - x1 - xu * v

            if ( i /= q ) then
              v = e(i+1)
            end if

           end do

          if ( u /= 0.0_dp ) then
            rv1(q) = u
          else
            rv1(q) = eps3
          end if

          rv2(q) = 0.0_dp
          rv3(q) = 0.0_dp
    !
    !  Back substitution.
    !
          do

            do i = q, p, -1
              rv6(i) = ( rv6(i) - u * rv2(i) - v * rv3(i) ) / rv1(i)
              v = u
              u = rv6(i)
            end do
    !
    !  Orthogonalize with respect to previous members of group.
    !
            do j = r - group, r - 1
              xu = dot_product ( rv6(p:q), z(p:q,j) )
              rv6(p:q) = rv6(p:q) - xu * z(p:q,j)
            end do

            norm = sum ( abs ( rv6(p:q) ) )

            if ( 1.0_dp <= norm ) then
              exit
            end if
    !
    !  Forward substitution.
    !
            if ( its == 5 ) then
              ierr = 4 * n + r
              return
            end if

            if ( norm == 0.0_dp ) then
              rv6(s) = eps4
              s = s + 1
              if ( q < s ) then
                s = p
              end if
            else
              xu = eps4 / norm
              rv6(p:q) = rv6(p:q) * xu
            end if
    !
    !  Elimination operations on next vector iterate.
    !
    !  If rv1(i-1) == e(i), a row interchange was performed earlier in the
    !  triangularization process.
    !
            do i = p + 1, q

              u = rv6(i)

              if ( rv1(i-1) == e(i) ) then
                u = rv6(i-1)
                rv6(i-1) = rv6(i)
              end if

              rv6(i) = u - rv4(i) * rv6(i-1)

            end do

            its = its + 1

          end do
    !
    !  Normalize so that sum of squares is 1 and expand to full order.
    !
          u = 0.0_dp
          do i = p, q
            u = pythag ( u, rv6(i) )
          end do

          xu = 1.0_dp / u

          z(1:n,r) = 0.0_dp
          z(p:q,r) = rv6(p:q) * xu

          x0 = x1

        end do

      end do

    End Subroutine Tsturm

End Module CFML_EisPack