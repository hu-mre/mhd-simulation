module subprog
  implicit none
 contains
   double precision function hmnsnr(dim,carray,ocarray) ! SNR of nth-harmonic signal
     integer, parameter :: n = 8192,m = 50
     integer    i, dim
     double precision seq(n),seq2(n)
     double precision :: sd=0d0,ave=0d0,sdr=800.0d0
     complex(kind(0d0)) carray(m,n)
     complex(kind(0d0)), optional :: ocarray(m,n)
     sd=0d0;ave=0d0;sdr=800.0d0

     seq(1:n) = dreal(carray(dim,1:n))
     ave = sum(seq(1:int(sdr))) / sdr
     do i = 1,int(sdr)
       sd = sd + (seq(i) - ave)**2.0d0
     end do
     sd = sqrt(1.0d0/sdr * sd)

     if(present(ocarray))then
       seq2 = dreal(ocarray(dim,1:n))
       hmnsnr = (maxval(seq2) - minval(seq2)) / (2.0d0*sd)
     else
       hmnsnr = (maxval(seq) - minval(seq)) / (2.0d0*sd)
     endif
   end function hmnsnr

  subroutine lpfilter(f,array,cf)
    use fftcf_int ! FFT
    use fftcb_int ! Inverse FFT
    integer, parameter :: n = 8192, m = 50
    integer i
    complex(kind(0d0)) suqwav(n), lpf(n), gauss(n), f(m,n)
    complex(kind(0d0)), optional :: array(m,n)
    double precision cf

    suqwav = (0.0,0.0)
    lpf    = (0.0,0.0)

    do i = 1,cf
      suqwav(i)     = 1.0
      suqwav(n-i+1) = 1.0
    end do
    do i = 1,n/2
      gauss(i) = dexp(((-1.0d0)*dble(i-1)**2.0d0) / (2.0d0*(5.0d0**2.0d0)))
      gauss(n-i+1) = gauss(i)
    end do

    call fftcf(n,suqwav,suqwav)
    call fftcf(n,gauss,gauss)

    do i=1,n
      lpf(i) = suqwav(i)/dble(n) * gauss(i)/dble(n)
    end do
    call fftcb(n,lpf,lpf)
    lpf = lpf/maxval(dble(lpf))

    if(present(array))then
      array(1,1:n) = f(1,1:n) * lpf(1:n)
    else
      f(1,1:n) = f(1,1:n) * lpf(1:n)
    end if
  end subroutine lpfilter

  subroutine bwfilter(harmonic,bwcf)
    use fftcf_int !FFT
    use fftcb_int !Inverse FFT
    use const_int !Const
    integer, parameter :: n=8192, m=50
    integer i,j,k
    double precision bwcf,pi
    complex(kind(0d0)) s(6), bwf(n) ,harmonic(m,n)
    complex(kind(0d0)) :: z = (0.0d0, 1.0d0)
    pi = const('PI')

    do k = 1,6
      s(k) = bwcf * zexp(z*pi*(2.0d0*dble(k)+6.0d0-1.0d0)/(2.0d0*6.0d0))
    end do

    do i = 1,m
      call fftcf(n,harmonic(i,1:n),harmonic(i,1:n))
      do j = 1,n/2
        bwf(j) = bwcf / ((z*dble(j)-s(1))*(z*dble(j)-s(2))*(z*dble(j)-s(3))*&
                        &(z*dble(j)-s(4))*(z*dble(j)-s(5))*(z*dble(j)-s(6)))
        bwf(n-j+1) = bwf(j)
      end do
    end do

    bwf = abs(bwf)
    bwf = bwf/maxval(dreal(bwf))

    do i = 1,n
      harmonic(1:m,i) = bwf(i) * harmonic(1:m,i)
      harmonic(1:m,n-i+1) = bwf(n-i+1) * harmonic(1:m,n-i+1)
    end do

    do i = 1,m
      call fftcb(n,harmonic(i,1:n),harmonic(i,1:n))
    end do

  end subroutine bwfilter
end module subprog


program mhdcalc
  !---------------------------
  ! Function and Subroutine
  ! IMSL Fortran Library
  !
  use subprog
  use const_int !Const
  use fftcf_int !FFT
  use bsjs_int !1st kind Bessel Function
  use fftcb_int !Inverse FFT
  use cerfe_int !complementary error function

  !--------------------------------
  !Veriable Declaration
  !
  implicit none
  integer, parameter ::  n = 8192, p = 128, m=50
  integer                i, j, t, seedsize
  integer,allocatable::  seed(:)
  double precision       pi, b, spec(n,p) ,nspec(n,p), x(p),alpha,snr1,cf,&
                        &bm, width, snr, bs(0:m+1), noise,temp(n),temp2(5,n),amplitude
  double precision ::    bmo = 0.060d0, dcf = 1.0d0 ,bwcf = 2604.0d0/10.0d0, &
                        &bo(3), sigma(3), ganmma(3),&
                        &hs=5.0d0,&
                        &wl(3) = (/0.0347d0,0.0337d0,0.0362d0/),&
                        &wg(3) = (/0.0449d0,0.0446d0,0.044d0/)
  complex(kind(0d0))     f_num(n), f_denom(n), f(m,n), array(m,n),&
                        &seq(n), nharmonic(m,n), harmonic(m,n), mharmonic(m,n), dn(m,n)
  complex(kind(0d0)) ::  z = (0.0d0,1.0d0)
  character(len=80) :: dummy

  !----------------------------------------------------------------------
  ! Calculate constants
  !
  pi = const('PI')
  ganmma = dsqrt(3.0d0)*wl/2.0d0
  sigma = wg/(2.0d0*sqrt(2.0d0*dlog(2.0d0)))
  bo = (/-1.61203d0, 0.00006400002d0, 1.724d0/) + 2.50d0

  !----------------------------------------------------------------------
  ! Open Files
  !
  open(11,file='./gradput/gamma-snr2.csv') ! Relation between noise and snr1
  open(12,file='./gradput/noisedharmonic.csv') ! Harmonics
  open(13,file='./gradput/non-noisedharmonics.csv') ! Harmonics (noise free)

  ! Type Modulation amplitude
  write(*,*)'alpha (modulation ratio):'
  read(*,*)alpha
  noise = 1.0d0
  bm = bmo*alpha

  !----------------------------------------------------------------------
  ! Get and put seed
  !
  call random_seed(size=seedsize)  ! get the seed size
  allocate(seed(seedsize))         ! allocate the memory space for containing the seed
  call random_seed(get=seed)
  call random_seed(put=seed)
  write(*,*)'seed:',seed

  !----------------------------------------------------------------------
  !spectroscopy
  !
  do  i = 1,n
    if(i<=2730)               t=1
    if(i>2730 .and. i<=5460)  t=2
    if(i>5460 .and. i<=8192)  t=3
    call random_number(x)
    do j = 1,p
      b         = (dble(i)*hs)/dble(n) + (bm/2.0d0) * dcos(dble(j) * 2.0d0 * pi / dble(p)) ! magnetic field
      nspec(i,j) =dreal(cerfe( ((b-bo(t)) + z*ganmma(t)) / (sigma(t)*dsqrt(2.0d0))))&
                &/ (sigma(t)*dsqrt(2.0d0*pi))
      spec(i,j) = nspec(i,j)   + noise * (x(j) - 0.50d0) ! Voigt function + noise
    end do
  end do

  !--------------------------------------------------------------------
  ! Phase Sensitive Detection
  !
  do i = 1,n
    seq(1:p) = spec(i,1:p) ! put the spectral data in array 'seq'
    call fftcf(p,seq,seq)
    do j = 1,m
      harmonic(j,i) = dreal(seq(j+1))/dble(p) ! put the Fourier coefficients in array 'harmonic'
    end do
  end do

  do i = 1,n
    seq(1:p) = nspec(i,1:p) ! put the spectral data in array 'seq'
    call fftcf(p,seq,seq)
    do j = 1,m
      nharmonic(j,i) = dreal(seq(j+1))/dble(p)
    end do
  end do

  !--------------------------------------------------------------------
  !Butterworth Filter
  !
  call bwfilter(harmonic,bwcf)
  call bwfilter(nharmonic,bwcf)

  do i = 1,n
    write(12,*)(dreal(harmonic(j,i)),',',j = 1,m)
    write(13,*)(dreal(nharmonic(j,i)),',',j = 1,m)
  end do

  !--------------------------------------------------------------------
  ! Output SNR(n=1)
  !
  write(11,*) alpha,',',hmnsnr(1,harmonic)



  close(11)
  close(12)
  close(13)

end program mhdcalc
