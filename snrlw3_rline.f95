module subprog
  implicit none
 contains
   double precision function hmnsnr(dim,carray,ocarray) ! SNR of n-th harmonic signal
     use rline_int
     integer, parameter :: n = 8192, m = 50
     integer    i, dim
     double precision seq(n),seq2(n),x(5201:5600),y(5201:5600),a0,a1
     double precision :: sd=0d0,ave=0d0,sdr=400.0d0
     complex(kind(0d0)) carray(m,n)
     complex(kind(0d0)), optional :: ocarray(m,n)
     sd=0d0;ave=0d0;sdr=400.0d0

     seq(1:n) = dreal(carray(dim,1:n))

     do i = 5201,5600
       x(i) = dble(i)
       y(i) = seq(i)
     end do

     call rline(x,y,a0,a1)

     do i =5201,5600
       y(i) = a1*i + a0
       sd = sd + (seq(i) - y(i))**2.0d0
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
  use const_int ! Const
  use fftcf_int ! FFT
  use bsjs_int  ! 1st kind Bessel Function
  use fftcb_int ! Inverse FFT
  use cerfe_int ! complementary error function

  !--------------------------------
  ! Variable Declaration
  !
  implicit none
  integer, parameter ::  n = 8192, p = 128, m=50
  integer                i, j, t, seedsize
  integer,allocatable::  seed(:)
  double precision       pi, b, spec(n,p) ,nspec(n,p), x(p),alpha,snr1,&
                        &bm, width, snr, bs(0:m+1), noise,temp(n),temp2(5,n),amplitude
  double precision ::    bmo = 0.060d0, dcf = 1.0d0 ,bwcf = 2604.0d0/10.0d0, &
                        &bo(3), sigma(3), gamma(3),&
                        &hs=5.0d0,cf=80,&
                        &wl(3) = (/0.0347d0,0.0337d0,0.0362d0/),&
                        &wg(3) = (/0.0449d0,0.0446d0,0.044d0/)
  complex(kind(0d0))     f_num(n), f_denom(n), f(m,n), array(m,n),&
                        &seq(n), nharmonic(m,n), harmonic(m,n), mharmonic(m,n), dn(m,n)
  complex(kind(0d0)) ::  z = (0.0d0,1.0d0)
  character(len=80) :: dummy

  !----------------------------------------------------------------------
  ! calculate constants
  !
  pi = const('PI')
  gamma = dsqrt(3.0d0)*wl/2.0d0
  sigma = wg/(2.0d0*sqrt(2.0d0*dlog(2.0d0)))
  bo = (/-1.61203d0, 0.00006400002d0, 1.724d0/) + 2.50d0

  !----------------------------------------------------------------------
  ! Open Files
  !
  open(11,file='./gradput/harmonics.csv')     ! n-th harmonic signals using phase-sensitive detection
  open(12,file='./gradput/reconstructed.csv') ! reconstructed first-derivative spectrum using MHD
  open(13,file='./gradput/nharmonics.csv')
  open(18,file='./gradput/snrlw2.csv')        ! filter passband, snr, linewdith
  open(27,file='./a3nakaoka.csv',action='read')
  open(28,file='./gradput/gamma-snr2.csv',action='read')
  open(29,file='./a1nakaoka.csv',action='read')  ! measured spectral data
  open(30,file='./a2nakaoka.csv',action='read')
  open(31,file='./3a1nakaoka.csv',action='read')
  open(32,file='./3a3nakaoka.csv',action='read')

  ! Type snr1 you want
  write(*,*) 'snr1'
  read(*,*) snr1
  read(28,*) alpha,noise

  noise = noise/snr1
  bm = bmo*alpha

  write(*,*)'gamma:',noise/8.52d0
  write(*,*)'alpha (modulation ratio):',alpha

  !----------------------------------------------------------------------
  !Get and put seed
  !
  call random_seed(size=seedsize)  ! get the seed size
  allocate(seed(seedsize))         ! allocate the memory space for containing the seed
  call random_seed(get=seed)
  call random_seed(put=seed)

  !----------------------------------------------------------------------
  !Spectroscopy
  !
  do  i = 1,n
    if(i<=2730)               t=1
    if(i>2730 .and. i<=5460)  t=2
    if(i>5460 .and. i<=8192)  t=3
    call random_number(x)
    do j = 1,p
      b         = (dble(i)*hs)/dble(n) + (bm/2.0d0) * dcos(dble(j) * 2.0d0 * pi / dble(p)) ! magnetic field
      nspec(i,j) =dreal(cerfe( ((b-bo(t)) + z*gamma(t)) / (sigma(t)*dsqrt(2.0d0))))&
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
  !Butterworth Filtering
  !
  call bwfilter(harmonic,bwcf)
  call bwfilter(nharmonic,bwcf)

  !--------------------------------------------------------------------
  !Read Measured Data
  !
  do i =1,n
    read(32,*)dummy,temp2(1,i),temp2(2,i),temp2(3,i),temp2(4,i),temp2(5,i)
  end do

  do i = 1,5
    mharmonic(i,1:n) = temp2(i,1:n)
  end do

  !-------------------------------------------------------------------------------
  ! Adjust the amplitudes of simulated harmonic signals (harmonic) and noise-free
  ! harmonic signals (nharmonic) to the measured harmonic signals (mharmonic)
  !
  temp(1:n) = temp2(1,1:n)
  amplitude = maxval(temp)
  temp(1:n) = dreal(harmonic(1,1:n))
  harmonic = harmonic * amplitude/maxval(temp)
  temp(1:n) = dreal(nharmonic(1,1:n))
  nharmonic(1:m,1:n) = nharmonic(1:m,1:n)*amplitude/maxval(temp)

  !--------------------------------------------------------------------
  ! Output SNR(n=1~5), both simulated and measured results
  !
  write(*,*)'output you want:',alpha,',',snr1
  write(*,*)'simulation output:,',alpha,',',hmnsnr(1,harmonic)
  write(*,*)'simulation'
  do i = 1,5
    write(*,*)hmnsnr(i,harmonic)
  end do
  write(*,*)'SNR for the measured data'
  do i =1,5
    write(*,*) hmnsnr(i,mharmonic)
  end do
  do i = 1,n
    write(11,*)(dreal(harmonic(j,i)),',',j = 1,m)
    write(13,*)(dreal(nharmonic(j,i)),',',j = 1,m)
  end do

  !---------------------------------------------------------------------
  ! FFT harmonic
  !
  do i = 1,m
    call fftcf(n,harmonic(i,1:n),harmonic(i,1:n))
    harmonic(i,1:n) = harmonic(i,1:n)/real(n)
  end do

  !---------------------------------------------------------------------
  !MHD Calculation
  !
  f_num   = (0.0d0,0.0d0)
  f_denom = (0.0d0,0.0d0)
  f       = (0.0d0,0.0d0)


  do i = 1,n/2
    if(i == 1)then
      bs(0) = 1
      bs(1:m+1) = 0
    else
      call bsjs(0.0d0, pi* bm*(dble(i))/hs, m+2, bs)
    end if
    do j = 1,10
      dn(j,i)        = (bm / (4.0 * dble(j))) * (z**(j-1)) * (bs(j-1) + bs(j+1))
      f_num(i)       = f_num(i)       + conjg(dn(j,i))     * harmonic(j,i)
      f_denom(i)     = f_denom(i)     + conjg(dn(j,i))     * dn(j,i)
      if(i /= 1) then
        if(mod(j,2) == 0) then
          f_num(n-i+2)   = f_num(n-i+2)   + conjg(dn(j,i))     * harmonic(j,n-i+2) * (-1)
          f_denom(n-i+2) = f_denom(n-i+2) + conjg(dn(j,i))     * dn(j,i)
        else
          f_num(n-i+2)   = f_num(n-i+2)   + conjg(dn(j,i))     * harmonic(j,n-i+2)
          f_denom(n-i+2) = f_denom(n-i+2) + conjg(dn(j,i))     * dn(j,i)
        end if
      end if
    end do
    f(1,i) = f_num(i) / f_denom(i)
    f(1,n-i+2) = f_num(n-i+2) / f_denom(n-i+2)
  end do

  !---------------------------------------------------------------------
  ! Output the SNR and the linewidth as a function of the passband
  !
  array(1,1:n) = f(1,1:n)
  do i = 1,100
    f(1,1:n) = array(1,1:n)
    call lpfilter(array,f,dcf*dble(i))

    !Inverse fft-------------------------------------------------------
    call fftcb(n,f(1,1:n),f(1,1:n))

    !Output the linewidth and the SNR-------------------------------------------------
    width = abs(maxloc(dreal(f(1,3500:4500)),1) - minloc(dreal(f(1,3500:4500)),1))*hs/dble(n)
    snr = hmnsnr(1,f)
    write(18,*)dcf*i,',',snr,',',width

    !Output Spectrum (when CutoffFrequency = cf) ----------------------
      if(i==int(cf))then
      write(12,*)'simulation',',','SNR:',',',snr,',','LW:',',',width
      do j =1,n
        write(12,*)dreal(f(1,j))
      end do
    end if

  end do

  close(11)
  close(12)
  close(13)
  close(18)
  close(22)
  close(23)
  close(24)
  close(25)
  close(26)
  close(27)
  close(28)
  close(29)
  close(30)
  close(31)
  close(32)


end program mhdcalc
