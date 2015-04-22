! Module to calculate the actual rate integral, eq. 7 and 8. in the paper
! 
! [1] http://arxiv.org/abs/1502.02005
!
! the subroutines GP_l and GM_l are called by rates.py 
! for every value of the gate voltage
! Note that the module has to be initialized first with the parameters of the system
! see subroutine INIT.
!
! use
!   make geval.so
! to compile a python object file with f2py
!
! dependencies:
!   + for the integration, qag routine from the quadpack is used, see quadpack.f95
!   + for the digamma function, CPSI from SLATEC is used, a single precision routine.
!     (therefore, one has to install SLATEC!)
!
!


module geval
  implicit none
  private
  public :: init, gm_l, gp_l, selfenergy
  
  real, allocatable :: test
  integer, parameter :: dp = selected_real_kind(15,307)      ! double precision
  integer, parameter :: sp = selected_real_kind(6,37)        ! single precision

  real(dp), parameter :: pi = 3.14159265_dp
  integer                          :: LEAD
  real(dp), dimension(0:7,0:3,0:23)        :: UPSTATES,DOWNSTATES
  integer, dimension(0:255)            :: MAP
  integer, dimension(0:255,0:255)  :: TARRAY
  integer, dimension(0:39,0:1)     :: TLIST

  real(dp), dimension(0:255)   :: EARRAY
  real(dp)                     :: Vg,beta,alpha,mu
  integer                          :: pm
  real(dp), dimension(0:1,0:3)     :: W0t,W0b
  real(dp), dimension(0:1,0:3)     :: WEt,WEb
  real(dp), dimension(0:1,0:3)     :: GAMMA
  real(dp), dimension(0:1)     :: V
  integer                          :: b,a
  integer                          :: so_ab

  ! INTEGRATION PARAMETERS AND RETURN VALUES
  real(dp) :: epsabs
  real(dp) :: epsrel
  real(dp) :: int_co
  integer      :: int_key

  real(dp) :: abserr
  integer :: ier,neval

contains
  subroutine selfenergy(x,b,a,se)
    ! calculates the selfenergy components Eq. 9 in [1]
    ! 
    ! Args:
    !   x: the running variable, energy
    !   b,a: the final state b and initial state a
    !   se: the resulting self energy component, real and imaginary part
    !
    !

    real(dp),intent(in)                      :: x
    integer, intent(in)                          :: b,a
    integer                                      :: j,qn
    real(dp),dimension(0:1),intent(out)          :: se
    real(dp)                                 :: Ecd,Ecu

    se = 0
    do qn = 0,3
       do j = 0,size(downstates,1)-1
          Ecd = downstates(j,qn,map(b))
          Ecu = upstates(j,qn,map(a))
          ! Ecd > 0 for out-tunneling transitions to excited states
          if (Ecd>0) then
             se(0) = se(0) + GAMMA(0,qn)*lorentzian(x,WEt(0,qn))*p_alpha1(x,Ecd-EARRAY(a),0,-1) + GAMMA(1,qn)*lorentzian(x,WEt(1,qn))*p_alpha1(x,Ecd-EARRAY(a),1,-1)
             se(1) = se(1) + GAMMA(0,qn)*lorentzian(x,WEt(0,qn))*alpha1(x,Ecd-EARRAY(a),0,-1) + GAMMA(1,qn)*lorentzian(x,WEt(1,qn))*alpha1(x,Ecd-EARRAY(a),1,-1)
             !write (*,*) "b",b," qn ",qn ," Ediff ", Ecd-EARRAY(a), "se ", se(0), se(1)
          end if
          ! Ecu > 0 for in-tunneling transitions to excited states
          if(Ecu>0) then
             se(0) = se(0) + GAMMA(0,qn)*lorentzian(x,WEb(0,qn))*p_alpha1(x,EARRAY(b)-Ecu,0,1) + GAMMA(1,qn)*lorentzian(x,WEb(1,qn))*p_alpha1(x,EARRAY(b)-Ecu,1,1)
             se(1) = se(1) + GAMMA(0,qn)*lorentzian(x,WEb(0,qn))*alpha1(x,EARRAY(b)-Ecu,0,1) + GAMMA(1,qn)*lorentzian(x,WEb(1,qn))*alpha1(x,EARRAY(b)-Ecu,1,1)
             !write (*,*) "a",a," qn ",qn ," Ediff ", EARRAY(b)-Ecu, "se ", se(0), se(1)
          end if
       end do
    end do

    return
  end subroutine selfenergy



  function integrand(x)
    ! the function that is passed to the integration routine 
    ! it has to be a single valued function,
    ! see quadpack.f90
    
    ! Args:
    !   x: running energy value 
    ! Returns:
    !   value of the integrand for this energy value
    !
    !
    implicit none
    ! the routine parameters
    real(dp),intent(in)                      :: x
    real(dp)                                 :: integrand
    real(dp)                                 :: Z1
    real(dp),dimension(0:1)                  :: se


    call selfenergy(x,b,a,se)
    
    ! The real part of the integrand is calculated
    Z1 = GAMMA(LEAD,so_ab)*f_p(x,mu-V(LEAD),pm*beta)

    se(0) = x - get_dE(b,a) + se(0) 

    integrand = Z1*se(1)/(se(1)**2 + se(0)**2)

    return
  end function integrand


  subroutine INIT(ext_LEAD,ext_Vb,ext_Vg,ext_W0t,ext_W0b,ext_WEt,ext_WEb,ext_beta,ext_alpha,ext_mu,ext_TLIST,ext_UPSTATES,ext_DOWNSTATES,ext_MAP,ext_EARRAY,ext_GAMMA,ext_Int_CO,ext_Int_acc_abs,ext_Int_acc_rel,ext_Int_gkp_key)
    ! the parameters of the system are passed to INIT once for every rate.
    ! They are stored in the modules global variables for all subsequent
    ! calls to GP_l and GM_l.
    ! 
    ! Args:
    !   ext_LEAD: lead index for the rate calculation (0/1 = left/right)
    !   ext_Vb: bias voltage
    !   ext_Vg: gate voltage. This value will be overwritten by calls to GM_l or GP_l
    !   ext_W0t: the bandwidth for the rate integral (upper bound)
    !   ext_W0b: the bandwidth for the rate integral (lower bound)
    !   ext_W0t: the bandwidth for the excited state integrals (upper bound)
    !   ext_W0b: the bandwidth for the excited state integrals (lower bound)
    !   ext_beta: the inverse temperature 1/k_BT
    !   ext_alpha: the gate conversion factor
    !   ext_mu: the chemical potential \mu_0
    !   ext_TLIST: the transtion list for all transtions between ground states
    !   ext_UPSTATES: the energies for all allowed excited states that can be
    !                 reached by an in-tunneling event.
    !                 Refer to cfg.py for details.
    !                 Note that the size is fixed due to the non-dynamical
    !                 nature of the interface python-fortran.
    !                 If you want to change the excited
    !                 state bandwidth you have to change NX in cnt.conf AND the size
    !                 of ext_UPSTATES and its local counterpart upstates!
    !   ext_DOWNSTATES: the energies for all allowed excited states that can be
    !                 reached by an in-tunneling event
    !   ext_MAP: the statemap to identify the entries in 
    !            UP/DOWNSTATES with the STATENUMBER
    !   ext_EARRAY: the energy array for all STATENUMBERs
    !   ext_GAMMA: the couplings
    !   ext_Int_CO: the part of the rate that is integrated 
    !               although its practically zero:
    !               the extend of the side of the fermi function that is calculated
    !               although it is exponentially surpressed
    !   ext_Int_acc_abs: absoulute value of the minimum precision 
    !                    for the numerical integration
    !   ext_Int_acc_rel: relative value of the minimum precision, see quadpack.f95
    !   ext_Int_gpk_key: another tuning variable for the numerical integration
    !
    integer, dimension(0:39,0:1), intent(in)    :: ext_TLIST
    real(dp), dimension(0:255),intent(in)   :: ext_EARRAY
    real(dp),intent(in)                     :: ext_beta,ext_alpha,ext_mu
    integer, intent(in)                         :: ext_LEAD

    integer, intent(in)                         :: ext_Int_gkp_key
    real(dp),intent(in)                     :: ext_Int_CO,ext_Int_acc_rel,ext_Int_acc_abs

    real(dp),intent(in)                     :: ext_Vb,ext_Vg
    real(dp), dimension(0:1,0:3),intent(in)     :: ext_W0t,ext_W0b
    real(dp), dimension(0:1,0:3),intent(in)     :: ext_WEt,ext_WEb

    real(dp), dimension(0:1,0:3),intent(in)     :: ext_GAMMA

    real(dp), dimension(0:7,0:3,0:23),intent(in)     :: ext_UPSTATES,ext_DOWNSTATES
    integer, dimension(0:255),intent(in)           :: ext_MAP

    if(allocated(test)) print *, 'allocated'
    ! INITIALIZE EXTERNAL VARIABLES
    LEAD = ext_LEAD

    TLIST = ext_TLIST

    EARRAY = ext_EARRAY

    beta = ext_beta
    alpha = ext_alpha
    mu = ext_mu

    GAMMA = ext_GAMMA

    LEAD = ext_LEAD

    V = ext_Vb*(/1,-1/)
    Vg = ext_Vg

    W0t = ext_W0t
    W0b = ext_W0b
    WEt = ext_WEt
    WEb = ext_WEb

    epsabs = ext_Int_acc_abs
    epsrel = ext_Int_acc_rel

    int_co = ext_Int_CO
    int_key = ext_Int_gkp_key

    UPSTATES = ext_UPSTATES
    DOWNSTATES = ext_DOWNSTATES
    MAP = ext_MAP


    return
  end subroutine INIT


  pure function f_p(my_e,my_mu,my_beta)
    ! the fermi function
    real(dp)                  :: f_p
    real(dp), intent(in)      :: my_e,my_mu,my_beta
    f_p = 1./(1.+EXP(my_beta*(my_e-my_mu)))
    return
  end function f_p

  function get_dE(n2,n1)
    ! the energy difference between two states given by their STATENUMBER
    ! with the gate voltage part (stored in global Vg)
    integer,intent(in) :: n2,n1
    real(dp) :: get_dE

    get_dE = (EARRAY(n2)-EARRAY(n1)) - alpha*Vg
    !write (*,*) "Charge diff. ",CHARGE(n2)-CHARGE(n1),"Energy difference ",get_dE
    return
  end function get_dE


  subroutine GP_l(ext_Vg,GPl)
    ! Interface routing for the rate calculation of in-tunneling rates.
    !
    !    Args:
    !      ext_Vg: the gate voltage for the rate calculation
    !      GPl: the result, a float array of size 39, the number of transitions
    !
    !
    real(dp),intent(out),dimension(0:39)    :: GPl
    real(dp),intent(in)                     :: ext_Vg

    ! LOCAL VARIABLES
    real(dp)                     :: int_val
    integer                          :: i

    Vg = ext_Vg

    GPl = 0

    pm = 1

    ! we loop through all transitions in TLIST and calculate the corresponding
    ! rate. Note that the order of TLIST therefore determines the order of the
    ! returned GP array.
    do i = 0,39
       b = TLIST(i,0)
       a = TLIST(i,1)
       so_ab = spinorbit_map(a,b)
       ! write (*,*) "transition ",i," so_num ", so_ab

       ! We leave the loop as soon as there is a zero element.
       if( b == 0 )  exit

       ! 
       call qag(integrand,-W0b(LEAD,so_ab),int_co,epsabs,epsrel,int_key,int_val,abserr,neval,ier)

       GPl(i) = 2*int_val*pi
       !GPl(i) = 2*pi*integrand(0._dp)
    end do
    return
  end subroutine GP_l

  subroutine GM_l(ext_Vg,GMl)
    ! Interface routing for the rate calculation of out-tunneling rates.
    !
    !    Args:
    !      ext_Vg: the gate voltage for the rate calculation
    !      GMl: the result, a float array of size 39
    !
    !
    real(dp), intent(out), dimension(0:39)    :: GMl
    real(dp), intent(in)                      :: ext_Vg

    ! LOCAL VARIABLES
    real(dp)                     :: int_val
    integer                          :: i

    Vg = ext_Vg

    GMl = 0

    ! IN- OR OUT-TUNNELING EVENT

    pm = -1

    do i = 0,39
       ! The order of the states a and b is swapped for the out-tunneling process
       !          
       b = TLIST(i,0)
       a = TLIST(i,1)
       so_ab = spinorbit_map(a,b)

       if( b == 0 )  exit

       call qag(integrand,-int_co,W0t(LEAD,so_ab),epsabs,epsrel,int_key,int_val,abserr,neval,ier)
       GMl(i) = 2*int_val*pi 
       !GMl(i) = 2*pi*integrand(0._dp)
    end do
    return
  end subroutine GM_l

  function spinorbit_map(mya,myb)
    ! Returns the index of the four vector of the electron that is transfered
    ! when going from mya to myb.
    ! Example:
    !   [0,1,0,0] -> [0,1,1,0]
    !   result is 2.
    !
    integer,intent(in) :: mya,myb
    integer :: spinorbit_map
    spinorbit_map = 3-log(abs(real(mya-myb)))/log(4.)
    return
  end function spinorbit_map



  function alpha1(e,dE,l,my_pm)
    ! Returns the imaginary part of the inner rate integrals
    !
    ! Args:
    !   e: the running variable
    !   dE: energy difference for the transitions
    !   l: the lead index
    !   my_pm: index whether this is an in (0) or out(1)-tunneling transiton.
    !
    !
    real(dp)             :: alpha1
    real(dp), intent(in)                 :: e,dE
    integer, intent(in)                  :: my_pm,l
    real(dp)                             :: arg

    arg = (e-dE)

    alpha1 = pi*f_p(arg,mu-V(l),my_pm*beta)!lorentzian(arg-mu+V(l),W)

    !write (*,*) "e ",get_dE(n2,n1)," and alpha ",alpha1

    return
  end function alpha1

  function p_alpha1(e,dE,l,my_pm)
    ! Returns the real or principal part of the inner rate integrals
    !
    ! Args:
    !   e: the running variable
    !   dE: energy difference for the transitions
    !   l: the lead index
    !   my_pm: index whether this is an in (0) or out(1)-tunneling transiton.
    !
    !
    real(dp)                              :: p_alpha1
    real(dp)                              :: arg
    real(dp), intent(in)                  :: e,dE
    complex(sp)                           :: digamma,CPSI
    integer, intent(in)                       :: my_pm,l

    arg = (e-dE)-mu+V(l)
    
    ! note that the routine CPSI from the Fortran library SLATEC is used to
    ! calculate the digamma function
    ! to this end, the values are converted to single-precision values!
    ! the overall accuracy is limited by this fact.
    digamma = CPSI(COMPLEX(0.5_sp,0.5_sp*real(arg*beta/pi)))
    p_alpha1 = dble(my_pm*(REALPART(digamma)))

    return
  end function p_alpha1

  subroutine testcpsi(x,y)
    real(sp), intent(in) :: x
    real(sp), intent(out) :: y
    complex(sp) :: CPSI

    y = REALPART(CPSI(COMPLEX(0.5_sp,x)))
  end subroutine testcpsi

  function lorentzian(e,W)
    real(dp), intent(in)  :: e,W
    real(dp)              :: lorentzian

    lorentzian = W**2/(e**2+W**2)
    return
  end function lorentzian

  function lorentzian_complex(e,W)
    complex(dp), intent(in)  :: e
    real(dp), intent(in)     :: W
    real(dp)                 :: lorentzian_complex

    lorentzian_complex = W**2/((REAL(e)**2+AIMAG(e)**2)+W**2)
    return
  end function lorentzian_complex

end module geval



      
