!! @author : Antoine Schneeberger 
!! @mail : antoine.schneeberger@lam.fr

module JFNK
    ! Jacobian Free Newton-Krylov solver module  
    ! this module will replace the Minpack module 
    ! at the end of the day. 
    USE MODCTE 

    implicit None 

    contains 


    ! function solver_jfnk(X,args,minbounds,maxbounds)
    !     ! API fucntion of the jfkn solver, is the front_end of the module
    !     ! This function will call succesively all needed funtion of the 
    !     ! of the module to performe a jacobian Free Newton-Krylov resolution

    ! end function 



    function jac_vec(N,func,u,v,N_args,args) 
        ! Compute the dot product of the Jacobian of 
        ! The equation system of function fun and vector vec at 
        ! the point x in the N dim space.
        ! The jacobian itself is not computed but the J·v 
        ! is approximated as (fun(u+eps*v) - fun(u))/(eps)
        ! eps being a small perturbation close to computation 
        ! precision. 
        !-------
        !Input :
        !-------
        ! 
        ! fun(double precision(N)) -> double precision(N) 
        !   : function for witch the J.v is evaluated 
        !
        ! u(N) : double precision : Point in function space where the J.v evaluated 
        !
        ! v(N) : double precision : Vector of the dot product J.v 
        !
        ! N : integer : number of points
        !
        ! N_args : interger : Number of arguments to pass in the function 
        !
        ! args(N_args) : double precision : arguments to pass in the function
        !
        !-------
        !Ouput :
        !-------
        !
        ! jac_vec(N) : The dot product J.v 


        !IN/OUT
        integer :: N  ! vectors size
        integer :: N_args  ! Arguments size
        external :: func ! callable subroutine computing the function to use 
        double precision, dimension(N) :: u ! point in the function where J.v is computed 
        double precision, dimension(N) :: v ! Vector used in the dot product J.v
        double precision, dimension(N) :: jac_vec ! J.v 
        double precision, dimension(N_args) :: args ! Vector of arguments

        !INTERNALS
        double precision :: eps ! Epsilon close to machine precision to compute derivative
        double precision :: norm_u ! L2 norm of U
        double precision :: norm_v ! L2 norm of V
        double precision, dimension(N) :: fvec, fvec_eps ! fvec = func(u) and fvec_eps = func(u + v*eps)

        ! Interface of the parametric subroutine func 
        interface 
        ! Func subroutine is defined as following :
            function func(i_N,i_u,i_N_args,i_args)
                
                ! i_N : number of equation 
                ! i_u : Point where func is evaluated 
                ! i_N_args : number of external arguments 
                ! i_args : number of external arguments 

                integer :: i_N , i_N_args 
                double precision , dimension(i_N)::  func , i_u
                double precision , dimension(i_N_args):: i_args 

            end function 
        end interface 

        norm_u = norm2(u)
        norm_v = norm2(v)

        if (norm_v == 0.0d0) then 
            ! if the perturbator vector have a null norm then the dot product is a null vector 
            jac_vec = v * 0.0d0
        else 
            ! Compute the eps used in the derivative, optimized to minimize machine error
            ! Taken from Knoll et al 2002
            eps = sqrt( (1.0d0+ norm_u) * p_machine_precision ) / norm_v
            
            !Compute J.v 
            fvec = func(N,u,N_args,args)
            fvec_eps = func(N,u+v*eps,N_args,args)
            jac_vec = (fvec_eps - fvec) / eps 
        end if 
        
        
    end function 

    function GMRES_given(N,func,U,dU0,N_args,args,tol,max_iter)
        ! Given rotation GMRES version, it include the computation of the Arnoldi basis and 
        ! the reduction of the Vk vector basis matrix. It return the du minimizing the 
        ! residual (||func(u) - Jdu||₂) for a system of dimension N. Where func(u) is
        ! the system we want to solve by JFNK, J is the jacobian of the system and du is  
        ! the newton step we are looking for. For more information on the method :
        ! Knoll et al 2009 and  Saad et al 1986.
        !
        !------
        !INPUTS
        !------
        !
        ! N : interger : number of equations in the system 
        !
        ! func(N,U,N_args,args) -> double precision(N) : function we are minimizing 
        !
        ! U(N) : double precision : Point where dU is searched 
        !
        ! dU0(N) : double precision : Initial guess of the newton step 
        !
        ! N_args : integer : number of arguments to pass in function func
        ! 
        ! args(N_args) : argument to pass in function func 
        !
        ! tol : double precision : Wanted tolerance on the minimizing of the residual 
        !
        ! max_iter : integer : Maximum number of iteration allowed to do the minimization
        !
        !-----
        !OUPUT
        !-----
        !
        ! GMRES_given(N) : double precision : newton step minimizing ||func(u) - Jdu||₂

        !IN/OUT 

        integer :: N  ! number of equations 
        integer :: N_args ! number of arguments 

        double precision, dimension(N) :: U ! point where the residual is minimized 
        double precision , dimension(N) :: dU0 ! Init guess of newton step 
        double precision , dimension(N_args) :: args ! args to pass in  func 

        double precision :: tol ! Wanted tolerance on the minimizing of the residual 
        integer  :: max_iter ! maximum number of iteration allowed to do the minimization
        
        double precision, dimension(N) :: GMRES_given 
        external :: func 

        ! Internal vars  

        double precision, dimension(N) :: fu ! Evaluation of func(u) that will pass through the rotation
        double precision, dimension(N) :: fu_init ! Evaluation of func(u)

        double precision, dimension(N) :: res ! residual
        double precision, dimension(N) :: res_norm ! norm of the residual  

        double precision, dimension(max_iter, N) :: Vk ! matrix of krylov space basis vectors
        double precision, dimension(N) :: Vk_estimate ! intermediary Krylov vector estimation 
        double precision :: Vk_estimate_norm ! norm of intermediary Krylov vector estimate

        double precision, dimension(max_iter+1, max_iter) :: Hess ! Hessenberg matrix, used to construct Vk 

        double precision, dimension(max_iter,1) :: Sn, Cn ! Given rotation coefficients vectors 
        double precision :: QR_temp ! Temporary given rotation of a given Hessenberg matrix term


        !Interface of the function used in the residual minimization
        interface 
            ! Func subroutine is defined as following :
            function func(i_N,i_u,i_N_args,i_args)
                
            ! i_N : number of equation 
            ! i_u : Point where func is evaluated 
            ! i_N_args : number of external arguments 
            ! i_args : number of external arguments 

            integer :: i_N , i_N_args 
            double precision , dimension(i_N)::  func , i_u
            double precision , dimension(i_N_args):: i_args 

        end function 
        end interface 

        ! Evaluation of func(u) 
        fu_init = func(N,U,N_args,args)
        fu = fu_init  ! copy of func(u) in the array that will pass through the rotation

        res = jac_vec(N,func,U,du0,N_args,args) - fu_init

    end function 
end module 