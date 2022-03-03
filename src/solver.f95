!! @author : Antoine Schneeberger 
!! @mail : antoine.schneeberger@lam.fr

module JFNK
    ! Jacobian Free Newton-Krylov solver module  
    ! this module will replace the Minpack module 
    ! at the end of the day. 
    USE MODCTE
    USE ENV 
    USE RESIDUES

    implicit None 

    contains 



    function  Test_fcn ( N, x, n_args, args )
        ! Test function to test solvers, constructed as stupulated in MINPACK
        ! Variables:  
        ! n : integer , IN : number of unknown and equations = 2 here
        ! x : double precision (n), IN/OUT : estimate of unknowns from previous iteration
        ! Test_fcn : double precision (n), OUT : residue of equation system 
        !-----------------
        ! Test function must solve eq system of type F[X] = 0 
        
        integer :: N, n_args 
        real( kind = 8 ) :: test_fcn(n)
        double precision , dimension(n_args):: args
        real ( kind = 8 ) ::  x(n)
    
        Test_fcn(1) = x(1) - 3.0d0 *x(2) -2.0d0/args(1)
        Test_fcn(2) = 3.0d0 * x(1) - 4.0d0 * x(2)
    end function


    function  Test_heat_eq ( N, V, n_args, args )
        ! Test function to test solvers, constructed as stupulated in MINPACK
        ! Variables:  
        ! n : integer , IN : number of unknown and equations = 2 here
        ! x : double precision (n), IN/OUT : estimate of unknowns from previous iteration
        ! Test_fcn : double precision (n), OUT : residue of equation system 
        !-----------------
        ! Test function must solve eq system of type F[X] = 0 
        
        integer :: N, n_args 
        double precision, dimension(n) :: Test_heat_eq
        double precision , dimension(n_args):: args
        double precision, dimension(n) ::  V

        !IN/OUT
        double precision :: V0,Vn
        double precision, dimension(n) :: dV, ddV

        Test_heat_eq = 0.0d0
        V0 = 800.0d0
        Vn = 200.0d0

        !$OMP parallel workshare 

        dV(2:N-1) = (V(3:N) - V(1:N-2)) / (2.0d0*args(1))

        !$OMP end parallel workshare 

        dV(1) = (V(2) - V(1))/args(1)
        dV(N) = (V(N) - V(N-1))/args(1)

        !$OMP parallel workshare 

        ddV(2:N-1) = (dV(3:N) - dV(1:N-2)) / (2.0d0*args(1))

        !$OMP end parallel workshare 

        ddV(1) = (dV(2) - dV(1))/args(1)
        ddV(N) = (dV(N) - dV(N-1))/args(1)
       
        dV(1) = (V(2) - V0)/(2.0d0*args(1)) 
        dV(N) = (Vn - V(N-1))/(2.0d0*args(1))
        

        ddV(1) = (V0 + V(2) - 2.0d0*V(1)) / (args(1)**2.0d0)
        ddV(N) = (Vn + V(N-1) - 2.0d0*V(N)) / (args(1)**2.0d0)

        !$OMP parallel workshare 

        Test_heat_eq = ddV + 1.0d4 * exp(-(args(2:)-0.5d0)**2.0d0/0.02d0) + &
        &1.0d4*exp(-(args(2:)-0.1d0)**2.0d0/0.02d0)     

        !$OMP end  parallel workshare 

    end function

    function  test_exp( N, V, n_args, args )
        ! Test function to test solvers, constructed as stupulated in MINPACK
        ! Variables:  
        ! n : integer , IN : number of unknown and equations = 2 here
        ! x : double precision (n), IN/OUT : estimate of unknowns from previous iteration
        ! Test_fcn : double precision (n), OUT : residue of equation system 
        !-----------------
        ! Test function must solve eq system of type F[X] = 0 
        
        integer :: N, n_args 
        double precision, dimension(n) :: test_exp
        double precision , dimension(n_args):: args
        double precision, dimension(n) ::  V

        !IN/OUT
        double precision, dimension(n) :: dV

        test_exp = 0.0d0

        dV(2:N-1) = (V(3:N) - V(1:N-2)) / (2.0d0*args(1))
        dV(1) = (V(2) - V(1))/args(1)
        dV(N) = (V(N) - V(N-1))/args(1)

        dV(1) = (1 - V(2))/(2.0d0*args(1)) 
        dV(N) = (exp(10.0d0) - V(N-1))/(2.0d0*args(1))
        

        test_exp = dV - 10.0d0*V

    end function

    function back_substitution(N,U,b) 
        ! Function that solve the system Ux = b 
        ! where U is a upper triangular matrix (N,N) 
        ! and b is a vector N. The reference can be found 
        ! here https://algowiki-project.org/en/Backward_substitution 
        !
        !--------
        !Inputs :
        !--------
        !
        ! N : integer : size of the matrix U and vector b
        !
        ! U(N,N) : double precision : upper triangular matrix 
        !
        ! b(N) : double precision : vector of value 
        !
        !--------
        !Output :
        !--------
        !
        ! back_substitution(N) : double precision : vector of the solutions

        !IN/OUT
        integer :: N
        double precision, dimension(N,N) :: U ! upper triangular matrix 
        double precision, dimension(N) :: b ! vector of value
        double precision, dimension(N) :: back_substitution ! solution of Ux = b
        
        ! iteration var
        integer :: i,j  
        back_substitution = 0.0d0
        ! Initialisation 
        back_substitution(N) = b(N)/U(N,N)

        !$OMP parallel DO
        do i=1,N-1
            back_substitution(i) = b(i)
        end do 
        !$OMP END PARALLEL DO

        ! The algo back propagate from N to 1 
        do j=N,1,-1
            !init accumator 
            back_substitution(j) = back_substitution(j)/U(j,j)
            !$OMP PARALLEL DO 
            do i=1,j-1
                !Do the sum 
                back_substitution(i) = back_substitution(i) - U(i,j)*back_substitution(j)

            end do 
            !$OMP END PARALLEL DO
            
 
        end do 

    end function 


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
        integer,intent(in) :: N  ! vectors size
        integer,intent(in) :: N_args  ! Arguments size
        external :: func ! callable subroutine computing the function to use 
        double precision, dimension(N),intent(in) :: u ! point in the function where J.v is computed 
        double precision, dimension(N),intent(in) :: v ! Vector used in the dot product J.v
        double precision, dimension(N):: jac_vec ! J.v 
        double precision, dimension(N_args),intent(in) :: args ! Vector of arguments

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
                double precision , dimension(i_N)::  i_u
                double precision , dimension(i_N_args):: i_args 
                double precision, dimension(i_N) :: func

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

        double precision, dimension(max_iter+1) :: fu ! Evaluation of func(u) that will pass through the rotation
        double precision, dimension(N) :: fu_init ! Evaluation of func(u)

        double precision, dimension(N) :: res ! residual
        double precision :: res_norm ! norm of the residual  

        double precision, dimension(max_iter, N) :: Vk ! matrix of krylov space basis vectors
        double precision, dimension(N) :: Vk_estimate ! intermediary Krylov vector estimation 

        double precision, dimension(max_iter+1, max_iter) :: Hess ! Hessenberg matrix, used to construct Vk 

        double precision, dimension(max_iter,1) :: Sn, Cs ! Given rotation coefficients vectors 
        double precision :: QR_temp ! Temporary given rotation of a given Hessenberg matrix term
        
        double precision, dimension(:),allocatable :: lmbd 

        !iteration vars 
        integer :: i,j,k,it

        !Interface of the function used in the residual minimization
        interface 
        ! Func subroutine is defined as following :
            function func(i_N,i_u,i_N_args,i_args)
                
                ! i_N : number of equation 
                ! i_u : Point where func is evaluated 
                ! i_N_args : number of external arguments 
                ! i_args : number of external arguments 

                integer :: i_N , i_N_args 
                double precision , dimension(i_N)::  i_u
                double precision , dimension(i_N_args):: i_args 
                double precision, dimension(i_N) :: func

            end function 
        end interface  
       
        ! Initialisation of empty arrays 
        Hess = 0.0d0
        Cs = 0.0d0
        Sn = 0.0d0
        Vk = 0.0d0
        fu = 0.0d0

        ! Evaluation of func(u) 
        fu_init = func(N,U,N_args,args)

        ! Residual with initial guess du0 
        res = jac_vec(N,func,U,du0,N_args,args) - fu_init
        res_norm = norm2(res)

        ! First vector of the krylov space is the normalized res 
        Vk(1,:) = res/res_norm
        ! First component of fu that goes through the given rotation
        fu(1) = res_norm

        
        !Creation of all krylov vector 
        do it=1,max_iter+1
            k=it 
            ! Estimation of Krylov vector k 
            Vk_estimate = jac_vec(N,func,U,Vk(k,:),N_args,args)

            do j=1,k+1
                !$OMP PARALLEL workshare  
                !Orthogonalisation of the estimated vector
                Hess(j,k) = dot_product(Vk(j,:),Vk_estimate)
                Vk_estimate = Vk_estimate - Hess(j,k) * Vk(j,:) 
                !$OMP END PARALLEL workshare
            end do 
            
            
            !The next element of the Hessenberg matrix 
            
            Hess(k+1,k) = norm2(Vk_estimate)
            
            ! If Hessenberg matrix is not singular and if we did not exceed 
            ! the maximum iteration number
            if ((Hess(k+1,k) .ne. 0.0d0) .and. (k .ne. max_iter-1) ) then 
                ! Add the basis vector k+1
                !$OMP Parallel workshare 
                Vk(k+1,:) = Vk_estimate/Hess(k+1,k)
                !$OMP end parallel workshare 
            else 
                k = k-1
                exit 
            end if 
            
            !If the residue is less than asked tolerance, exit the loop
            if (abs(Hess(k+1,k)) < tol) exit 
            
            ! We apply k Given Rotations to the Hessenberg matrix 
            do i=1,k-1
                
                ! Since we need the value of H[i,k] to compute the value
                ! of the QR of  H[i+1,k]
                QR_temp = Cs(i,1) * Hess(i,k) + Sn(i,1) * Hess(i+1,k) 
                Hess(i+1,k) = -1.0d0 * Sn(i,1) * Hess(i,k) + Cs(i,1) * Hess(i+1,k)
                Hess(i,k) = QR_temp
            end do 
            
            ! Computation of the k+1 th Given Rotation 
            QR_temp = sqrt(Hess(k,k)*Hess(k,k) + Hess(k+1,k)*Hess(k+1,k))
            Cs(k,1) = Hess(k,k)/QR_temp
            Sn(k,1) = Hess(k+1,k)/QR_temp

            ! Do the k+1 th Given Rotation transforming the Hessenberg Matrix into
            ! pure upper Triangular matrix (easy to invert)            
            Hess(k,k) = Cs(k,1) * Hess(k,k) + Sn(k,1) * Hess(k+1,k)

            ! Since the transformation give an upper triangular matrix, this term is 
            ! nullified by the k+1 th Given rotation
            Hess(k+1,k) = 0.0d0
            
            ! We apply do the k+1 th rotation to fu to stay in the same basis
            fu(k+1) =  -1.0d0 * Sn(k,1) * fu(k)
            fu(k) = Cs(k,1) * fu(k)
        end do 


        ! Allocate the space for lmbd, now that we know the its size 
        ALLOCATE(lmbd(k)) 
        write(*,*) k
        write(*,*) 'GMRES res', Hess(k,k)
        !Value du in the krylov space 
        lmbd = back_substitution(k,Hess(:k,:k),fu(:k))
        
        GMRES_given = matmul(transpose(Vk(:k,:)),lmbd)

        DEALLOCATE(lmbd)

    end function 


    function solve_JFNK(N,func,U0,N_args,args,tol,max_iter)
        ! Solve an equation system F(X) = 0 using a Jacobian Free Newton Krylov 
        ! method, as described in Knoll et al 2002. This method is a derivative of a newton 
        ! method ui+1 = ui + J(ui)^(-1) * f(ui). However in our case, we can not compute the Jacobian 
        ! matrix explicitly. The vector dui is computed by minimizing the residu ||f(ui) - J(ui)dui||
        ! with the given rotation GMRES method. 
        !------
        !INPUTS
        !------
        !
        ! N : interger : number of equations in the system 
        !
        ! func(N,U,N_args,args) -> double precision(N) : function we are minimizing 
        !
        ! U0(N) : double precision : Initial guess
        !
        ! N_args : integer : number of arguments to pass in function func
        ! 
        ! args(N_args) : argument to pass in function func 
        !
        ! tol : double precision : Wanted tolerance on the minimizing of the function  
        !
        ! max_iter : integer : Maximum number of iteration allowed to do the minimization
        !
        !--------
        !Output :
        !--------
        !
        ! solve_JFNK(N) : double precision : solution of F(X) = 0
        
        !IN/OUT
        integer :: N ! number of equations 
        integer :: N_args ! number of arguments 
        double precision, dimension(N) :: U0 ! Solution guesses 
        double precision, dimension(N_args) :: args !arguments to be passed to the function 
        double precision :: tol ! Wanted tolerance 
        integer :: max_iter ! maximum number of iteration
        external :: func  ! Function containing the system to solve 

        double precision, dimension(N) :: solve_JFNK, solve_JFNK_test ! solution, and test solution for line search 

        ! Internals 
        double precision :: res, res_test ! residual of the  function, and test residual for line search
        double precision, dimension(N) :: du0 ! initial gmres step guess 
        double precision, dimension(N) :: du ! newton step preformed
        double precision :: line_coef = 1.0d0 ! coef for line search 

        !Interface of the function used in the func minimization
        interface 
        ! Func subroutine is defined as following :
            function func(i_N,i_u,i_N_args,i_args)
                
                ! i_N : number of equation 
                ! i_u : Point where func is evaluated 
                ! i_N_args : number of external arguments 
                ! i_args : number of external arguments 

                integer :: i_N , i_N_args 
                double precision , dimension(i_N)::  i_u
                double precision , dimension(i_N_args):: i_args 
                double precision, dimension(i_N) :: func

            end function 
        end interface      

        du0 = 0.0d0 ! fill du0 with 0s
        solve_JFNK = U0 ! initialize the solution with initial guess
        res = norm2(func(N,U0,N_args,args))

        open(unit=130,file=trim(env_datapath)//'/res.dat',status='new')

        do while (res > tol )
            ! Find the newton step with grmes given 
            du = GMRES_given(N,func,solve_JFNK,du0,N_args,args,1.0d-10,max_iter)

            !update guess with newton step 
            ! Since the step might overshoot the convergence point 
            ! a line search of the good mixing ratio is done

            solve_JFNK_test = solve_JFNK + du

            res_test = norm2(func(N,solve_JFNK_test,N_args,args))
            
            do while (res_test > res + 0.1d0 *res ) 
                
                line_coef = line_coef/2.0

                solve_JFNK_test = solve_JFNK + line_coef * du 

                res_test = norm2(func(N,solve_JFNK_test,N_args,args))

            end do 

            write(*,*) 'line_coef' , line_coef 

            line_coef = 1.0d0

            res = res_test 
            solve_JFNK = solve_JFNK_test

            write(130,*) res
            
            write(*,*) 'res', res
            write(*,*) 'x', solve_JFNK

        end do 

        close(unit=130)

    end function 

    function run_test()
        ! Function that test the JFNK solver with a known equation 
        ! system and write the results in datapath
        !
        integer,PARAMETER :: N=5000 
        integer :: i
        double precision, dimension(N) :: r
        double precision :: dr 
        double precision, dimension(N) :: guess
        double precision,dimension(N) :: X
        integer :: run_test 

        forall(i = 1: N) r(i) = float(i)/(float(N+1))
        dr = 1/float(N)

    
        guess = 1.0d0

        X=solve_JFNK(N,Test_heat_eq,guess,N+1,[dr,r],3.0d-5,500)

        open(unit=20,file=trim(env_datapath)//'/heat.dat', status='new')
        do i=1,N
            write(20,*) r(i), X(i)
        end do 
        close(unit=20)
        run_test = 1

    end function 

end module 