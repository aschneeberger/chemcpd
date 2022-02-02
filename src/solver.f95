!! @author : Antoine Schneeberger 
!! @mail : antoine.schneeberger@lam.fr

module JFNK
    ! Jacobian Free Newton-Krylov solver module  
    ! this module will replace the Minpack module 
    ! at the end of the day. 

    USE PHYCTE 

    implicit None 

    contains 


    ! function solver_jfnk(X,args,minbounds,maxbounds)
    !     ! API fucntion of the jfkn solver, is the front_end of the module
    !     ! This function will call succesively all needed funtion of the 
    !     ! of the module to performe a jacobian Free Newton-Krylov resolution

    ! end function 


    subroutine GMRES()
        ! Performe a General minimization of residues
        ! of the form || Ax - f || in the definition of L2 norme 

    end subroutine 

    function jac_vec(N,func,u,v) 
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
        !-------
        !Ouput :
        !-------
        !
        ! jac_vec(N) : The dot product J.v 


        !IN/OUT
        integer :: N  ! vectors sizes 
        external :: func ! callable subroutine computing the function to use 
        double precision, dimension(N) :: u ! point in the function where J.v is computed 
        double precision, dimension(N) :: v ! Vector used in the dot product J.v
        double precision, dimension(N) :: jac_vec ! J.v 

        !INTERNALS
        double precision :: eps ! Epsilon close to machine precision to compute derivative
        double precision :: norm_u ! L2 norm of U
        double precision :: norm_v ! L2 norm of V

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
            !jac_vec = (func(u + eps*v) - func(u)) / eps 
        end if 
        
        
    end function 

    
end module 