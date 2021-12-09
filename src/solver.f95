!! @author : Antoine Schneeberger 
!! @mail : antoine.schneeberger@lam.fr

module JFNK
    ! Jacobian Free Newton-Krylov solver module  
    ! this module will replace the Minpack module 
    ! at the end of the day. 

    implicit None 

    contains 

    public 

    function solver_jfnk(X,args,minbounds,maxbounds)
        ! API fucntion of the jfkn solver, is the front_end of the module
        ! This function will call succesively all needed funtion of the 
        ! of the module to performe a jacobian Free Newton-Krylov resolution

    end function 

    private 

    subroutine GMRES()
        ! Performe a General minimization of residues
        ! of the form || Ax - f || in the definition of L2 norme 

    end subroutine 

    function jac_vec(fun,vec,x) 
        ! Compute the dot product of the Jacobian of 
        ! The equation system of function fun and vector vec at 
        ! the point x in the N dim space.
        ! The jacobian itself is not computed but the J·v 
        ! is approximated as (fun(x+eps*vec) - fun(x))/(eps)
        ! eps being a small perturbation close to computation 
        ! precision. 

    end function 

    
end module 