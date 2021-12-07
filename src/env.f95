Module ENV 
    ! Module made for working environment management 
    ! It containes all path variables and will aim 
    ! to contains parallelisation status. 
    ! it will work like a classe, wich is sketchy.
    ! All public variables must begin by env_ 

    USE MODCTE

    implicit none 

    CHARACTER (len=:), allocatable :: env_path, env_datapath 

    contains 

    subroutine init_env()
        ! Will initialize env variables. 
        ! It will initialize the data path.
        
        ! INTERNALS
        character(len=:), ALLOCATABLE :: path       ! Current working path  
        character(len=:), ALLOCATABLE :: datapath   ! Generated datapath
        
        ! Get the working directory 
        call getcwd(path) 

        ! From the working directory and the asked data directory name, 
        ! create the absolute data path 
        datapath = path//"/"//Trim(p_datadir)
        
        ! Create the data directory by a call to the terminal 
        call system('mkdir '//DATAPATH)

    end subroutine 

    subroutine get_datapath()
        !Get the datapath 
        
        

    end subroutine 

    function listdir(path) 
        !Function that return a string containing all directories and files
        !-----
        !INPUT
        !-----
        !
        ! PATH : character*n : path to the target location (absloute)
        !
        !------
        !OUTPUT
        !------
        !
        ! listdir : character*n : string with all the elements in the path

        !IN/OUT
        CHARACTER(len=:), ALLOCATABLE :: path        
        CHARACTER(len=:), ALLOCATABLE :: listdir 

        !The strategie is to write the ouput ofa ls command to a tmp file and then read it
        write(*,*) path
        call execute_command_line("ls "//Trim(path)//" >> listdir.tmp")

        open(unit=60,file="lisdir.tmp",status='read')

    end function 

end module