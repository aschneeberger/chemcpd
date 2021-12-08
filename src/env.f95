Module ENV 
    ! Module made for working environment management 
    ! It containes all path variables and will aim 
    ! to contains parallelisation status. 
    ! it will work like a classe, wich is sketchy.
    ! All public variables must begin by env_ 

    USE MODCTE

    implicit none 

    CHARACTER (len=:), allocatable ::  env_datapath 

    contains 

    subroutine init_env()
        ! Will initialize env variables. 
        ! It will initialize the data path.
        
        ! INTERNALS
        character(len=255) :: path       ! Current working path  
        character(len=:), ALLOCATABLE :: datapath   ! Generated datapath
        character(len=:), ALLOCATABLE :: dir_elems   ! Generated datapath
        
        ! Get the working directory 
        call getcwd(path) 

        !Add to the path the relative area where all data folders are stored

        ! From the working directory and the asked data directory name, 
        ! create the absolute data path 
        datapath = Trim(path)//"/"//Trim(p_datadir)

        ! Check if p_datadir already exist
        dir_elems = listdir(path)

        ! the index method return 0 if there is no folder of the name p_datadir 
        if (index(dir_elems,Trim(p_datadir)) == 0) then 

            ! Create the data directory if it does not exist 
            ! by a call to the terminal 
            call execute_command_line('mkdir '//DATAPATH)

        else 

            ! Verify data files are present in p_datadir. If there are,
            ! The programme Stop and send an Error.
            dir_elems = listdir(datapath) 

            if (index(dir_elems,".dat") /= 0) stop '[ENV ERROR] Data files already present in the data directory'

        end if 

        ! Setup the env_datapath variable in the module
        env_datapath = datapath

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
        ! listdir : character*n : string with all the elements in the path sparated by a comma

        !IN/OUT
        CHARACTER(len=*) :: path  !The size of path is inferred from the input var 
        CHARACTER(len=:) , ALLOCATABLE :: listdir      
        CHARACTER(len=255) :: file


        !INTERNALS
        integer :: stat !temp file reading status

        !The strategie is to write the ouput ofa ls command to a tmp file and then read it
        write(*,*) path
        call execute_command_line("ls "//Trim(path)//" >> "//Trim(path)//"/listdir.tmp")

        !Open the listdir.tmp file 
        open(unit=70,file=Trim(path)//"/listdir.tmp",status='old',action='read')

        !Parse the files in a character array 
        do   
            read(70,'(A)',iostat=stat) file
            !verify if we are at end of file, if so the loop stop
            IF(IS_IOSTAT_END(stat)) exit

            !concatenate all files with a delimiter 
            listdir=listdir//Trim(file)//','

        end do  
        close(unit=70)

        !garbage collection, the tmp file is not needed anymore
        call execute_command_line("rm "//Trim(path)//"/listdir.tmp")
        
    end function 

end module