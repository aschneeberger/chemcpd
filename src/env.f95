!! @author : Antoine Schneeberger 
!! @mail : antoine.schneeberger@lam.fr
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
        character(len=255) :: datapath   ! Generated datapath
        character(len=255) :: dir_elems   ! Generated datapath
        
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

    subroutine write_file(fname, values, N_rows, N_cols, colnames)
        ! Subroutine parsing, formating and writing datas in a file
        !-----
        !INPUT
        !-----
        !
        ! fname : character : name of the  file 
        ! values  : double precision : data to write 
        ! N_rows : integer : number of rows  
        ! N_cols : integer : number of columns
        ! colnames : character : names of the columns formatted as "col1 col2 col3" 

        !IN/OUT
        character(len=*),intent(in) :: fname ! file name 
        character(len=*),intent(in) :: colnames ! columns names 

        integer,intent(in) :: N_rows, N_cols 
        double precision, dimension(:),intent(in) :: values

        !Internals
        integer :: i,j ! iterators 
        character(len=255) :: out_format 

        ! Check if Number of values is coherent with the asked numbers of rows and cols 
        if (N_rows * N_cols .ne. size(values)) then 

            WRITE(*,*) "[ENV] ERROR In write_file subroutine, values array size not equal to N_rows * N_cols"
            WRITE(*,*) "File :" filename
            WRITE(*,*) "Header :", colnames
            WRITE(*,*) "array length :" , size(values)
            stop
            
        end if 


        ! Create the file 
        open(unit=300,file=Trim(env_datapath)//'/'//Trim(fname),status='new')

        ! Write the columns names 
        write(300,'(A)') colnames
        ! Format and write values

        101 format(*(g0, ", "))
        do i=1,N_rows

            ! write each line of the file 

                write(300,101) (values(i+(j-1)*N_rows) , j=1,N_cols)
                
            write(300,*)

        end do 

        close(unit=300)
        
    end subroutine 

end module