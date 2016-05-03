program project
        use omp_lib
        implicit none
        
        character(len=3)                   :: saw_length_arg, num_threads_arg
        integer                            :: seed
        common /RAND/ SEED
        integer                            :: NUM_THREADS
        
        !$omp threadprivate(/RAND/)
        integer, parameter                 :: nmc = 10000
        integer                            :: j, k, l, m, saw_length
        real, dimension(:,:), allocatable  :: pos
        double precision                   :: rand, start, finish, distance, r_squared, mseted, time                    
        
        start = omp_get_wtime()

        ! Command line arguments that specify the length of the SAW
        ! and the number of OpenMP threads to use
        call get_command_argument(1, saw_length_arg)
        call get_command_argument(2, num_threads_arg)
        read (saw_length_arg, *) saw_length
        read (num_threads_arg, *) NUM_THREADS

        print*, NUM_THREADS
        ! Initial seed value for PRNG
        seed = 1

        allocate(pos(saw_length,3))

        mseted = 0
        !$omp parallel num_threads(NUM_THREADS)
        seed = omp_get_thread_num() + 1
        !$omp do schedule(dynamic), private(distance,r_squared,pos,rand), reduction(+:mseted)  
        do j = 1, nmc
                pos(:,:) = 0
                10 do k = 2, saw_length
                        rand = get_randnum()
                        if (0 <= rand .and. rand < (1./6)) then
                                pos(k,:) = pos(k-1,:) + (/ 1, 0, 0 /)        
                        else if ((1./6) <= rand .and. rand < (2./6)) then
                                pos(k,:) = pos(k-1,:) + (/ -1, 0, 0 /)
                        else if ((2./6) <= rand .and. rand < (3./6)) then
                                pos(k,:) = pos(k-1,:) + (/ 0, 1, 0 /)
                        else if ((3./6) <= rand .and. rand < (4./6)) then
                                pos(k,:) = pos(k-1,:) + (/ 0, -1, 0 /)
                        else if ((4./6) <= rand .and. rand < (5./6)) then
                                pos(k,:) = pos(k-1,:) + (/ 0, 0, 1 /)
                        else
                                pos(k,:) = pos(k-1,:) + (/ 0, 0, -1 /)
                        end if

                        do l = 1, k
                                do m = l + 1, k
                                        if (abs(pos(l,1) - pos(m,1)) < 1 .and. abs(pos(l,2) - pos(m,2)) < 1 &
                                                    .and. abs(pos(l,3) - pos(m,3)) < 1)  then
                                                goto 10
                                        end if
                                end do
                        end do
                end do
                distance = sqrt(pos(saw_length,1) ** 2 + pos(saw_length,2) ** 2 + pos(saw_length,3) ** 2)
                r_squared = distance ** 2
                mseted = mseted + r_squared
        end do
        !$omp end do
        !$omp end parallel
        mseted = mseted / nmc

        finish = omp_get_wtime()

        time = finish - start

        print*, "SAW Length: ", saw_length, "Mean squared displacement: ", mseted, "Total time: ", time
        
        contains

        double precision function get_randnum() result(randnum)
            integer :: seed
            common /RAND/ SEED
            
            seed = seed * 65539
            if (seed < 0) seed = (seed + 1) + 2147483647
            randnum = seed * 0.4656613e-9
        end function
 
end program
