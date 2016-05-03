program project
        use omp_lib
        implicit none
        
        character(len=5)                   :: saw_length_arg, num_threads_arg
        integer, parameter                 :: nmc = 2
        integer                            :: i, j, k, l, m, saw_length, NUM_THREADS, seed
        real, dimension(:,:), allocatable  :: pos
        double precision                   :: rand, start, finish, distance, r_squared, mseted, time                    
        common /RAND/ SEED
        
        !$omp threadprivate(/RAND/)
         
        start = omp_get_wtime()

        ! Command line arguments that specify the length of the SAW
        ! and the number of OpenMP threads to use
        call get_command_argument(1, saw_length_arg)
        call get_command_argument(2, num_threads_arg)
        read (saw_length_arg, *) saw_length
        read (num_threads_arg, *) NUM_THREADS
        allocate(pos(saw_length,3))

        open(20,file="position.stdout")

        print*, "Number of threads: ", NUM_THREADS

        ! Initial seed value for PRNG
        seed = 1

        mseted = 0

        !$omp parallel num_threads(NUM_THREADS)
        seed = omp_get_thread_num() + 1
        !$omp do schedule(dynamic), private(distance,r_squared,pos,rand), reduction(+:mseted)  
        do j = 1, nmc
                pos(1,:) = 0
                call dimer(saw_length,pos)

                ! Sample SAW position for visualization
                if (j == nmc / 2) then
                    do i = 1, saw_length
                        write(20,*) pos(i,:)
                    end do
                end if

                ! Calculate the end to end distance, square it, and add to the
                ! running average
                distance = sqrt(pos(saw_length,1) ** 2 + pos(saw_length,2) ** 2 + pos(saw_length,3) ** 2)
                r_squared = distance ** 2
                mseted = mseted + r_squared
        end do
        !$omp end do
        !$omp end parallel

        ! Finish calculating the mean squared end to end distance
        mseted = mseted / nmc

        ! Calculate the total time
        finish = omp_get_wtime()
        time = finish - start

        print*, "SAW Length: ", saw_length, "Mean squared displacement: ", mseted, "Total time: ", time

        contains

!--------------------------------------------------------------------------------------------------------------------
!
!
!                                         Parallel Random Number Generator
!
!
!--------------------------------------------------------------------------------------------------------------------

        double precision function get_randnum() result(randnum)
            integer :: seed
            common /RAND/ SEED
            
            seed = seed * 65539
            if (seed < 0) seed = (seed + 1) + 2147483647
            randnum = seed * 0.4656613e-9
        end function

!----------------------------------------------------------------------------------------------------------------------
!
!
!                                       Non-Reversal SAW Generator
!
!
!----------------------------------------------------------------------------------------------------------------------

        subroutine nrssamp(saw_length, pos)
            implicit none
            
            real,dimension(:,:), allocatable :: pos
            integer                          :: saw_length, k, l, m
            double precision                 :: rand

            10 do k = 2, saw_length
                    ! Add a point in a random direction to the existing walk
                    25 rand = get_randnum()
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

                    ! Eliminate immediate reversals
                    if (k > 2 .and. pos(k,1) == pos(k-2,1) .and. pos(k,2) == pos(k-2,2) &
                    .and. pos(k,3) == pos(k-2,3)) go to 25 

                    ! Check to see if the new point overlaps with an
                    ! existing point
                    do l = 1, k
                            do m = l + 1, k
                                    if (abs(pos(l,1) - pos(m,1)) < 1 .and. abs(pos(l,2) - pos(m,2)) < 1 &
                                                 .and. abs(pos(l,3) - pos(m,3)) < 1)  then
                                            goto 10
                                    end if
                            end do
                    end do
            end do

            
        end subroutine

!------------------------------------------------------------------------------------------------------------------------
!
!
!                                            Dimerization SAW Generator
!
!
!------------------------------------------------------------------------------------------------------------------------

        subroutine dimer(saw_length, pos)
            implicit none

            real, dimension(:,:), allocatable :: pos
            real, dimension(:,:), allocatable :: N1_pos, N2_pos
            integer                           :: saw_length, k, l, m, N_0, N_1, N_2
            double precision                  :: rand

            allocate(N1_pos(saw_length,3))
            allocate(N2_pos(saw_length,3))

            N_0 = 30
            N_1 = huge(0)
            N_2 = huge(0)

            if (saw_length <= N_0) then
                call nrssamp(saw_length, pos)
            else
                N_1 = saw_length / 2
                N_2 = saw_length - N_1
                
                N1_pos(:,:) = 0
                N2_pos(:,:) = 0
                40 call dimer(N_1, N1_pos)
                call dimer(N_2, N2_pos)

                ! Shift N_2 so that the first point of N_2 lines up with last point of N_1
                rand = get_randnum()
                do k = 1, N_2
                    if (0 <= rand .and. rand < (1./6)) then
                        N2_pos(k,:) = N2_pos(k,:) + N1_pos(N_1,:) + (/ 1, 0, 0 /) 
                    else if ((1./6) <= rand .and. rand < (2./6)) then
                        N2_pos(k,:) = N2_pos(k,:) + N1_pos(N_1,:) + (/ -1, 0, 0 /) 
                    else if ((2./6) <= rand .and. rand < (3./6)) then
                        N2_pos(k,:) = N2_pos(k,:) + N1_pos(N_1,:) + (/ 0, 1, 0 /) 
                    else if ((3./6) <= rand .and. rand < (4./6)) then
                        N2_pos(k,:) = N2_pos(k,:) + N1_pos(N_1,:) + (/ 0, -1, 0 /) 
                    else if ((4./6) <= rand .and. rand < (5./6)) then
                        N2_pos(k,:) = N2_pos(k,:) + N1_pos(N_1,:) + (/ 0, 0, 1 /) 
                    else
                        N2_pos(k,:) = N2_pos(k,:) + N1_pos(N_1,:) + (/ 0, 0, -1 /) 
                    end if
                end do

                ! Concatenate
                do k = 1, N_1
                    pos(k,:) = N1_pos(k,:)
                    !print*, pos(k,:), k, "N_1"
                end do

                do k = N_1 + 1, N_2 + N_1
                    pos(k,:) = N2_pos(k-N_1,:)
                    !print*, pos(k,:), k, "N_2"
                end do
            
                ! Check to see if the new point overlaps with an
                ! existing point
                do l = 1, saw_length
                        do m = l + 1, saw_length
                                if (abs(pos(l,1) - pos(m,1)) < 1 .and. abs(pos(l,2) - pos(m,2)) < 1 &
                                    .and. abs(pos(l,3) - pos(m,3)) < 1)  then
                                        N1_pos(:,:) = 0
                                        N2_pos(:,:) = 0
                                        !print*, "HELLO"
                                        goto 40
                                end if
                        end do
                end do
            end if
        end subroutine
 
end program
