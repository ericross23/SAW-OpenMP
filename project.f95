program project
        use omp_lib
        implicit none
        
        integer, dimension(:), allocatable :: saw_length
        integer                            :: seed
        common /RAND/ SEED
        integer, parameter                 :: NUM_THREADS = 4
        integer                            :: max_threads
        
        !$omp threadprivate(/RAND/)
        integer, parameter                 :: nmc = 1000
        double precision                   :: randnum
        integer                            :: i, j, k, l, m, n_saw, counter
        real, dimension(:,:), allocatable  :: pos
        real                               :: rand, start, finish, avg_chains, distance, r_squared, mseted                    
        

        n_saw = 96
        seed = 1
        max_threads = omp_get_max_threads()

        print*, max_threads


        allocate(saw_length(n_saw))
        allocate(pos(n_saw+3,3))

        saw_length(1) = 4
        do i = 2, n_saw
                saw_length(i) = saw_length(i-1) + 1        
        end do

        write(*,'(A11,A14,A20,A10)') "SAW length ", "Avg # SAW's ", "Mean square displacement ", "time"
        !$omp parallel num_threads(NUM_THREADS)
        seed = omp_get_thread_num() + 1
        !$omp do private(pos), schedule(dynamic), reduction(+:avg_chains,mseted), collapse(1) 
        do i = 1, n_saw
                avg_chains = 0
                mseted = 0
                start = omp_get_wtime()
                do j = 1, nmc
                        pos(:,:) = 0
                        counter = 1
                        10 do k = 2, saw_length(i) + 1
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
                                                        counter = counter + 1 
                                                        goto 10
                                                end if
                                        end do
                                end do
                        end do
                        avg_chains = avg_chains + counter
                        distance = sqrt(pos(saw_length(i),1) ** 2 + pos(saw_length(i),2) ** 2 + pos(saw_length(i),3) ** 2)
                        r_squared = distance ** 2
                        mseted = mseted + r_squared
                end do
                avg_chains = avg_chains / nmc
                mseted = mseted / nmc
                finish = omp_get_wtime()
                print*, saw_length(i),  avg_chains, mseted, finish - start        
        end do
        !$end omp do
        !$omp end parallel

        contains

        double precision function get_randnum() result(randnum)
            integer :: seed
            common /RAND/ SEED
            
            seed = seed * 65539
            if (seed < 0) seed = (seed + 1) + 2147483647
            randnum = seed * 0.4656613e-9
        end function
 
end program
