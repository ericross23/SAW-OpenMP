program project

        implicit none
        
        integer, dimension(:), allocatable :: saw_length
        integer                            :: i, j, k, l, m, nmc, n_saw, counter
        real, dimension(:,:), allocatable  :: pos
        real                               :: rand, start, finish, avg_chains, distance, r_squared, mseted

        
        call init_random_seed()

        n_saw = 96
        nmc = 1000

        allocate(saw_length(n_saw))
        allocate(pos(n_saw+3,3))

        saw_length(1) = 4
        do i = 2, n_saw
                saw_length(i) = saw_length(i-1) + 1        
        end do

        write(*,'(A11,A14,A20,A10)') "SAW length ", "Avg # SAW's ", "Mean square displacement ", "time"

        do i = 1, n_saw
                avg_chains = 0
                mseted = 0
                call cpu_time(start)
                do j = 1, nmc
                        k = 2
                        pos(:,:) = 0
                        counter = 1 
                        do while (k <= saw_length(i))
                                10 call random_number(rand)
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
                                                        k = 2
                                                        counter = counter + 1 
                                                        goto 10
                                                end if
                                        end do
                                end do
                                k = k + 1
                        end do
                        avg_chains = avg_chains + counter
                        distance = sqrt(pos(saw_length(i),1) ** 2 + pos(saw_length(i),2) ** 2 + pos(saw_length(i),3) ** 2)
                        r_squared = distance ** 2
                        mseted = mseted + r_squared
                end do
                avg_chains = avg_chains / nmc
                mseted = mseted / nmc
                call cpu_time(finish)
                print*, saw_length(i),  avg_chains, mseted, finish - start        
        end do
 
contains
        subroutine init_random_seed()
                implicit none
                integer                          :: i, n, clock
                integer,allocatable,dimension(:) :: seed

                call random_seed(size = n)
                allocate(seed(n))

                call system_clock(count=clock)

                seed = clock + 37 * (/ (i-1, i = 1, n) /)
                call random_seed(put = seed)

                deallocate(seed)
        end subroutine
end program
