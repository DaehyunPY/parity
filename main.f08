module input
    implicit none
    
    double precision, parameter ::      pi = 2.d0*acos(0.d0)
    complex(kind(0.d0)), parameter ::   im = cmplx(0.d0, 1.d0, kind(0.d0))
end module input





module md_h1
    use input
	implicit none 

    double precision, parameter ::     x_a = 10.d0 !-a~a
	integer, parameter ::              x_n = 2**8 !-a~a
    double precision, parameter ::     x_d = (2.d0*x_a)/dble(2*x_n)
    double precision, save ::          x(-x_n:x_n)

    double precision, parameter ::     k_a = 2.d0*pi/(2.d0*x_a)*(dble(2*x_n)/2.d0) !-a~a
    integer, parameter ::              k_n = x_n !-a~a
    double precision, parameter ::     k_d = 2.d0*pi/(2.d0*x_a)
    double precision, save ::          k(-k_n:k_n)

    complex(kind(0.d0)), save ::       psi_1x(-x_n:x_n), psi_2x(-x_n:x_n), psi_1k(-k_n:k_n), psi_2k(-k_n:k_n) 
    double precision, parameter ::     mu1 = 2.d0, sigma1 = 0.1d0, mu2 = -2.d0, sigma2 = 0.1d0 
    double precision, parameter :: 	   phase1 = 0.d0, pimentum1 = -30.d0, phase2 = 0.d0, pimentum2 = 30.d0
    complex(kind(0.d0)), save ::       psi_tx(-x_n:x_n, -x_n:x_n), psi_tk(-k_n:k_n, -k_n:k_n) 
contains





function gauss(x, mu, sigma) result(f)
    double precision, intent(in) :: x, mu, sigma
    double precision :: f

        f = 1.d0/(sigma*(2.d0*pi)**0.5d0) &
                *exp(-(x -mu)**2.d0/(2.d0*sigma**2.d0))        

end function gauss





subroutine hamilton
	double precision :: norm, parity
    integer :: i, j, l
    
    x(:) = 0.d0
    do i = -x_n, x_n
        x(i) = x_d*dble(i)
    enddo

    k(:) = 0.d0
    do i = -k_n, k_n
        k(i) = k_d*dble(i)
    enddo

    psi_1x(:) = 0.d0
    psi_2x(:) = 0.d0
    do i = -x_n, x_n
        psi_1x(i) = gauss(x(i), mu1, sigma1)*x_d**0.5d0 &
	                    *exp(im*pimentum1*x(i) +im*phase1)
        psi_2x(i) = gauss(x(i), mu2, sigma2)*x_d**0.5d0 &
	                    *exp(im*pimentum2*x(i) +im*phase2)
    enddo
    psi_1k(:) = 0.d0
    psi_2k(:) = 0.d0

    psi_tx(:, :) = 0.d0
    norm = 0.d0
    do i = -x_n, x_n
    	do j = -x_n, x_n
	    	psi_tx(i, j) = psi_1x(i)*psi_2x(j) -psi_1x(j)*psi_2x(i)
	    	norm = norm +dble(conjg(psi_tx(i, j))*psi_tx(i, j))
	    enddo
	enddo
	psi_tx(:, :) = psi_tx(:, :)/norm**0.5d0
	psi_tk(:, :) = 0.d0

	norm = 0.d0
	parity = 0.d0
    do i = -x_n, x_n
    	do j = -x_n, x_n
	    	norm = norm +dble(conjg(psi_tx(i, j))*psi_tx(i, j))
	    	parity = parity +dble(conjg(psi_tx(j, i))*psi_tx(i, j))
	    enddo
	enddo
	write(*, *) int(0.d0), 0.d0*100.d0, norm, parity
end subroutine hamilton
end module md_h1










module md_h2
    use input
    implicit none

    double precision, parameter ::  t_a     = 0.1d0 !jikan hatten no ookisa
    integer, parameter ::           t_n     = 2**8
    double precision, parameter ::  t_d     = t_a/dble(t_n)
    double precision, save ::       t(0:t_n)
contains






subroutine time
    integer :: i
    
    t(:) = 0.d0
    do i = 0, t_n
        t(i) = t_d*dble(i)
    enddo

end subroutine time
end module md_h2






















program main 
    use md_h1
    use md_h2
    use fftsg
    implicit none

    integer, parameter :: file_psi_1x  = 11
    integer, parameter :: file_psi_1k  = 12
    integer, parameter :: file_psi_2x  = 21
    integer, parameter :: file_psi_2k  = 22
    integer, parameter :: file_psi_tx  = 31
    integer, parameter :: file_psi_tk  = 32
    integer, parameter :: file_psi_tx_zero  = 41
    integer, parameter :: file_psi_tk_zero  = 42

    integer, parameter :: t_plot      = 20
    integer, parameter :: x_plot      = 50
    integer, parameter :: t_pp        = t_n/t_plot  
    integer, parameter :: x_pp        = x_n/x_plot  

    double precision :: norm, parity
    integer :: i, j, j1, j2
    character(3) :: ch3

    call hamilton
    call time

!     open(file_psi_1x, file = 'output/psi_1x.d')
!     open(file_psi_1k, file = 'output/psi_1k.d')
!     open(file_psi_2x, file = 'output/psi_2x.d')
!     open(file_psi_2k, file = 'output/psi_2k.d')


!     do i = 0, t_n
!     	if(i /= 0) then 
! 	        psi_1k(:) = psi_1x(:)
! 	        call fft(-1, psi_1k(-k_n +1:k_n))
! 	        do j = -k_n, k_n
! 	            psi_1k(j) = exp(-im*1.d0/2.d0*k(j)**2.d0*t_d)*psi_1k(j)
! 	        enddo
! 	        psi_1x(:) = psi_1k(:)
! 	        call fft(+1, psi_1x(-x_n +1:x_n))
! 	    endif

! 	    if(mod(i, t_pp) == 0) then
! 		    do j = -x_n +1, x_n
! 		    	if(mod(j, x_pp) == 0) then
! 				    write(11, *) t(i), x(j), dble(conjg(psi_1x(j))*psi_1x(j))
! 				endif
! 		    enddo
! 		    write(11, *)
! 		endif
!     enddo

!     close(file_psi_1x)
!     close(file_psi_1k)
!     close(file_psi_2x)
!     close(file_psi_2k)





open(file_psi_tx_zero, file = 'output/psi_x_zero.d')
open(file_psi_tk_zero, file = 'output/psi_k_zero.d')

    do i = 0, t_n
    	if(i /= 0) then 
	        psi_tk(:, :) = psi_tx(:, :)
	        call fft2d(-1, psi_tk(-k_n +1:k_n, -k_n +1:k_n))
	        do j1 = -k_n, k_n
	        	do j2 = -k_n, k_n
		            psi_tk(j1, j2) = exp(-im*1.d0/2.d0*k(j1)**2.d0*t_d) &
			            				*exp(-im*1.d0/2.d0*k(j2)**2.d0*t_d)*psi_tk(j1, j2)
		        enddo
	        enddo
	        psi_tx(:, :) = psi_tk(:, :)
	        call fft2d(+1, psi_tx(-x_n +1:x_n, -x_n +1:x_n))

	        norm = 0.d0 
	        parity = 0.d0
		    do j1 = -x_n +1, x_n
		    	do j2 = -x_n +1, x_n
			    	norm = norm +dble(conjg(psi_tx(j1, j2))*psi_tx(j1, j2))
			    	parity = parity +dble(conjg(psi_tx(j2, j1))*psi_tx(j1, j2))
			    enddo
			enddo
			psi_tx(:, :) = psi_tx(:, :)/norm**0.5d0
			write(*, *) int(i/t_pp), dble(i)/dble(t_n)*100.d0, norm, parity
	    endif

	    if(mod(i, t_pp) == 0) then
	    	open(file_psi_tx, file = 'output/psi_x_'//ch3(int(i/t_pp))//'.d')
	    	open(file_psi_tk, file = 'output/psi_k_'//ch3(int(i/t_pp))//'.d')
		    do j1 = -x_n +1, x_n
			    do j2 = -x_n +1, x_n
			    	if(mod(j1, x_pp) == 0 .and. mod(j2, x_pp) == 0) then
					    write(file_psi_tx, *) x(j1), x(j2), dble(conjg(psi_tx(j1, j2))*psi_tx(j1, j2))
					endif
					if(j1 == 0) then 
						write(file_psi_tx_zero, *) t(i), x(j2), dble(conjg(psi_tx(j1, j2))*psi_tx(j1, j2))
					endif
				enddo
			    if(mod(j1, x_pp) == 0) write(file_psi_tx, *)
			    if(j1 == 0) write(file_psi_tx_zero, *) 
		    enddo
		    close(file_psi_tx)
		    close(file_psi_tk)
		endif
    enddo

    close(file_psi_tx_zero)
    close(file_psi_tk_zero)
end program main