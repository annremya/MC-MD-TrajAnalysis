	Program Orientation paramter
	Implicit none

	Integer i, j, k, a, N, Np, Nt 
	Integer nhis, ngr, ind, nhis1,frames
	Double precision box2
	Double precision box, dummy
	Double precision rx(500), ry(500), rz(500)
	Double precision ex(500,1), ey(500,1), ez(500,1)
	Double precision fx(500,1), fy(500,1), fz(500,1)
	Double precision rxij, ryij, rzij
	Double precision rxi, ryi, rzi, fi, fj
	Double precision rij,rij2,f
	Double precision deln, pi,dist,theta
	Double precision p(1000),ndot,rcut
	Character He1(500)*3, He2(500)*3, dum


	N = 500
	Np = 1
!	frames = 500
	pi = 3.141519
	nhis = 500
	rcut = 2.0d0
	ind = 0
	Print *, 'Enter the number of frames: '
	read *, frames

	!--------reading vmd file-------------------------------	
	do 15 k = 1,frames
	   open(11, file = 'runningprod.xyz')
	   Read(11,*) Nt
	   Read(11,*)
!	   Read(11,1100) dum, dummy, dummy, dummy
	   Read(11,1100) dum, box2, dummy, dummy
!	   Read(11,1100) dum, dummy, dummy, dummy
!	   Read(11,1100) dum, dummy, dummy, dummy
!	   Read(11,1100) dum, dummy, dummy, dummy
!	   Read(11,1100) dum, dummy, dummy, dummy
!	   Read(11,1100) dum, dummy, dummy, dummy 
!	   Read(11,1100) dum, dummy, dummy, dummy
	   Do i = 1, N
	      Read(11,1100) He1(i),rx(i), ry(i), rz(i)
	      do j = 1, Np
		   Read(11, 1100) He2(i) , ex(i,j), ey(i,j), ez(i,j)
	      enddo
	    Enddo
!	    print *, 'box2', box2

	    Do i = 1, N
	       Do j = 1, Np
		    fx(i,j) = 20.0*(ex(i,j)-rx(i))	!unit patch vector for 50  patch coverages
		    fy(i,j) = 20.0*(ey(i,j)-ry(i))	!since He2 is written as rx+(ex/20.0)
	       Enddo
	    Enddo
	 !--------initialise P(n1.n2)------------------------------
	   If (k.eq.1) then
	      box = box2*2.0
	      ngr = 0
	      deln = 2.0/(2.0*dble(nhis))
	      print *,box, box2, deln
	      do i = 1, 2*nhis
	         p(i) = 0
		enddo
	   Endif
	!-----------sample P(n1.n2)---------------------------------
	   ngr = ngr+1
	   do i = 1, N-1
	     	do 16 j = i+1, N
           	   rxij = rx(j) - rx(i)             	
           	   ryij = ry(j) - ry(i) 

           	   rxij = rxij - box * anint(rxij/box)
          	   ryij = ryij - box * anint(ryij/box)
		   rij2  = (rxij**2) + (ryij**2)
		   rij  = sqrt(rij2)
	   	   If (rij .lt. rcut) then         
	     	      If (rij .lt. 1.0d0) then
	    	         print *, rij, 'Overlapping particle'
		         go to 16
            	      elseif (rij .lt. 1.05d0) then
		         do a = 1, Np
		          ndot = (fx(i,a)*fx(j,a))+(fy(i,a)*fy(j,a))
			  theta = acos(real(ndot))
			  theta = (180.0/3.14)*theta
!			  print *, 'ndot',ndot
!			  print *, 'theta',theta
		         enddo
		         ind = int((real(ndot)+1.0)/deln)+1	
!			 print *, 'ind', ind
		         p(ind) = p(ind)+2
		      Endif
		   Endif
16		continue
	   enddo

         If(mod (k, 10000) .eq. 0) then
	   print *, k
	 Endif

15 	continue
	   print*,'ngr is', ngr
	   do i = 1, 2*nhis
!	      print *,  p(i)
	   enddo

	!-----------calculate P(n1.n2)-----------------------------
	do i =1, (2*nhis)
	   dist = (deln*(dble(i)+0.5))-1.0
!	   print *, i, dist
	   p(i) = p(i)/(dble(ngr*N)) 
	   open(17,file='ndot.dat')
	   write(17, *) dist, p(i)
	enddo


1100	Format(2X, 1A, 3(F14.7, 2X))


	Stop
	End program Orientation paramter




