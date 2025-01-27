!----------------------------------------------------------------
!       Energy calculation 
!----------------------------------------------------------------
      

	Program chiparameter
	Implicit none

	Integer i, j, N,switch,Nn, p, num, k, Nt
	Integer frames
	Double precision rx(500), ry(500), rz(500)
	Double precision rxij, ryij, rzij
	Double precision rxi, ryi, rzi
	Double precision rij,rij2, box, box2
	Double precision ex(500), ey(500), ez(500)
	Double precision rcut, chi_fin
	Double precision angi,chi,chitot
	Double precision sinni,cosni,thetaij
	Double precision realp, imagp, dummy
	Character He1(500)*3, He2(500)*3, dum

	N = 500
!	frames = 500
!	switch = 4
	chi_fin = 0.0
	Print *, 'Enter the number of frames: '
	read *, frames
	print *, 'Enter the vale of switch: '
	read *, switch
	if (switch .eq. 4) then
	    p = 4
	    Nn = 8
	elseif (switch .eq. 6) then
	    p =6
	    Nn = 6
	endif

!--------reading vmd file-------------------------------	
	Do 15 k = 1,frames
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
	      Read(11,1100) He1(i), rx(i), ry(i), rz(i)
	      Read(11,1100) He2(i), ex(i), ey(i), ez(i)
	   Enddo	
! 	open(11, file='config6_4.txt')
!	     do j = 1, N
!	        read(11, *) rx(j), ry(j)!, rz(j)
!	        read(11, *) ex(j), ey(j)!, ez(j)
!	        print *,
!	        print *, 'rx =', rx(j), 'ry =', ry(j)!, 'rz =', rz(j)
!	        print *, 'ex =', ex(j), 'ey =', ey(j)!, 'ez =', ez(j)
!	     enddo
	

	   rcut = 2.0d0
	   box = 2.0*box2

	   chi = 0.0d0
	   chitot = 0.0


	   Do i = 1, N
!	     print *, 'i----------',i
	     num = 0
	     realp = 0.0
	     imagp = 0.0
             Do 16 j = 1, N
	        If (j .ne. i)	then
!	           print *, 'j',j
                   rxij = rx(j) - rx(i)             	
                   ryij = ry(j) - ry(i) 

                   rxij = rxij - box * anint(rxij/box)
                   ryij = ryij - box * anint(ryij/box)

!		   print *, 'rxij, ryij', rxij, ryij

                   rij2 = (rxij*rxij) + (ryij*ryij) !+ (rzij*rzij)
                   rij  = sqrt(rij2)
            
	           If (rij .lt. rcut) then         
	              If (rij .lt. 1.0d0) then
!	                  print *, rij, 'Overlapping particle'
		          go to 16
                      Elseif (rij .lt. 1.06d0) then
		          num = num + 1
!	                  print *,rij
!                         angi = (rxij*ex(i))+(ryij*ey(i))!+(rzij*ez(i))  !rxij.eai
                          angi = (rxij*1.0)+(ryij*0.0)!+(rzij*ez(i))  !rxij.eai
                          thetaij = acos(angi/rij)
		          cosni = cos(p*thetaij)
		          sinni = sin(p*thetaij)
!		          print *, 'theta =' , thetaij*180/3.14
		          realp = realp + cosni
		          imagp = imagp + sinni
                      Endif 
                   Endif
	        Endif
16	     continue


	     If (num .gt. 0) then
!	         print *, 'number of nearest neighbors = ', num
	         realp = realp/dble(num)
	         imagp = imagp/dble(num)
	         chi = sqrt((realp**2)+(imagp**2))
!	         print *, 'chi = ', chi
	         chitot = chitot + chi
	     Endif
!	     print *, 'num =' , num

	   Enddo
!	   print *, 'chitot = ', chitot
	   chitot = chitot/N
	   chi_fin = chi_fin + chitot
!	   print *, 'chitot = ', chitot
!	   print *, 'p', p
	   open(13, file ='chi.dat')
	   write(13,*) k, chitot 

           If(mod (k, 10000) .eq. 0) then
	      print *, k
	   Endif

15	Continue
 
	chi_fin = chi_fin/dble(frames)
	print *, 'Value of chiparameter is = ', chi_fin
	write(13,*) 'Value of chiparameter is = ', chi_fin


1100	Format(2X, 1A, 3(F14.7, 2X))

	Stop
	End program chiparameter




