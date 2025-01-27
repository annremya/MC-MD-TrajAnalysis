	Program rdf
	Implicit none

	Integer i, j, k, a, m, N, Np, Nt
	Integer xi, yi, zi, nb, nv, frames
	Integer nhis, ngr, switch, ig, nhis1
	Double precision box3, box2, b, q, x0, y0, a1, a2
	Double precision l, fact, box(24), box1
	Double precision rx(500), ry(500), rz(500)
	Double precision ex(500,1), ey(500,1), ez(500,1)
!	Double precision fx(500,1), fy(500,1), fz(500,1)
	Double precision delg, pi,dist, dummy
	Double precision g(1000),dx,dy,d, ro, const
	Double precision rlow, rup, nideal, vb

	Character He1(500)*3, He2(500)*3, dum


	N = 500
	Np = 1
!	frames = 500
	pi = 3.141519
	nhis = 500
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
	 !--------initialise rdf------------------------------
	   If (k.eq.1) then
	      box1 = box2*2.0
	      ngr = 0
	      delg = box1/(2*nhis)
	      print *,box1, box2
	      do i = 1, nhis
	         g(i) = 0
	      enddo
	   Endif
	!-----------sample rdf---------------------------------
	   ngr = ngr+1
	   do a = 1, N-1
	     	do m = a+1, N
		   dx = rx(a) - rx(m)
		   dy = ry(a) - ry(m)		    
		   if(dx.gt.box2) dx = dx-box1
		   if(dx.lt.-box2) dx = dx+box1
		   if(dy.gt.box2) dy = dy-box1
		   if(dy.lt.-box2) dy = dy+box1
		   d = sqrt((dx**2)+(dy**2))
		   if (d .lt. box2) then
			ig = int(d/delg)+1
			g(ig) = g(ig)+2
!			print *, ig
		   endif
		enddo
	   enddo

         If(mod (k, 10000) .eq. 0) then
	   print *, k
	 Endif

15 	continue
	   print*,'ngr is', ngr

	!-----------calculate rdf-----------------------------
!	nhis1 = int(box2/delg)+1
	ro = dble(N)/box1**2
	const = pi * ro
	do i =1, nhis
	   dist = delg*(dble(i)+0.5)
	   rlow = dble(i)*delg
	   rup = rlow + delg
	   vb = ((rup)**2-(rlow)**2)
	   nideal = const*vb
	   g(i) = g(i)/(ngr*N*nideal) 
	   open(16,file='rdf.dat')
	   write(16, *) dist, g(i), nideal
	enddo

!-------------------------------------------------------------------------
!	    Do i = 1, N
!	       Do j = 1, Np
!		    fx(i,j) = 5.0*(ex(i,j)-rx(i))
!		    fy(i,j) = 5.0*(ey(i,j)-ry(i))
!		    fz(i,j) = 5.0*(ez(i,j)-rz(i))
!	       Enddo
!	    Enddo
!        !--------writing vmd file------------------------------
!	   open(12, file = 'running2.xyz', ACCESS = 'append')
!	   write(12,*) Nt
!	   write(12,*)
!	   write(12,1100) C(1), box(1),box(2),box(3)
!	   write(12,1100) C(2), box(4),box(5),box(6)
!	   write(12,1100) C(3), box(7),box(8),box(9)
!	   write(12,1100) C(4), box(10),box(11),box(12)
!	   write(12,1100) C(5), box(13),box(14),box(15)
!	   write(12,1100) C(6), box(16),box(17),box(18)
!	   write(12,1100) C(7), box(19),box(20),box(21)
!	   write(12,1100) C(8), box(22),box(23),box(24)
!	   Do i = 1, N
!	      write(12,1100) 'F', rx(i), ry(i), rz(i)
!	      do j = 1, Np
!		write(12, 1100) 'He', ex(i,j), ey(i,j), ez(i,j)
!!		   write(12, 1100) 'He2' , rx(i)+fx(i,j)/20.0,
!!     &               ry(i)+fy(i,j)/20.0, rz(i)+fz(i,j)/20.0
!	      enddo
!	   Enddo
!-------------------------------------------------------------------------------


1200	Format(2X, A2, 3(F14.7, 2X)) 
1100	Format(2X, 1A, 3(F14.7, 2X))	   
	!--------------------------------------------------------

	!-------------------writing finalconfig.txt file--------------
!	open(13, file ='finalconfig.txt')
!	Do i = 1, N
!	   write(13,*) rx(i), ry(i), rz(i)
!	   do j = 1, Np
!		write(13, *) fx(i,j), fy(i,j), fz(i,j)
!	   enddo
!	Enddo
	!---------------------------------------------------------

	Stop
	End program rdf




