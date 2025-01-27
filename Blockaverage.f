!-------------------------------------------------------------   
        Program blockav
        implicit none

        integer nblock, i, j, ndata, ii, opbl
	Parameter (ndata = 2)
        double precision av(ndata), sav(ndata), bdata(100000, ndata)
        double precision arav(ndata), step
	
	
        open (22, file='fort.51')
	nblock = 100000
        Do j=1, ndata
           arav(j) = 0.0
        Enddo
        Do i=1, nblock
           Read(22,*) step, (bdata(i,j), j=1, ndata)
           Do ii=1, ndata
              arav(ii) = arav(ii) + bdata(i,ii)
           Enddo
        Enddo
        Do j=1, ndata
           arav(j) = arav(j)/real(nblock)
        Enddo
        Do While (nblock .ge. 4)
           nblock = nblock/2
           i = 1
           Do j = 1, ndata
              av(j) = 0
              sav(j) = 0
           Enddo
           Do ii = 1, nblock
              Do j = 1, ndata
               bdata(ii, j) = ( bdata(i,j) + bdata(i+1,j) ) /2.0
               av(j) = av(j) + bdata(ii, j)
               sav(j) = sav(j) + bdata(ii, j)*bdata(ii, j)
              Enddo
              i = i + 2
           Enddo
           Do j = 1, ndata
              av(j) = av(j)/nblock
              sav(j) = (sav(j)/nblock) - av(j)*av(j)
           Enddo
	      print *, av(2), SQRT(sav(2))
           Do j = 1, ndata
              av(j) = SQRT(sav(j)/(nblock-1))
              sav(j) = av(j)/SQRT(2.*(nblock-1.))
           Enddo
           opbl = opbl + 1
           write (19, 90000) opbl, (arav(j), av(j), sav(j), j=1, ndata)
!            write (19, 90000) opbl, (av(j), sav(j), j=1, ndata)

       Enddo
90000  Format (1x, i9, 30(1x,f8.4))
       Stop
       End
