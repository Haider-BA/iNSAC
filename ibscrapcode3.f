! ------------------------------------------------------------
!      add IB terms (do this after you've computed residual)
! ------------------------------------------------------------

       if(ibadd.eq.1) then      ! input flag to add or not add IB stuff

       do j=2,jj
       do i=2,ii

        if(hh(i,j).eq.1.0) then

         if(distg(i,j).le. 0.0) then   ! interior cell

           ubc = 0.0
           vbc = 0.0
           pbc = 0.0

         else                          ! band cell

         uave = 0.0
         vave = 0.0
         pave = 0.0
         do k=1,nintpts

           iq = i+idif(k)
           jq = j+jdif(k)

           uave = uave + weight(i,j,k)*u(iq,jq)
           vave = vave + weight(i,j,k)*v(iq,jq)
           pave = pave + weight(i,j,k)*p(iq,jq)

         enddo

         xnxd = xnxs(itag(i,j,lpri(i,j)),lpri(i,j))
         xnyd = xnys(itag(i,j,lpri(i,j)),lpri(i,j))
         udotn = uave*xnxd + vave*xnyd

         unor = udotn*xnxd
         vnor = udotn*xnyd
         utan = uave - unor
         vtan = vave - vnor

         ubc = utan*acoef(i,j) + unor*bcoef(i,j)
         vbc = vtan*acoef(i,j) + vnor*bcoef(i,j)
         pbc = pave

         endif

         dtv = vol(i,j)/dt(i,j)                  ! dt is your (local) time step
         res(i,j,1) = dtv*(p(i,j)-pbc)/beta(i,j)**2 ! beta is your AC parameter
         res(i,j,2) = dtv*rho*(u(i,j)-ubc) 
         res(i,j,3) = dtv*rho*(v(i,j)-vbc) 

       endif

       enddo
       enddo

       endif

