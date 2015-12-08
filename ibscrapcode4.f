
c ---- update to get provisional velocity

      do j=2,jj
      do i=2,ii
        vdtl = vol(i,j)/deltat(i,j)  !local 
        vdtg = vol(i,j)/physdt       !physical
        bta2 = beta(i,j)**2
        t11 = ((1.0-hh(i,j))*vdtl + hh(i,j)*vdtg)/bta2
        t22 = ((1.0-hh(i,j))*vdtl + vdtg)*rho

        deltap = -res(i,j,1)/t11
        deltau = -res(i,j,2)/t22
        deltav = -res(i,j,3)/t22

        p(i,j) = p(i,j) + deltap
        u(i,j) = u(i,j) + deltau
        v(i,j) = v(i,j) + deltav
      enddo
      enddo

c ---- IB additions to residual: note that physical time step is used here

         dtv = vol(i,j)/physdt
         res(i,j,1) = dtv*(p(i,j)-pbc)/bta**2
         res(i,j,2) = dtv*rho*(u(i,j)-ubc) 
         res(i,j,3) = dtv*rho*(v(i,j)-vbc) 

