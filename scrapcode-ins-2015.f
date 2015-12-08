
!! Code scraps:  incompressible Navier-Stokes solver: 2015 
!! Dr Jack R. Edwards Jr.


! ------------------------------------------------------------
!      "j" direction flux balance (inviscid only) 
!
!      Note: beta2 stores beta**2
!      Note: const scales dissipation terms
!      Note: rho is density
! ------------------------------------------------------------

       do i=2,ii   

       do j=1,jj

        area = sqrt(del(i,j,2,1)**2 + del(i,j,2,2)**2)
        xnx = del(i,j,2,1)/area
        xny = del(i,j,2,2)/area
        uave = 0.5*(u(i,j)+u(i,j+1))
        vave = 0.5*(v(i,j)+v(i,j+1))
        pave = 0.5*(p(i,j)+p(i,j+1))
        udotn = uave*xnx + vave*xny
        btave2 = 0.5*(beta2(i,j)+beta2(i,j+1))
        eiga = 0.5*(abs(udotn) + sqrt(udotn**2 + 4.0*btave)
        delp = p(i,j)-p(i,j+1)                 ! Note: first order model
c       delp = pleft-pright                    ! higher-order model
        delu = u(i,j)-u(i,j+1)
        delv = v(i,j)-v(i,j+1)

        dissip  = 0.5*area*const*eiga*delp/btave2   ! multiply by 0.5 to make const < 1.0
        xmassflux = area*rho*udotn 
        h(j,1) = xmassflux + dissip
        h(j,2) = xmassflux*uave + xnx*area*pave 
        h(j,3) = xmassflux*vave + xny*area*pave 

       enddo

c      add viscous flux calculation here 

       do j=1,jj

       h(j,2) = h(j,2) - 
       h(j,3) = h(j,3) - 

       enddo

       do j=2,jj
         res(i,j,1) = res(i,j,1) + h(j,1) - h(j-1,1)                         
         res(i,j,2) = res(i,j,2) + h(j,2) - h(j-1,2) 
         res(i,j,3) = res(i,j,3) + h(j,3) - h(j-1,3) 
       enddo

       enddo

c----- minimum modulus averaging

       function xmd(a,b) = 0.5*(sign(1.0,a)+sign(1.0,b))*min(abs(a),abs(b))

c----- left and right states using minimum modulus averaging ('j' direction)

       do j=2,jj
         du(j) = xmd(p(i,j+1)-p(i,j),p(i,j)-p(i,j-1))
       enddo
       du(1) = du(2)
       du(jj+1) = du(jj)

       pleft = p(i,j) + 0.5*du(j)
       pright = p(i,j+1) - 0.5*du(j+1)

c ----- local time step (average areas to cell centers)

       do j=2,jj
       do i=2,ii
        a11ave = 0.5*(del(i,j,1,1) + del(i-1,j,1,1))
        a12ave = 0.5*(del(i,j,1,2) + del(i-1,j,1,2))
        b11ave = 0.5*(del(i,j,2,1) + del(i,j-1,2,1))
        b12ave = 0.5*(del(i,j,2,2) + del(i,j-1,2,2))
        areaaave = sqrt(a11ave**2 + a12ave**2)
        areabave = sqrt(b11ave**2 + b12ave**2)
        vdotna = (u(i,j)*a11ave + v(i,j)*a12ave)/areaaave
        vdotnb = (u(i,j)*b11ave + v(i,j)*b12ave)/areabave
        eiga = 0.5*areaaave*(abs(vdotna) + sqrt(vdotna**2 + 4.0*beta2(i,j)))
        eigb = 0.5*areabave*(abs(vdotnb) + sqrt(vdotnb**2 + 4.0*beta2(i,j)))
        voldt = eiga + eigb
        dt(i,j) = cfl*vol(i,j)/voldt     !cfl is cfl number
       enddo
       enddo
        
c ---- solution update

      do j=2,jj
      do i=2,ii
        dtv = dt(i,j)/vol(i,j)
        deltap = -beta2(i,j)*dtv*res(i,j,1)
        deltau = -dtv*res(i,j,2)/rho
        deltav = -dtv*res(i,j,3)/rho
        p(i,j) = p(i,j) + deltap
        u(i,j) = u(i,j) + deltau
        v(i,j) = v(i,j) + deltav
      enddo
      enddo

! Don't have here:
! - Zeroing out of mass-flux terms in inviscid flux
! - Boundary conditions, ie, ghost cell values
! - Calculation of beta
! - Viscous fluxes
