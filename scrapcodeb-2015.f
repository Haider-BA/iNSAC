! A few code scraps (PLEASE CHECK FOR ERRORS!!!!)

! ------------------------------------------------------------
!      "j" direction flux balance (parallelogram CV)
! ------------------------------------------------------------

       do i=2,ii

       do j=1,jj
        f1 = u(i,j+1)
        f3 = u(i,j)
        f2 = 0.25*(u(i,j)+u(i,j+1)+u(i+1,j+1)+u(i+1,j))
        f4 = 0.25*(u(i,j)+u(i,j+1)+u(i-1,j+1)+u(i-1,j))
        volave = 0.5*(vol(i,j) + vol(i,j+1)) 
        deln(1) = yc(i,j+1) - yc(i,j)
        deln(2) = -(xc(i,j+1)-xc(i,j))
        gradf(1:2) = (f1-f3)*area(i,j,2,1:2) + (f2-f4)*deln(1:2)     ! Green's theorem
        gradf(1:2) = gradf(1:2)/volave                               !divide by average volume
        h(j) = gradf(1)*area(i,j,2,1) + gradf(2)*area(i,j,2,2)       !grad u * n * area
       enddo

       do j=2,jj
         res(i,j) = res(i,j) + h(j)-h(j-1)    ! note augmentation of residual
       enddo

       enddo

! ---------------------------------------------------------------------------
!      "j" direction flux balance (normal / tangential gradient decomposition
! ---------------------------------------------------------------------------

       do i=2,ii

       do j=1,jj
        f1 = u(i,j+1)
        f3 = u(i,j)

        delt(1) = xc(i,j+1) - xc(i,j)          
        delt(2) = yc(i,j+1) - yc(i,j)
        dels = sqrt(delt(1)**2 + delt(2)**2)   
        tau(1) = delt(1)/dels                             !unit vector pointing from j to j+1
        tau(2) = delt(2)/dels

        areaface = sqrt(area(i,j,2,1)**2+area(i,j,2,2)**2)
        xn(1) = area(i,j,2,1)/areaface                    !normal vector to mesh face
        xn(2) = area(i,j,2,2)/areaface
        taudotn = tau(1)*xn(1) + tau(2)*xn(2)             !dot product of tau and n

        gradave(1:2) = 0.5*(grad(i,j,1:2) + grad(i,j+1,1:2)) !average gradient vector
        graddottau = (gradave(1)*tau(1) + gradave(2)*tau(2))*taudotn
        gradf(1:2) = (f1-f3)*taudotn*xn(1:2)/dels + gradave(1:2) - graddottau*xn(1:2) !gradient
        h(j) = gradf(1)*area(i,j,2,1) + gradf(2)*area(i,j,2,2)   !grad u * n * area
       enddo

       do j=2,jj
         res(i,j) = res(i,j) + h(j)-h(j-1)    ! note augmentation of residual
       enddo

       enddo

! -----------------------------------------------------------------------------
!      Matrix calculations  (based on 'thin layer form' of parallel CV model
!      Note that you only have to do this once (before the iterations)
! ----------------------------------------------------------------------------

      do j=2,jj

      do i=1,ii
        volave = 0.5*(vol(i,j) + vol(i+1,j)) 
        ahalf(i) = (area(i,j,1,1)**2 + area(i,j,1,2)**2)/volave
      enddo

      do i=2,ii
        a(i,j,1) = ahalf(i-1)
        a(i,j,5) = ahalf(i)
        a(i,j,3) = -(ahalf(i)+ahalf(i-1))
      enddo

      enddo
 
      do i=2,ii
 
      do j=1,jj
        volave = 0.5*(vol(i,j) + vol(i,j+1))
        bhalf(j) = (area(i,j,2,1)**2 + area(i,j,2,2)**2)/volave
      enddo

      do j=2,jj
        a(i,j,2) = bhalf(j-1)
        a(i,j,4) = bhalf(j)
        a(i,j,3) = a(i,j,3) -(bhalf(j)+bhalf(j-1))
      enddo

      enddo


! ------------------------------------------------------------
!      Cell centered gradients using Green's theorem
! ------------------------------------------------------------

      do j=2,jj
      do i=2,ii
        f1 = 0.5*(u(i,j)+u(i+1,j))
        f2 = 0.5*(u(i,j)+u(i,j+1))
        f3 = 0.5*(u(i,j)+u(i-1,j))
        f4 = 0.5*(u(i,j)+u(i,j-1))
        grad(i,j,:) = f1*area(i,j,1,:) - f3*area(i-1,j,1,:) + f2*area(i,j,2,:) - f4*area(i,j-1,2,:)
        grad(i,j,:) = grad(i,j,:)/vol(i,j)
      enddo
      enddo
       

