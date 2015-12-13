
c     Routine for cell classification

      subroutine genib(xs,ys,xnxs,xnys,xc,yc,dist,distg,hh,weight,acoef,bcoef,
     c                 nlist,itag,lpri,idif,jdif,ii,jj,nlistmx,nintpts,nobjs,power)
      implicit real(a-h,o-z)

c ---- real variables passed from outside

      dimension xs(nlistmx,nobjs),ys(nlistmx,nobjs)
      dimension xnxs(nlistmx,nobjs),xnys(nlistmx,nobjs)
      dimension xc(1:ii+1,1:jj+1)
      dimension yc(1:ii+1,1:jj+1)
      dimension dist(1:ii+1,1:jj+1,nobjs)
      dimension distg(1:ii+1,1:jj+1)
      dimension hh(1:ii+1,1:jj+1)
      dimension weight(1:ii+1,1:jj+1,nintpts)
      dimension acoef(1:ii+1,1:jj+1)
      dimension bcoef(1:ii+1,1:jj+1)

c ---- integer variables passed from outside

      dimension nlist(nobjs)
      dimension idif(nintpts), jdif(nintpts)
      dimension itag(1:ii+1,1:jj+1,nobjs)
      dimension lpri(1:ii+1,1:jj+1)

c ---------------------------------------------------------
c ---- Step 1:  Classify Cells
c ---------------------------------------------------------

c ---- Step 1a:  determine local signed distance

      do n=1,nobjs
      do j=1,jj+1
      do i=1,ii+1
        dist(i,j,n) = 1e12
        do k=1,nlist(n)
          dd = (xc(i,j)-xs(k,n))**2 + (yc(i,j)-ys(k,n))**2
          if(dd.lt.dist(i,j,n)) then
           itag(i,j,n) = k   !tag for nearest surface point
           dist(i,j,n) = dd        
          endif
        enddo
        dist(i,j,n) = sqrt(dist(i,j,n))
        dotp = (xc(i,j)-xs(itag(i,j,n),n))*xnxs(itag(i,j,n),n) 
     c       + (yc(i,j)-ys(itag(i,j,n),n))*xnys(itag(i,j,n),n)
        dist(i,j,n) = dist(i,j,n)*sign(1.0,dotp)
      enddo
      enddo
      enddo


c ---- Step 1b:  determine global signed distance

      do j=1,jj+1
      do i=1,ii+1
        distg(i,j) = 1e6
        lpri(i,j) = 1
        do n=1,nobjs
          distg(i,j) = min(distg(i,j),dist(i,j,n))
          if(distg(i,j).eq.dist(i,j,n)) lpri(i,j) = n
        enddo
      enddo
      enddo

c ---- Step 1c:  determine Heaviside function

      do j=2,jj
      do i=2,ii
        icflg = 0
        do k=1,nintpts
          iq = i+idif(k)
          jq = j+jdif(k)
          if(distg(i,j).gt.0.0 .and. distg(iq,jq).lt.0.0) icflg = 1
        enddo
        if(icflg.eq.1) hh(i,j) = 1.0
        if(distg(i,j).lt.0.0) hh(i,j) = 1.0
      enddo
      enddo

      hh(1,:) = hh(2,:)
      hh(ii+1,:) = hh(ii,:)
      hh(:,1) = hh(:,2)
      hh(:,jj+1) = hh(:,jj)

c -------------------------------------------------------------
c ---- now you have four arrays that define cell classification
c
c      distg(i,j) = global signed distance function
c      hh(i,j) = Heaviside function (1.0 if interior or band; 0.0 otherwise)
c      lpri(i,j) = tag that tells which object is closest
c      itag(i,j,l) = index for nearest point on object 'l'
c -------------------------------------------------------------

c ---------------------------------------------------------
c ---- Step 2:  Determine interpolation data
c ---------------------------------------------------------
	! deldsum1 is same as deld1
	! deldsum1d is same as deld2
	
       do j=2,jj
       do i=2,ii
        if(hh(i,j).eq.1.0 .and. distg(i,j).gt.0.0) then
        
			weight(i,j,1:nintpts) = 0.0 
	
			xnxd = xnxs(itag(i,j,lpri(i,j)),lpri(i,j))
			xnyd = xnys(itag(i,j,lpri(i,j)),lpri(i,j))
			xcd = xc(i,j)
			ycd = yc(i,j)

			deldsum1 = 0.0
			deldsum1d = 0.0
			do k=1,nintpts

				iq = i + idif(k)
				jq = j + jdif(k)
				xcp = xc(iq,jq)
				ycp = yc(iq,jq)
				dx = xcp - xcd
				dy = ycp - ycd
				distance = sqrt(dx**2 + dy**2)
				deld = dx*xnxd + dy*xnyd
				dcross = sqrt(distance**2-deld**2)

				if(deld.gt.0.0.and.hh(iq,jq).eq.0.0) then   ! consider only field cells
					dcross = 1.0/(dcross + 1e-12)
					weight(i,j,k) = dcross
					deldsum1 = deldsum1 + dcross
					deldsum1d = deldsum1d + dcross*deld
				endif

			enddo

			if(deldsum1.eq.0.0) then                    ! consider band and field cells

				deldsum1 = 0.0
				deldsum1d = 0.0
				do k=1,nintpts
					iq = i + idif(k)
					jq = j + jdif(k)
					xcp = xc(iq,jq)
					ycp = yc(iq,jq)
					dx = xcp - xcd
					dy = ycp - ycd
					distance = sqrt(dx**2 + dy**2)
					deld = dx*xnxd + dy*xnyd
					dcross = sqrt(distance**2-deld**2)

					if(deld.gt.0.0.and.hh(iq,jq).ge.0.0) then   
						dcross = 1.0/(dcross + 1e-12)
						weight(i,j,k) = dcross
						deldsum1 = deldsum1 + dcross
						deldsum1d = deldsum1d + dcross*deld
					endif

				enddo

			endif

			if(deldsum1.ne.0.0) then

				weight(i,j,1:nintpts) = weight(i,j,1:nintpts)/deldsum1
				deld = deldsum1d/deldsum1

			else

				write(6,*) 'no interpolation point found'
				weight = 0.0
				deld = 0.0

			endif


			dratio = deld/distg(i,j)
			acoef(i,j) = (1.0/(1.0+dratio))**power
			diste = distg(i,j) + 0.5*deld
			distd = 0.5*distg(i,j)
			term = (distd/diste)**power/dratio
			bcoef(i,j) = term/(1.0 + term)             ! attempts to satisfy continuity
	c       bcoef(i,j) = distg(i,j)/(deld+distg(i,j))  ! linear behavior

        endif

       enddo
       enddo

       return
       end
