c ---- Note: ii = imx, jj = jmx, standard cell-centered FV ordering

c ---- open file containing immersed boundaries

      nintpts = 8   ! number of neighbors of a cell

c     open(10,file='ibsurf.dat',status='old')

c     read(10,*) nobjs
c     nlistmx = 0.0
c     do n=1,nobjs
c       read(10,*) nlist(n)
c       nlistmx = max(nlistmx,nlist(n))
c     enddo

c     allocate(xs(nlistmx,nobjs))
c     allocate(ys(nlistmx,nobjs))
c     allocate(xnxs(nlistmx,nobjs))
c     allocate(xnys(nlistmx,nobjs))

c     do n=1,nobjs
c     do m=1,nlist(n)
c       read(10,*) xs(m,n),ys(m,n),xnxs(m,n),xnys(m,n)
c     enddo
c     enddo

c     close(10)

c     Here, I set up a couple of circles as IBs

      nobjs = 2
      allocate(nlist(nobjs))
      nlistmx = 2001
      power = 1.0/1.0
      nlist(1) = 2001
      nlist(2) = 2001

      allocate(xs(nlistmx,nobjs))
      allocate(ys(nlistmx,nobjs))
      allocate(xnxs(nlistmx,nobjs))
      allocate(xnys(nlistmx,nobjs))

      rad = 0.125
      xcnt = 0.75
      ycnt = 0.75
      pi = 4.0*atan(1.0)
      do k=1,2001
       thet = 2.0*pi*float(k-1)/2000.0
       xs(k,1) = rad*cos(thet) + xcnt
       ys(k,1) = rad*sin(thet) + ycnt
       xnxs(k,1) = cos(thet)
       xnys(k,1) = sin(thet)
      enddo

      rad = 0.125
      xcnt = 0.5
      ycnt = 0.5
      pi = 4.0*atan(1.0)
      do k=1,2001
       thet = 2.0*pi*float(k-1)/2000.0
       xs(k,2) = rad*cos(thet) + xcnt
       ys(k,2) = rad*sin(thet) + ycnt
       xnxs(k,2) = cos(thet)
       xnys(k,2) = sin(thet)
      enddo

c ---- variable allocations

      allocate(dist(1:ii+1,1:jj+1,nobjs))
      allocate(distg(1:ii+1,1:jj+1))
      allocate(hh(1:ii+1,1:jj+1))
      allocate(weight(1:ii+1,1:jj+1,nintpts))
      allocate(acoef(1:ii+1,1:jj+1))
      allocate(bcoef(1:ii+1,1:jj+1))
      allocate(itag(1:ii+1,1:jj+1,nobjs))
      allocate(lpri(1:ii+1,1:jj+1))
      allocate(idif(nintpts))
      allocate(jdif(nintpts))

c ---- index-shifts

      idif(1) = -1
      jdif(1) = 1
      idif(2) = 0
      jdif(2) = 1
      idif(3) = 1
      jdif(3) = 1
      idif(4) = 1
      jdif(4) = 0
      idif(5) = 1
      jdif(5) = -1
      idif(6) = 0
      jdif(6) = -1
      idif(7) = -1
      jdif(7) = -1
      idif(8) = -1
      jdif(8) = 0

c ---- call statement for IB routine (do this AFTER you've done all processing
c      of the mesh)

      call genib(xs,ys,xnxs,xnys,xc,yc,dist,distg,hh,weight,acoef,bcoef,
     c                 nlist,itag,lpri,idif,jdif,ii,jj,nlistmx,nintpts,nobjs,power)

