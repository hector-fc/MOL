C     ==========================================
C     
C     A simple program  Runge Kutta
C     MOL:  Method of line
C     Hector FC
C      
C     ==========================================
      
      program  molpack
      implicit none
      integer nmax
            
      parameter (nmax = 10000)
C     scalars
      double precision dx,dt,tmax,cf
      integer n,iprint
      
      
C     array
      double precision x(nmax),u(nmax),v(nmax),unew(nmax),vnew(nmax)
      
C     local arrays
      real dum(2)
      real time      
C     data statements
      data dum /0.0, 0.0/ 
      
      open(12,FILE='molpack.out')
C     *******************************************
C     Initialization 
C     *******************************************      
      call ini(x,u,v,dx,dt,tmax,n,iprint)      
      cf = dt/(dx*dx)      
      write(*,fmt=1000) dx,dt,cf
      
C     *******************************************
C     call solver 
C     *******************************************
      time = dtime(dum)      
      call rkmet(x,u,unew,v,vnew,dx,dt,tmax,n,iprint)
C     call imrkmet(x,u,unew,v,vnew,dx,dt,tmax,n,iprint)

      time = dtime(dum)
      time = dum(1)
      write(*, fmt=1010) time
      write(12, fmt=1010) time
      write(12, 1020)  dx, dt, tmax, cf      
      close(12)
      
C     Non -executable statements
      
 1000 format(/,1X,'dx:=',F8.2,2X,'dt:=',D8.2,2X,'cf:=',F8.2)      
 1010 format(/,1X,'Total CPU time in seconds:',F8.2)
 1020 format(/,1X,'dx=',F8.4,1X,'dt=',F8.4,1X,'T=',F8.4,1X,'cf=',F8.2) 
      
      stop      
      end
      
C     ==========================================
C     ==========================================
            
      subroutine  ini(x,u,v,dx,dt,tmax,n,iprint)
      implicit none
C     scalar
      integer n, iprint
      double precision dx,dt,tmax,uini,vini
C     arrays
      double precision x(n),u(n),v(n)
      
C     scalar
      double precision a,b,rho,leta1,leta2,reta1,reta2,lx0,rx0,
     &     lele,rele,lc,rc
      
      integer  i
      character*36  Label
      character(len=4) soliton
      
      common /noise/ rho
      common /probp/ leta1,leta2,reta1,reta2,lx0,rx0,lele,rele,
     &     lc,rc
      
C     *******************************************
      
      open(8,FILE='molpack.ini',FORM='FORMATTED')
      
      read(8,*) Label, n
      write(6,*) Label,'n:      ',n
      read(8,*) Label, a
      write(6,*) Label,'a:      ',a
      read(8,*) Label, b
      write(6,*) Label,'b:      ',b
      read(8,*) Label, tmax
      write(6,*) Label,'tmax:   ',tmax
      read(8,*) Label, dt
      write(6,*) Label,'dt:     ',dt
      read(8,*) Label, rho
      write(6,*) Label,'rho:    ',rho      
      read(8,*) Label, iprint
      write(6,*) Label,'info:  ',iprint
      read(8,*) Label, soliton 
      write(6,*) Label,'soliton:',soliton 
      close(8)
                
      if (soliton.eq.'K') then
         leta1 =  1
         leta2 =  1
         lx0   = -12
         lele  = -1 
         lc    = 0.5
         reta1 = 1                
         reta2 = 0                
         rx0   = 12
         rele  = 0.0
         rc    = 0.0
      else if(soliton.eq.'KaK') then
         leta1 =  1
         leta2 =  1
         lx0   = -12
         lele  = -1 
         lc    = 0.5
         reta1 = 1                
         reta2 = -1                
         rx0   = 12
         rele  = 0.0
         rc    =-0.5
      else if(soliton.eq.'KK') then     
         leta1 =  1
         leta2 =  1
         lx0   = -12
         lele  = -1 
         lc    = 0.5
         reta1 = 1                
         reta2 = 1                
         rx0   = 12
         rele  = 0.0
         rc    =-0.5
      else if(soliton.eq.'aKK') then     
         leta1 =  -1
         leta2 =  1
         lx0   = -12
         lele  = -1 
         lc    = 0.5
         reta1 = 1                
         reta2 = 1                
         rx0   = 12
         rele  =-1.0
         rc    =-0.5                    
      else if(soliton.eq.'aKaK') then     
         leta1 =  1
         leta2 = -1
         lx0   = -12
         lele  =  1 
         lc    = 0.5
         reta1 = 1                
         reta2 = -1                
         rx0   = 12
         rele  = 0.0
         rc    =-0.5                    
      else if(soliton.eq.'sKaK') then     
         leta1 = -1
         leta2 =  1
         lx0   = -12
         lele  = -1 
         lc    = 0.5
         reta1 = 1                
         reta2 = 1                
         rx0   = 12
         rele  =-1.0
         rc    =-0.5                                     
      end if
      
C     ********************************************
      
      dx = abs(b-a)/dble(n-1)
            
      x(1) = a 
      u(1) = uini( 0.0d0, x(1))
      v(1) = vini( 0.0d0, x(1))      
      do i = 2,n
         x(i) = x(i-1 ) + dx 
         u(i) = uini(0.0d0, x(i))
         v(i) = vini(0.0d0, x(i))
      enddo

C     ******************************************
C     files 
      open(11,FILE='energy.txt',STATUS='REPLACE')     
      close(11)

            
      return
      end
                 
C     =========================================
      subroutine  imrkmet(x,u,unew,v,vnew,dx,dt,tmax,n,iprint)      
      implicit none
      
C     scalar
      integer n,iprint
      double precision dt,dx,tmax 
      
C     arrays 
      double precision x(n), u(n),unew(n),v(n),vnew(n),
     &     teng(n),tpot(n),qp(n),qn(n),gp(n),gn(n)  
      
C     local arrays 
      double precision um1(n),um2(n),um3(n), vm1(n),vm2(n),vm3(n) 
      
C     local variables
      double precision tc,err, errinf
      integer iter, i, iplot, info, csol, eplot, qplot, gplot 
      
      iplot = 5000
      eplot = 6000
      qplot = 7000
      gplot = 8000
      
C     arrays 
      iter = 0
      tc  = 0.0 
      csol  = int(0.5*tmax/dt) 
      
C     plot initial value                                 
C      call plotc(x, u, v, n, iplot)
C     iplot = iplot + 1
      
C     first call rk4
      call rk4(u,um3,v,vm3,dx,dt,n,info)            
      do i =1,n
         u(i) = um3(i)
         v(i) = vm3(i)
      enddo
      tc = tc + dt
      iter =iter + 1
      
C     second call rk4
      call rk4(u,um2,v,vm2,dx,dt,n,info)            
      do i =1,n
         u(i) = um2(i)
         v(i) = vm2(i)
      enddo
      tc = tc + dt
      iter = iter + 1
      
C     thrid call rk4
      call rk4(u,um1,v,vm1,dx,dt,n,info)
      do i =1,n
         u(i) = um1(i)
         v(i) = vm1(i)
      enddo
      tc = tc + dt
      iter = iter + 1
      
C     four call rk4
      call rk4(u,unew,v,vnew,dx,dt,n,info)
      do i =1,n
         u(i) = unew(i)
         v(i) = vnew(i)
      enddo
      tc = tc + dt 
      iter =iter + 1
      
C     four order method

      call plotc(x, u, v, n, iplot)
      iplot = iplot + 1
      
      call energy(tc,u,v,teng,tpot,qp,qn,gp,gn, dx,n,info )
      
      call plotc(x, teng, tpot, n, eplot)
      eplot = eplot + 1
      call plotc(x, qp, qn, n, qplot)
      qplot = qplot + 1
      call plotc(x, gp, gn, n, gplot)
      gplot = gplot + 1
      
      do while (tc.lt.tmax)         
         
                                    
         call  imrk4(u,um1,um2,um3,unew,v,vm1,vm2,vm3,vnew,
     &        dx,dt,n,info)

         
         if (iprint.eq.0 ) then
            if( mod(iter,50).eq.0) then
               err = errinf( tc, x, u, n)               
            endif
         endif
         
         if (iter.eq.csol) then            
            call plotc(x, u, v, n, iplot)
            iplot = iplot + 1            
            call plotc(x, teng, tpot, n, eplot)
            eplot = eplot + 1
            call plotc(x, qp, qn, n, qplot)
            qplot = qplot + 1
            call plotc(x, gp, gn, n, gplot)
            gplot = gplot + 1            
         endif

         if (mod(iter,10).eq.0 ) then
            call energy(tc,u,v,teng,tpot,qp,qn,gp,gn,dx,n,info )                        
         endif
                  
         do i=1,n
            um3(i) = um2(i)                      
            vm3(i) = vm2(i)
            
            um2(i) = um1(i)                      
            vm2(i) = vm1(i)
            
            um1(i) = u(i)                      
            vm1(i) = v(i)
            
            u(i) = unew(i)                      
            v(i) = vnew(i)            
         enddo
         
         tc = tc + dt 
         iter = iter + 1
         
      enddo

      
      call plotc(x, u, v, n, iplot)
      call plotc(x, teng, tpot, n, eplot)
      call plotc(x, qp, qn, n, qplot)
      call plotc(x, gp, gn, n, gplot)
      
      return
      end

      
C     ==========================================      
      subroutine  rkmet(x,u,unew,v,vnew,dx,dt,tmax,n,iprint)
      implicit none
      external  evalf
C     scalar
      integer n,iprint
      double precision dt,dx,tmax 
      
C     arrays 
      double precision x(n), u(n),unew(n),v(n),vnew(n),
     &     teng(n),tpot(n),qp(n),qn(n),gp(n),gn(n)  
      
C     local variables
      double precision tc,err, errinf
      integer iter, i, iplot, info, csol,eplot, qplot, gplot 
      
      iplot = 5000
      eplot = 6000
      qplot = 7000
      gplot = 8000
      
C     arrays 
      iter = 0
      tc  = 0.0 
      csol  = int(0.5*tmax/dt) 
      
C     plot initial value                                 
      call plotc(x, u, v, n, iplot)
      iplot = iplot + 1
      
      call energy(tc,u,v,teng,tpot,qp,qn,gp,gn,dx,n,info )   ! new gamm
      
      call plotc(x, teng, tpot, n, eplot)
      eplot = eplot + 1
      call plotc(x, qp, qn, n, qplot)
      qplot = qplot + 1
      call plotc(x, gp, gn, n, gplot)
      gplot = gplot + 1

      
      do while (tc.lt.tmax)         

         call rk4(u,unew,v,vnew,dx,dt,n,info)

         
         if (iprint.eq.0 ) then
            if( mod(iter,50).eq.0) then
            err = errinf( tc, x, u, n)
            endif
         endif
         
         if (iter.eq.csol) then            
            call plotc(x, u, v, n, iplot)
            iplot = iplot + 1            
            call plotc(x, teng, tpot, n, eplot)
            eplot = eplot + 1
            call plotc(x, qp, qn, n, qplot)
            qplot = qplot + 1
            call plotc(x, gp, gn, n, gplot)
            gplot = gplot + 1            
         endif
         
         if (mod(iter,10).eq.0 ) then
            call energy(tc,u,v,teng,tpot,qp,qn,gp,gn,dx,n,info )            
         endif
         
         do i=1,n
            u(i) = unew(i)                      
            v(i) = vnew(i)            
         enddo
         
         tc = tc + dt 
         iter = iter + 1
         
      enddo

      call plotc(x, u, v, n, iplot)
      call plotc(x, teng, tpot, n, eplot)
      call plotc(x, qp, qn, n, qplot)
      call plotc(x, gp, gn, n, gplot)
      
      return
      end
      
C     ==========================================
      subroutine energy(t,u,v,teng,tpot,qp,qn,gp,gn,dx,n,info) 
      implicit none
C     scalars      
      double precision t,dx  
      integer  n, info
      double precision   u(n),v(n)
            
C     local scalars 
      integer i
      double precision teng(n),tpot(n),qp(n),qn(n),gp(n),gn(n)
      double precision funv, int1, int2, int3, int4, dx2, w20,
     &     dfunv, intgp, intgn
      
      external funv, dfunv
     
      info  = 0
      dx2 = dx*dx 
            
      teng(1) = 0.5*( v(1) )**2 + 0.5*((u(2)-u(1))/dx)**2 +
     &     funv( u(1))

      tpot(1)  = v(1)* ( (u(2)-u(1))/dx)
           
      w20 = ( u(1) - 2.0*u(2) + u(3) )/dx2
      
C     *******************************************      
      qp(1) = 2.0 * ( 2.0*w20 - dfunv(u(1))) * funv(u(1))
      
      qn(1) = 4.0 * ((v(2)-v(1))/dx)*funv(u(1))
      
C     *********************************************
      gp(1) = ( v(1)**2 - ((u(2)-u(1))/dx)**2) * dfunv(u(1))*v(1)
      
      gn(1) = ( v(1)**2 - ((u(2)-u(1))/dx)**2)*      
     &     dfunv(u(1)) * ( (u(2)-u(1))/dx )   
C     *********************************************
      do i = 2,n-1

         w20 = ( u(i-1) - 2.0*u(i) + u(i+1) )/dx2
                           
         teng(i)= 0.5*( v(i))**2 + 0.5*( (u(i+1)-u(i) )/dx)**2         
     &        + funv(u(i))
         
         tpot(i) = v(i) * (u(i+1)-u(i))/dx

C     *************************************         
         qp(i) = 2.0 * ( 2.0*w20 - dfunv(u(i)) ) * funv(u(i))
         
         qn(i) = 4.0 * ( (v(i+1) - v(i-1))/(2*dx)) * funv(u(i))

C     *************************************
         
         gp(i) = ( v(i)**2 - ( (u(i+1)-u(i-1))/(2*dx) ) **2) *         
     &         dfunv(u(i))*v(i)
      
         gn(i)  = ( v(i)**2 - (( u(i+1)-u(i-1))/(2*dx) )**2) *      
     &        dfunv(u(i)) * ((u(i+1)-u(i-1))/(2*dx))   
                  
      enddo

      
      teng(n) =  0.5*(v(n))**2  +
     &     0.5*((u(n) - u(n-1))/dx)**2 + funv(u(n)) 
         
      tpot(n) = v(n) * (u(n)-u(n-1))/dx
      
      w20 = ( u(n-2) - 2.0*u(n-1) + u(n) )/dx2
      
      qp(n) = 2.0 * ( 2.0*w20 - dfunv(u(1)) )*funv(u(n))      
      qn(n) = 4.0 * ( (v(n)-v(n-1))/dx)*4.0*funv(u(n))
      
      gp(n) = ( v(n)**2 - ((u(n)-u(n-1))/dx)**2) * dfunv(u(n))*v(n)      
      gn(n) = ( v(n)**2 - ((u(n)-u(n-1))/dx)**2)*      
     &     dfunv(u(n)) * ( (u(n)-u(n-1))/dx )   
      
            
      call itrap(teng,n,dx,int1)
      call itrap(tpot,n,dx,int2)
      call itrap(qp,n,dx,int3)
      call itrap(qn,n,dx,int4) 
      call itrap(gp,n,dx,intgp) 
      call itrap(gn,n,dx,intgn) 
      
      open(11, FILE='energy.txt',ACCESS='APPEND',FORM='FORMATTED')      
      write(11,*) t, int1, int2, int3, int4, intgp, intgn 
      close(11)
      
      return
      end

      
C     ==========================================
      
      
      subroutine itrap(fvec,n,dx,out)
      implicit none
C     scalars
      integer n
      double precision dx,out
      double precision fvec(n)
C     local scalars
      double precision  temp
      integer i
      
      temp = 0.0
      do i  = 2,n-1
         temp = temp + fvec(i)
      enddo
      
      out = dx*temp + dx*(fvec(1)+fvec(n))*0.5
      
      return
      end
      
C     ==========================================      
      subroutine isimp(fvec,n,dx,out)
      implicit none
C     scalars
      integer n
      double precision dx,out
      double precision fvec(n)
C     local scalars
      double precision  temp
      integer i
      
      temp = 0.0 
      do i  = 2,n-1
         temp = temp + 2.0*fvec(i) + fvec(i+1)
      enddo            
      out = dx*( 2.0*temp +  (fvec(1) - fvec(n)) )/3.0      
      return
      end
      
C     ==========================================


      
      subroutine  plotc(x, u, v, n, uci ) 
      implicit none
C     scalars
      integer n, uci      
C     array
      double precision x(n),u(n),v(n)
C     local scalar
      integer i
      character*36 Label
      
      write(Label,50) uci
 50   format(I4,'.txt')
      
      write(6,*)'Creating a file:',Label
      open(10,FILE=Label,FORM='FORMATTED')
      
      do i = 1,n
         write(10,*)  x(i), u(i), v(i)         
      enddo
      
      close(10)
            
      return
      end
      
C     ==========================================
            
      subroutine imrk4(u,um1,um2,um3,unew,v,vm1,vm2,vm3,vnew,
     &     dx,dt,n,info)
      
      implicit none
C     scalars
      integer n, info
      double precision dt,dx, tol
      double precision u(n),um1(n),um2(n),um3(n),unew(n),
     &     v(n),vm1(n),vm2(n),vm3(n),vnew(n)

      
C     local scalar
      integer i, iiter
C     local arrays

      double precision unew0(n), vnew0(n), f1new0(n), f2new0(n),f1m1(n),
     &     f2m1(n), f1m2(n),f2m2(n), f1m3(n),f2m3(n), f1(n),f2(n)
      
      double precision norm

      iiter = 0
      tol = 0.0001
      norm = 1
      

      
      call funf1(u, v, dx, n, f1)
      call funf2(u, dx, n, f2)
      
      call funf1(um1, vm1, dx, n, f1m1)
      call funf2(um1, dx, n, f2m1)
      
      call funf1(um2, vm2, dx, n, f1m2)
      call funf2(um2, dx, n, f2m2)
      
      call funf1(um3, vm3, dx, n, f1m3)
      call funf2(um3, dx, n, f2m3)
      
         
      do i=1,n
         unew0(i) = u(i) + dt * ( 55.0 * f1(i)
     &        - 59.0*f1m1(i) + 37.0*f1m2(i) - 9.0*f1m3(i))/24.0
         
         vnew0(i) = v(i) + dt*( 55.0 * f2(i)

     &        - 59.0*f2m1(i) + 37.0*f2m2(i) - 9.0*f2m3(i))/24.0            
      enddo

      
C     Iner loop
      
      do while ( norm > tol)
         
         call funf1(unew0, vnew0, dx, n, f1new0)
         call funf2(unew0, dx, n, f2new0)
         
         do i = 1,n
            unew(i) = u(i) + dt*( 9.0 * f1new0(i)
     &           + 19.0 * f1(i)- 5.0 * f1m1(i) + f1m2(i) )/24.0
            
            vnew(i) = v(i) + dt*( 9.0 * f2new0(i)
     &           + 19.0*f2(i) - 5.0*f2m1(i) + f2m2(i) )/24.0
         enddo

         info = 1
         
         call normF(unew,vnew,unew0, vnew0,n,info,norm)

         norm = norm/(2.0*dble(n))
                                             
         do i =1,n
            unew0(i) = unew(i)
            vnew0(i) = vnew(i)
         enddo
         iiter = iiter + 1
         
         if( iiter.ge.1000) then
            
            exit 
         endif         
      enddo
      
      return                     
      end

      
C     ==========================================      
C     ==========================================
      
      subroutine mrk4(u,um1,um2,um3,unew,v,vm1,vm2,vm3,vnew,
     &     dx,dt,n,info)
      
      implicit none
C     scalars
      integer n, info
      double precision dt,dx
      double precision u(n),um1(n),um2(n),um3(n),unew(n),
     &     v(n),vm1(n),vm2(n),vm3(n),vnew(n)

      
C     local scalar
      integer i
C     local arrays

      double precision  f1m1(n),f2m1(n), f1m2(n),f2m2(n),
     &     f1m3(n),f2m3(n), f1(n),f2(n)

      info = 1
      
      call funf1(u, v, dx, n, f1)
      call funf2(u, dx, n, f2)
      
      call funf1(um1, vm1, dx, n, f1m1)
      call funf2(um1, dx, n, f2m1)
      
      call funf1(um2, vm2, dx, n, f1m2)
      call funf2(um2, dx, n, f2m2)
      
      call funf1(um3, vm3, dx, n, f1m3)
      call funf2(um3, dx, n, f2m3)
      
      do i=1,n
         unew(i) = u(i) + dt * ( 55.0 * f1(i)
     &        - 59.0*f1m1(i) + 37.0*f1m2(i) - 9.0*f1m3(i))/24.0
         
         vnew(i) = v(i) + dt*( 55.0 * f2(i)
     &        - 59.0*f2m1(i) + 37.0*f2m2(i) - 9.0*f2m3(i))/24.0            
      enddo
                  
      return                     
      end

      
C     ==========================================                
C     ==========================================
      
      subroutine rk4(u,unew,v,vnew,dx,dt,n,info)
      implicit none
C     scalars
      double precision dx,dt
C     arrays 
      double precision  u(n), unew(n), v(n), vnew(n)
      integer  n, info

C     local arrays
      double precision k1f1(n), k1f2(n), k1f3(n), k1f4(n),
     &     k2f1(n), k2f2(n), k2f3(n), k2f4(n),uk(n),vk(n)

C     local scalar
      
      integer i
      info = 0
      
C     Get four  sample value of the revivative 

C     stage 1

      call funf1(u, v, dx, n, k1f1 )      
      call funf2(u, dx, n, k2f1 )
      
C     stage 2
      do i=1,n
         uk(i) = u(i) +  0.5 * dt * k1f1(i)
         vk(i) = v(i) +  0.5 * dt * k2f1(i)
      enddo
      
      call funf1(uk, vk, dx, n, k1f2)
      call funf2(uk, dx, n, k2f2)
      
C     stage 3
      do i=1,n
         uk(i) = u(i) +  0.5 * dt * k1f2(i)
         vk(i) = v(i) +  0.5 * dt * k2f2(i)
      enddo
                  
      call funf1(uk, vk, dx, n, k1f3)
      call funf2(uk, dx, n, k2f3)
C     stage 4
      
      do i=1,n
         uk(i) = u(i) +  dt * k1f3(i)
         vk(i) = v(i) +  dt * k2f3(i)


      enddo
            
      call funf1(uk, vk, dx, n, k1f4)                
      call funf2(uk, dx,n, k2f4)

C     Combine them to estimate the solution at time  tnew = t + dt 
      
      do i=1,n
         unew(i) = u(i) + dt*( k1f1(i) + 2.0 * k1f2(i) + 2.0*k1f3(i) +
     &        k1f4(i) )/( 6.0D0)
         
         vnew(i) = v(i) + dt*( k2f1(i) + 2.0 * k2f2(i) + 2.0*k2f3(i) +
     &        k2f4(i) )/( 6.0D0)   
      enddo
            
      return
      end
           
C     ==========================================      
      subroutine  funf1( u,v,dx,n,fone)
      implicit none
C     scalars
      double precision dx
      integer n
C     arrays
      double precision u(n),v(n),fone(n) 
C     local scalars
      integer i
      
      fone(1) =  ( -11*u(1) + 18*u(2) -9*u(3) + 2*u(4) )/( 6.0*dx)      
      do i=2,n-1
         fone(i) = v(i)
      end do
      fone(n) =-(-2*u(n-3) +9*u(n-2) -18*u(n-1)+11*u(n))/(6.0*dx)      
      return
      end
      
C     ==========================================      
      subroutine  funf2(u,dx,n,ftwo)
      implicit none
C     scalars
      double precision  dx, dfunv
      integer n
C     arrays
      double precision u(n),ftwo(n)
C     local escalar
      double precision dx2
      integer i

      external dfunv
            
      dx2 = dx*dx
      
      ftwo(1) = (2*u(1) - 5*u(2) + 4*u(3) - u(4))/ dx2  -
     &     dfunv( u(1) )
      
      ftwo(2) = ( 11*u(1) - 20*u(2) + 6.*u(3) + 4*u(3) - u(5))      
     &     /(12.0*dx2)   - dfunv( u(2))
      
      do i = 3,n-2
         ftwo(i)=(-u(i-2)+ 16.0*u(i-1)- 30.0*u(i)+ 16.0*u(i+1)-u(i+2))         
     &        /(12.0*dx2)  -  dfunv(u(i))
         
      enddo
      
      ftwo(n-1)=(-u(n-4) + 4*u(n-3) + 6*u(n-2)- 20*u(n-1)+ 11*u(n))
     &     /(12.0*dx2) - dfunv(u(n-1))
      
      ftwo(n)=(-u(n-3 ) + 4*u(n-2) - 5*u(n-1) + 2*u(n))
     &     /dx2 - dfunv(u(n))
      
C     ************************************************
      
      return
      end

      
C     ==========================================      
C     potencial  function
      
      double precision function funv(x) 
      implicit none
C     scalar
      double precision x,rho      
C     local scalar
      double precision  absin      
      common /noise/ rho
            
C     funv = 2.0-2.0 *cos(2.0*x)
      absin = ( abs( sin( x*0.5) ))**rho      
      funv  = (64.0/( rho**2))*( tan(x*0.5)*( 1.0 - absin))**2
      
      return
      end
      
C     ==========================================
      double precision function dfunv(x) 
      implicit none
      
C     scalars      
      double precision rho,x
      
C     local scalar
      double precision  absin, temp

C     external funv
      common /noise/ rho
           
      absin  = ( abs(sin( x*0.5 )))**rho
      temp   = (1.0 - absin)/( (cos(0.5*x))**2) - rho*absin
      dfunv  = (64.0/(rho**2) )* tan( x*0.5 )*(1.0 - absin )*temp
      
C     fundv =  4.0*sin(2*x) 
      return
      end
      
C     ==========================================     
      double precision function phi(t,x,rho,eta1,eta2,x0,ele,c)
      implicit none               
C     scalars
      double precision t,x,rho,eta1,eta2,x0,ele,c
      
C     local scalars      
      double precision  pi,gam,bgam,temp

      
      pi   = 4.0 * atan(1.0d0)      
      gam  = 1.0 / sqrt(1.0 - c*c)
      bgam = 2.0 * sqrt(2.0)*eta1 * gam *(x - t*c - x0)      
      temp = 0.5 + 0.5 *tanh( bgam )
C     temp = 1.0 -  1.0/(1.0 + exp(2*bgam))      
      phi  = 2.0 * eta2 * asin( (temp**(1.0/rho)) ) + ele*pi      
      return
      end

      
C     ==========================================     
      double precision function dphi(t,x,rho,eta1,eta2,x0,ele,c)
      implicit none
      
      double precision phi
      external phi
      
C     scalars 
      double precision t,rho,x,eta1,eta2,x0,ele,c
      
C     local scalars
      double precision h
      
      h = 0.001
      
      dphi = (phi(t-2*h,x, rho,eta1,eta2,x0,ele,c) -
     &     4 * phi(t-h,x, rho,eta1,eta2,x0,ele,c) +
     &     3 * phi(t, x, rho,eta1,eta2,x0,ele,c))/(2.0*h)
      
      return      
      end
      
C     ==========================================     
      double precision function uini(t,x)
      implicit none

      double precision phi,t,x
      external phi
      
      double precision rho,leta1,leta2,reta1,reta2,lx0,rx0,lele,rele,
     &     lc,rc
      
      common /noise/ rho
      common /probp/ leta1,leta2,reta1,reta2,lx0,rx0,lele,rele,lc,rc

      if (x.lt.0.0) then 
         uini =  phi(t,x,rho,leta1,leta2,lx0,lele,lc )
      else         
         uini =  phi(t,x,rho,reta1,reta2,rx0,rele,rc )
      endif
      
      return
      end

      
C     ==========================================           
      double precision function vini(t,x)
      implicit none
      double precision dphi,t,x
      double precision rho,leta1,leta2,reta1,reta2,lx0,rx0,lele,rele,
     &     lc,rc
      
      external dphi

      common /noise/rho
      common /probp/leta1,leta2,reta1,reta2,lx0,rx0,lele,rele,lc,rc
      
      if(x.lt.0.0) then
         vini = dphi( t,x,rho,leta1,leta2,lx0,lele,lc)
      else
         vini = dphi( t,x,rho,reta1,reta2,rx0,rele,rc)
      endif
      
      return
      end
C     ==========================================     
      double precision function errinf(t,x,u,n)
      implicit none
      double precision t,err,phi,rho
      double precision x(n),u(n)
      double precision leta1,leta2,reta1,reta2,lx0,rx0,lele,rele,lc,rc
      integer i,n
      external phi
      
      common /noise/ rho
      common /probp/leta1,leta2,reta1,reta2,lx0,rx0,lele,rele,lc,rc
                          
      err = 0
C     phi(t,x,rho,eta1,eta2,x0,ele,c)
      do i=1,n         
         err = max( err, abs(u(i) -
     &        phi( t,x(i),rho,leta1,leta2,lx0,lele,lc))
     &        )
      enddo
      
      errinf = err 
      write(*,fmt=2000) err, t
      
 2000 format(1X,D7.1,2X,F10.2)      
      
      return      
      end           
C     ==========================================     

      subroutine  normF(u,v,u0,v0,n,info,norm)
      implicit none
      
      integer  n, info      
      double precision  er1,er2, norm
      double precision u(n),u0(n),v(n),v0(n)
C     local  scalars
      integer i
      
      er1 = 0.0
      er2 = 0.0
      if ( info.eq.0)then
         do i = 1,n
            er1 = max(er1,abs(u(i)-u0(i)))  
            er2 = max(er2,abs(v(i)-v0(i)))
         enddo
         norm = max(er1,er2)
         
      elseif(info.eq.1)then
         do i =1,n
            er1 = er1 + abs(u(i)-u0(i))  
            er2 = er2 + abs(v(i)-v0(i))
         enddo
         norm = er1 +  er2
      endif

      return
      end

      

      
      








      
      
      
