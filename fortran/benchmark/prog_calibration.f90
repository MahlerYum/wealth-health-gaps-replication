Program prog_calibration
    
 use mod_parameters
 use mod_spline
 use sort_all
 
 implicit none

 
 Real(8) :: x(nparam), y(nmom)

 Real(8) :: fr, fe, fc,aaa ,bbb,std, corlogearnsloghealth, fval, inipval(2)
 Real(8) :: thlthrsc, mdc1, mdc2, mdc3, suba11, suba33, mcdcut, weight_mat(nmom,nparam), weightvec(nmom)

 Real(8), allocatable ::xbar(:), xr(:), xe(:), xc(:),temp1(:), temp2(:)
 Real(8), allocatable :: mr(:), me(:), mc(:), tempm(:),  modelmm(:)
 Real(8), allocatable :: txx(:,:), yy(:), xx(:,:),tmoment(:,:),outxx(:,:)
 Real(8), allocatable :: moment(:,:), unit(:,:), ovalue(:), tovalue(:),parar(:)
 Integer, allocatable :: sorindx(:)

 Integer :: i1, j1, k1, i2, j2, k2 , i3, j3, k3,step, dim, momdim, testmode, evalnum, i, exper_run
 
 evalnum = 1
 
 ! Load Weight Matrix
    open(18, file="Weights.txt")
    do i = 1,nmom
        read(18,*) weight_mat(i,:)
    enddo
    close(18)
    
 if (parm_est>0) then
     weightvec = weight_mat(:,parm_est)
 end if
 
 If (calind == 0 ) then ! Obtain the solution based on a given set of parameters
                          
     x=(/nuh_10, nuh_20, nuh_30, nuh_40, nuu_10, nuu_20, nuu_30, nuu_40, nuadd0, xiHSh_10, &
         xiHSh_20, xiHSh_30, xiHSh_40, xiHSu_10, xiHSu_20, xiHSu_30, xiHSu_40, xiCLh_10, xiCLh_20, xiCLh_30, &
         xiCLh_40, xiCLu_10, xiCLu_20, xiCLu_30, xiCLu_40, chi_10, chi_20, sigx0, penre0, betamean0, &
         betadev0, bb0, conp0, psi0, ytHS_s0, ytHS_sq0, ytCL_s0, ytCL_sq0, wagep_HS0, wagep_CL0, y1_HS0, y1_CL0 /)
     
     if (run_experiments == 0) then   
        
        exper_run = 0 
        call sub_main(x,y,evalnum,exper_run)   
        
     else
            write(*, '(x)')
            write(*, '(a)') ' Running all experiments'
         
        do i = 1,15
             exper_run = i
             call sub_main(x,y,evalnum,exper_run)
        enddo
     endif
     
         
 else if (calind == 1) then ! Do a calibration
     dim = nparam
     exper_run = 0

     allocate(xbar(dim), xr(dim), xe(dim), xc(dim), temp1(dim), temp2(dim), &
         txx(dim,dim+1), xx(dim,dim+1), &
         unit(dim,dim), sorindx(dim+1), &
         outxx(dim,dim+1), parar(dim))

     allocate(mr(nmom), me(nmom), mc(nmom), tempm(nmom), &
             tmoment(nmom,dim+1), moment(nmom,dim+1), &
             modelmm(nmom))

     allocate(ovalue(dim+1), tovalue(dim+1))

     If (calstep == 1) then ! Given initial parameter, generate a simplex.                      
         x=(/nuh_10, nuh_20, nuh_30, nuh_40, nuu_10, nuu_20, nuu_30, nuu_40, nuadd0, xiHSh_10, &
         xiHSh_20, xiHSh_30, xiHSh_40, xiHSu_10, xiHSu_20, xiHSu_30, xiHSu_40, xiCLh_10, xiCLh_20, xiCLh_30, &
         xiCLh_40, xiCLu_10, xiCLu_20, xiCLu_30, xiCLu_40, chi_10, chi_20, sigx0, penre0, betamean0, &
         betadev0, bb0, conp0, psi0, ytHS_s0, ytHS_sq0, ytCL_s0, ytCL_sq0, wagep_HS0, wagep_CL0, y1_HS0, y1_CL0 /)
         
         call transfrm(dim,1,x,txx(:,1))
         
         unit = 0.0d0
         Do i1=1,dim
             unit(i1,i1)=0.050d0
         end do

         Do i1=2,dim+1
             txx(:,i1)=txx(:,1)+unit(:,i1-1)
         end do

        open (unit=19,file=calfileplace//"initial_simplex.csv",action="write",status="replace")   
            do i1=1,dim+1
                call transfrm(dim,2,txx(:,i1),outxx(:,i1))
                write (19,"(20000(f0.8,',',:))") outxx(:,i1)
            end do
        close(19)

         Do i1= 1, dim+1
             call transfrm(dim,2,txx(:,i1),x)
             call sub_main(x,y,evalnum,exper_run)
             tmoment(:,i1) = y
             tovalue(i1) = 0.0d0
             
             if (efweight==1) then
                 if (parm_est==0) then
                     tovalue(i1)= sum((tdatase**(-2.0d0))*(y-tdatamm)**2.0d0)
                 else
                     tovalue(i1)= sum(weightvec*(tdatase**(-2.0d0))*(y-tdatamm)**2.0d0)
                 end if                 
             else
                 tovalue(i1)= sum(log(y/tdatamm)**2.0d0)
             end if
             
             open (unit=19,file=calfileplace//"moment.csv",action="write",status="replace")
             Do j1=1,dim+1
                 write (19,"(20000(f0.16,',',:))") tmoment(:,j1)
             end do
             close(19)
             open (unit=19,file=calfileplace//"objective fnc value.csv",action="write",status="replace")
             write (19,"(20000(f0.16,',',:))") tovalue
             close(19)
         end do

     else if (calstep == 2) then ! Read a stored set of moments and parameters.
         open (unit=19,file=calfileplace//"moment.csv",status="old")
         Do i1=1,dim+1
             read (19,*) tmoment(:,i1)
         end do
         close(19)

        open (unit=19,file=calfileplace//"parameters.csv",status="old")   
        Do i1=1,dim+1
            read (19,*) txx(:,i1)
        end do
        close(19)
        Do i1=1,dim+1
            call transfrm(dim,1,txx(:,i1),temp1)
            txx(:,i1)=temp1
        end do


         Do i1=1,dim+1
             modelmm=tmoment(:,i1)
             
             tovalue(i1)= 0.0d0
             if (efweight==1) then
                 if (parm_est==0) then                     
                     tovalue(i1) = sum((tdatase**(-2.0d0))*(modelmm-tdatamm)**2.0d0)                        
                 else
                     tovalue(i1) = sum(weightvec*(tdatase**(-2.0d0))*(modelmm-tdatamm)**2.0d0)   
                 end if
                 
             else
                 tovalue(i1) = sum(log(modelmm/tdatamm)**2.0d0)
             end if
             
         end do

         open (unit=19,file=calfileplace//"objective fnc value.csv",action="write",status="replace")
         write (19,"(20000(f0.16,',',:))") tovalue
         close(19)

     end if
     
     Do while ((tovalue(1) > 1.0d-2) .and. (evalnum < 20000))
         !!!!!!!!!!!!!!!!!!! Calculate the centroid!!!!!!!!!!!!!!!!!!!!!!!
        call bsort(dim+1,tovalue,ovalue,sorindx)
        print*, sorindx
        do j1 = 1, dim+1
            xx(:,j1)     = txx(:,sorindx(j1))
            moment(:,j1) = tmoment(:,sorindx(j1))
        end do
    
        tovalue=ovalue
        txx=xx
        tmoment=moment
    
    
        open (unit=19,file=calfileplace//"moment.csv",action="write",status="replace")   
        Do i1=1,dim+1
            write (19,"(20000(f0.16,',',:))") moment(:,i1)
        end do
        close(19)
    
        open (unit=19,file=calfileplace//"parameters.csv",action="write",status="replace")   
        Do i1=1,dim+1
            call transfrm(dim,2,txx(:,i1),temp1)
            outxx(:,i1)=temp1
            write (19,"(20000(f0.16,',',:))") outxx(:,i1)
        end do
        close(19)
    
        open (unit=19,file=calfileplace//"objective fnc value.csv",action="write",status="replace")   
            write (19,"(20000(f0.16,',',:))") ovalue
        close(19)
    
        Do i1=1, dim
            aaa=0.
            Do j1=1, dim
                aaa=aaa+xx(i1,j1)
            end do
            xbar(i1)=dble(1.0d0/dim)*aaa
        end do
         print*, 'reflection1'
        xr=xbar+rho_nm*(xbar-xx(:,dim+1))
        call transfrm(dim,2,xr,x)
        call sub_main(x,y,evalnum,exper_run)
        mr=y
        
        if (efweight==1) then
            if (parm_est==0) then
                fr = sum((tdatase**(-2.0d0))*(mr-tdatamm)**2.0d0)
            else
                fr = sum(weightvec*(tdatase**(-2.0d0))*(mr-tdatamm)**2.0d0)
            end if
            
        else
            fr = sum(log(mr/tdatamm)**2.0d0)
        end if
        
        If ((fr>=ovalue(1)).and.(fr<ovalue(dim))) then
            print*, 'reflection2'
            txx(:,dim+1) = xr
            tovalue(dim+1)=fr
            tmoment(:,dim+1)=mr
        else if (fr<ovalue(1)) then
            print*, 'expansion'
            xe=xbar+chi_nm*(xr-xbar)
            call transfrm(dim,2,xe,x)
            call sub_main(x,y,evalnum,exper_run)
            me=y
            if (efweight==1) then
                if (parm_est==0) then
                    fe = sum((tdatase**(-2.0d0))*(me-tdatamm)**2.0d0)
                else
                    fe = sum(weightvec*(tdatase**(-2.0d0))*(me-tdatamm)**2.0d0)
                end if
                
            else
                fe=sum(log(me/tdatamm)**2.0d0)
            end if
                        
            If (fe<fr) then
                print*, 'reflection with expansion'
                txx(:,dim+1)=xe
                tovalue(dim+1)=fe
                tmoment(:,dim+1)=me
            else
                print*, 'reflection without expansion'
                txx(:,dim+1)=xr
                tovalue(dim+1)=fr
                tmoment(:,dim+1)=mr
            end if
        
        else if (fr>=ovalue(dim)) then
            print*, 'contraction'
            If ((fr>=ovalue(dim)).and.(fr<ovalue(dim+1))) then
                print*, 'outside contraction'
                xc=xbar+gamma_nm*(xr-xbar)
                call transfrm(dim,2,xc,x)
                call sub_main(x,y,evalnum,exper_run)
                mc=y
                if (efweight==1) then
                    if (parm_est==0) then
                        fc = sum((tdatase**(-2.0d0))*(mc-tdatamm)**2.0d0)
                    else
                        fc = sum(weightvec*(tdatase**(-2.0d0))*(mc-tdatamm)**2.0d0)
                    end if
                    
                else
                    fc= sum(log(mc/tdatamm)**2.0d0)
                end if
                
                If (fc<=fr) then
                    print*, 'outside contraction without shrink'
                    txx(:,dim+1)=xc
                    tovalue(dim+1)=fc
                    tmoment(:,dim+1)=mc
                else
                    Do i1=2, dim+1
                        print*, 'outside contraction with shrink, i1=', i1
                        temp1=xx(:,1)+sigma_nm*(xx(:,i1)-xx(:,1))
                        call transfrm(dim,2,temp1,x)
                        call sub_main(x,y,evalnum,exper_run)
                        tempm=y
                        if (efweight==1) then
                            if (parm_est==0) then
                                fval = sum((tdatase**(-2.0d0))*(tempm-tdatamm)**2.0d0)
                            else
                                fval = sum(weightvec*(tdatase**(-2.0d0))*(tempm-tdatamm)**2.0d0)
                            end if                            
                        else
                            fval= sum(log(tempm/tdatamm)**2.0d0)
                        end if
                        
                        txx(:,i1)=temp1
                        tovalue(i1)=fval
                        tmoment(:,i1)=tempm
                        
                        open (unit=19,file=calfileplace//"contmoment.csv",action="write",status="replace")   
                        Do j1=1,i1
                            write (19,"(20000(f0.16,',',:))") tmoment(:,j1)
                        end do
                        close(19)
    
    
                        open (unit=19,file=calfileplace//"contpara.csv",action="write",status="replace")   
                        Do j1=1,i1
                            call transfrm(dim,2,txx(:,j1),temp1)
                            outxx(:,j1)=temp1
                            write (19,"(20000(f0.16,',',:))") outxx(:,j1)
                        end do
                        close(19)
    
                        open (unit=19,file=calfileplace//"contobjective fnc value.csv",action="write",status="replace")   
                            write (19,"(20000(f0.16,',',:))") tovalue(1:i1)
                        close(19)

                    end do
                end if
            else if  (fr>=ovalue(dim+1)) then
                print*, 'Inside contraction'
                xc=xbar-gamma_nm*(xbar-xx(:,dim+1))
                call transfrm(dim,2,xc,x)
                call sub_main(x,y,evalnum,exper_run)
                mc=y
                if (efweight==1) then
                    if (parm_est==0) then
                        fc = sum((tdatase**(-2.0d0))*(mc-tdatamm)**2.0d0)
                    else
                        fc = sum(weightvec*(tdatase**(-2.0d0))*(mc-tdatamm)**2.0d0)
                    end if                    
                else
                    fc = sum(log(mc/tdatamm)**2.0d0)
                end if
                
                If (fc<ovalue(dim+1)) then
                    print*, 'Inside contraction without shrink'
                    txx(:,dim+1)=xc
                    tovalue(dim+1)=fc
                    tmoment(:,dim+1)=mc
                else
                    Do i1=2, dim+1
                        print*, 'inside contraction with shrink, i1=', i1
                        temp1=xx(:,1)+sigma_nm*(xx(:,i1)-xx(:,1))
                        call transfrm(dim,2,temp1,x)
                        call sub_main(x,y,evalnum,exper_run)
                        tempm=y
                        if (efweight==1) then
                            if (parm_est==0) then
                                fval = sum((tdatase**(-2.0d0))*(tempm-tdatamm)**2.0d0)
                            else
                                fval = sum(weightvec*(tdatase**(-2.0d0))*(tempm-tdatamm)**2.0d0)
                            end if                            
                        else
                            fval = sum(log(tempm/tdatamm)**2.0d0)
                        end if
                        
                        txx(:,i1)=temp1
                        tovalue(i1)=fval
                        tmoment(:,i1)=tempm
                        
                        open (unit=19,file=calfileplace//"contmoment.csv",action="write",status="replace")   
                        Do j1=1,i1
                            write (19,"(20000(f0.16,',',:))") tmoment(:,j1)
                        end do
                        close(19)
    
                        open (unit=19,file=calfileplace//"contpara.csv",action="write",status="replace")   
                        Do j1=1,i1
                            call transfrm(dim,2,txx(:,j1),temp1)
                            outxx(:,j1)=temp1
                            write (19,"(20000(f0.16,',',:))") outxx(:,j1)
                        end do
                        close(19)
    
                        open (unit=19,file=calfileplace//"contobjective fnc value.csv",action="write",status="replace")   
                            write (19,"(20000(f0.16,',',:))") tovalue(1:i1)
                        close(19)
                        
                        
                    end do
                end if
            end if
        end if
     
     std = sqrt((sum(ovalue**2)-sum(ovalue)**2/size(ovalue))/(size(ovalue)-1)) 
     print*, std
    end do
                
                
                

        
        
        
     
     
 end if
 
    
end program prog_calibration
    