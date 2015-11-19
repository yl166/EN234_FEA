 !    Subroutines for basic 3D linear elastic elements

!==========================SUBROUTINE el_linelast_3dbasic_dynamic ==============================
subroutine el_gurson_3dbasic_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords, length_coord_array, &                               ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                       ! Input variables
    n_state_variables, initial_state_variables, &                                                        ! Input variables
    updated_state_variables,element_residual,element_deleted)                                            ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only :  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only :  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only :  dNdy => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only : dNbardy => vol_avg_shape_function_derivatives_3D
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine

    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)

    logical, intent( inout )      :: element_deleted                                        ! Set to .true. to delete element

    ! Local Variables
    integer      :: n_points,kint,i,j,a,c

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant, deter,det           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
 !   real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    real (prec)  ::  E, xnu,yy,e0,mm,q1,q2,q3,fn,en,sn,fc,ff              ! Material properties
    real (prec)  ::  F_incre(3,3),F_mid(3,3),F_mid_inv(3,3),JJ,L_incre(3,3),Lbar_incre(3,3),LL_incre(3,3)
    real (prec)  ::  eta,eta_incre,volume,e_incre(3,3),w_incre(3,3),II(3,3),R_incre(3,3),I_w_inv(3,3),temp(3,3)
    !
    !     Subroutine to compute element force vector for a linear elastodynamic problem
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0


    E = element_properties(1)
    xnu = element_properties(2)
    yy = element_properties(3)
    e0 = element_properties(4)
    mm = element_properties(5)
    q1 = element_properties(6)
    q2 = element_properties(7)
    q3 = element_properties(8)
    fn = element_properties(9)
    en = element_properties(10)
    sn = element_properties(11)
    fc = element_properties(12)
    ff = element_properties(13)

    volume = 0.d0
    eta = 0.d0
    eta_incre = 0.d0
    dNbardy=0.d0
    F_mid = 0.d0
    F_incre = 0.d0

do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        F_incre = 0.d0
        do i=1,3
            do j=1,3
                do a=1,n_nodes
                    F_incre(i,j)=F_incre(i,j)+dNdx(a,j)*dof_increment(3*(a-1)+i)
                end do
            end do
        end do
        F_mid = 0.d0
        F_mid(1,1) = 1.d0
        F_mid(2,2) = 1.d0
        F_mid(3,3) = 1.d0
        do i=1,3
            do j=1,3
                do a=1,n_nodes
                    F_mid(i,j)=F_mid(i,j)+dNdx(a,j)*(dof_total(3*(a-1)+i)+0.5d0*(dof_increment(3*(a-1)+i)))
                end do
            end do
        end do

        call invert_small(F_mid,F_mid_inv,JJ)
        !write (*,*) 1,jj
        volume = volume+w(kint)*determinant
        !JJ = JJ+JJ*w(kint)*determinant
        L_incre = matmul(F_incre,F_mid_inv)
        eta = eta + JJ*w(kint)*determinant
        eta_incre = eta_incre + JJ*(L_incre(1,1)+L_incre(2,2)+L_incre(3,3))*w(kint)*determinant
        dNdy = matmul(dNdx,F_mid_inv)
        dNbardy = dNbardy + JJ*dNdy*w(kint)*determinant
     end do
        eta = eta/volume
        eta_incre = eta_incre/(volume*eta)
        dNbardy = dNbardy/(volume*eta)

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        F_incre = 0.d0
         do i=1,3
             do j=1,3
                 do a=1,n_nodes
                     F_incre(i,j)=F_incre(i,j)+dNdx(a,j)*dof_increment(3*(a-1)+i)
                 end do
             end do
         end do
         !write (*,*) 1,F_incre
        F_mid = 0.d0
        F_mid(1,1) = 1.d0
        F_mid(2,2) = 1.d0
        F_mid(3,3) = 1.d0
        do i=1,3
            do j=1,3
                do a=1,n_nodes
                    F_mid(i,j)=F_mid(i,j)+dNdx(a,j)*(dof_total(3*(a-1)+i)+0.5d0*(dof_increment(3*(a-1)+i)))
                end do
            end do
        end do
        call invert_small(F_mid,F_mid_inv,JJ)
        !write(*,*) 2,F_mid_inv
        Lbar_incre=0.d0
        Lbar_incre = matmul(F_incre,F_mid_inv)
        !write (*,*) 'Lbar_incre',Lbar_incre
        do i=1,3
            do j=1,3
                if (i==j) then
                    Lbar_incre(i,j) = Lbar_incre(i,j)+(eta_incre-&
                    (Lbar_incre(1,1)+Lbar_incre(2,2)+Lbar_incre(3,3)))/3.d0
                end if
            end do
        end do

!
        e_incre= (Lbar_incre+transpose(Lbar_incre))/2.d0
        w_incre= (Lbar_incre-transpose(Lbar_incre))/2.d0
        II=0.d0
        II(1,1)=1.d0
        II(2,2)=1.d0
        II(3,3)=1.d0
        temp=II-w_incre/2.d0
       ! write (*,*) 'I=',II
        !write (*,*) 'temp=',temp
        call invert_small(temp,I_w_inv,det)
        !write (*,*) 2, det
        R_incre = matmul(I_w_inv,(II+w_incre/2.d0))

        call gurson(element_properties,n_properties,8,initial_state_variables&
        (((kint-1)*8+1):((kint-1)*8+8)),updated_state_variables&
        (((kint-1)*8+1):((kint-1)*8+8)),e_incre,R_incre,stress)




        if(updated_state_variables((kint-1)*8+7)>ff) then
            element_deleted = .true.
        else
            element_deleted = .false.
        end if

        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)+1.d0/3.d0*(dNbardy(1:n_nodes,1)-dNdy(1:n_nodes,1))
        B(1,2:3*n_nodes-1:3) = 1.d0/3.d0*(dNbardy(1:n_nodes,2)-dNdy(1:n_nodes,2))
        B(1,3:3*n_nodes:3)   = 1.d0/3.d0*(dNbardy(1:n_nodes,3)-dNdy(1:n_nodes,3))
        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)+1.d0/3.d0*(dNbardy(1:n_nodes,2)-dNdy(1:n_nodes,2))
        B(2,3:3*n_nodes:3)   = 1.d0/3.d0*(dNbardy(1:n_nodes,3)-dNdy(1:n_nodes,3))
        B(2,1:3*n_nodes-2:3) = 1.d0/3.d0*(dNbardy(1:n_nodes,1)-dNdy(1:n_nodes,1))
        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)+1.d0/3.d0*(dNbardy(1:n_nodes,3)-dNdy(1:n_nodes,3))
        B(3,2:3*n_nodes-1:3) = 1.d0/3.d0*(dNbardy(1:n_nodes,2)-dNdy(1:n_nodes,2))
        B(3,1:3*n_nodes-2:3) = 1.d0/3.d0*(dNbardy(1:n_nodes,1)-dNdy(1:n_nodes,1))
        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)

        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

    end do

    return
end subroutine el_gurson_3dbasic_dynamic

!
!==========================SUBROUTINE el_linelast_3dbasic ==============================

subroutine gurson(element_properties,n_properties,n_state_variables,initial_state_variables, &
updated_state_variables,dstrain,dRot,stress1)
     use Types
     use ParamIO
     use Globals, only : TIME, DTIME
     use Element_Utilities, only : invert_small

     implicit none

     integer, intent( in ) :: n_properties
     integer, intent( in ) :: n_state_variables

     real (prec), intent( in ) :: element_properties(n_properties)
     real (prec), intent( in ) :: initial_state_variables(n_state_variables)
     real (prec), intent( in ) :: dstrain(3,3)
     real (prec), intent( in ) :: dRot(3,3)
     real (prec), intent( out ) :: stress1(6)
     real (prec), intent( out ) :: updated_state_variables(n_state_variables)
    integer :: i,j,k,l
    real (prec)  :: stress0(6), ematrix, Vf,phisq,eps,pp,sigma,yy
    real (prec)  :: taomat(3,3), taomatdiv(3,3), dstraindiv(3,3), pstar, Sstar(3,3), phi, sigmastar, fstar
    real (prec)  :: ffbar, tol, taomat1(3,3),e,xnu,e0dot,mm,q1,q2,q3,fn,een,ssn,fc,ff,dee,dev
    real (prec)  :: ddeevv(2),dphidse,dphidp,dphiddee,dphiddev,dtempddee,dtempddev,se,matrix(2,2)
    real (prec)  :: matrix_inv(2,2),df1ddee,df1ddev,df2ddee,df2ddev,dpddev,dseddee,dt,f1,f2,line(2)
    real (prec)  :: matrix_det,p,temp,error,vf1,ematrix1,dematrix

    E = element_properties(1)
    xnu = element_properties(2)
    YY=element_properties(3)
    e0dot=element_properties(4)
    mm=element_properties(5)
    q1=element_properties(6)
    q2=element_properties(7)
    q3=element_properties(8)
    fn=element_properties(9)
    een=element_properties(10)
    ssn=element_properties(11)
    fc=element_properties(12)
    ff=element_properties(13)

    stress0 = initial_state_variables(1:6) ! Stress at start of increment
    ematrix = initial_state_variables(7)
    Vf = initial_state_variables(8)

!
    taomat = 0.d0
    taomat(1,1) = stress0(1)
    taomat(2,2) = stress0(2)
    taomat(3,3) = stress0(3)
    taomat(1,2) = stress0(4)
    taomat(2,1) = stress0(4)
    taomat(1,3) = stress0(5)
    taomat(3,1) = stress0(5)
    taomat(2,3) = stress0(6)
    taomat(3,2) = stress0(6)

!
    taomatdiv = 0.d0
    do i = 1,3
        do j = 1,3
        if(i==j)then
        taomatdiv(i,j) = taomat(i,j)-(taomat(1,1)+taomat(2,2)+taomat(3,3))/3.d0
        else
        taomatdiv(i,j)=taomat(i,j)
        end if
        end do
    end do

    dstraindiv= 0.d0
    do i = 1,3
        do j = 1,3
            if (i == j) then
            dstraindiv(i,j) = dstrain(i,j)-(dstrain(1,1)+dstrain(2,2)+dstrain(3,3))/3.d0
            else
            dstraindiv(i,j) = dstrain(i,j)
            end if
        end do
    end do

    pstar = (taomat(1,1)+taomat(2,2)+taomat(3,3))/3.d0+E/(3.d0*(1.d0-2.d0*xnu))* &
                     (dstrain(1,1)+ dstrain(2,2)+ dstrain(3,3))

   Sstar = E/(1.d0+xnu)*dstraindiv + matmul(dRot,matmul(taomatdiv,transpose(dRot)))
   !write(*,*)'drot',drot
  !  write(*,*)'sstar',sstar

!
    ffbar = (q1+dsqrt(q1**2.d0-q3))/q3

    if (Vf < fc) then
    fstar = Vf
    else if (fc < Vf .AND. Vf < ff) then
    fstar = fc + (ffbar-fc)/(ff-fc)*(Vf-fc)
    else if (Vf > ff) then
    fstar = ffbar
    end if

!

    sigmastar = dsqrt(3.d0/2.d0*(Sstar(1,1)**2.d0+Sstar(2,2)**2.d0+Sstar(3,3)**2.d0+ &
                       2.d0*Sstar(1,2)**2.d0+2.d0*Sstar(2,3)**2.d0+2.d0*Sstar(1,3)**2.d0))
  !write(*,*) 's1',sigmastar
    Phisq = sigmastar**2.d0/YY**2.d0+2.d0*q1*fstar*cosh(3.d0/2.d0*q2*pstar/YY)-(1.d0+q3*fstar**2.d0)

!
   eps = 10.d0**(-8.d0)
   taomat1=0.d0
    if (Phisq < eps) then
        do i = 1,3
            do j = 1,3
            if (i == j) then
            taomat1(i,j) = Sstar(i,j)+pstar
            else
            taomat1(i,j) = Sstar(i,j)
            end if
            end do
         end do
       dematrix=0.d0
       ematrix1=ematrix
    else

!

        dee = 0.d0
        dev = 0.d0
        ddeevv=1.d0
        error=1.d0
        do while(error>eps)
        dt=dtime
        !write(*,*)dt
        se = sigmastar - 1.5d0*E/(1.d0+xnu)*dee
   !     write(*,*)'s2',sigmastar
        p = pstar - E/3.d0/(1.d0-2.d0*xnu)*dev
        phi = se**2.d0/YY**2.d0 + 2.d0*q1*fstar*cosh(1.5d0*q2*p/Yy) &
             - (1.d0+q3*fstar**2.d0)
        if (phi<0.d0) then
                dee=dee/10.d0
                dev=dev/10.d0
                cycle
            endif
        phi=dsqrt(phi)
        dphidse = se / phi / Yy**2.d0
        dphidp = 1.5d0 * q1 * q2 * fstar / phi / Yy * sinh(1.5d0 * q2 * p / Yy)
        dseddee = -1.5d0*e/(1+xnu)
        dpddev= -e/(3.d0*(1-2.d0*xnu))
        dphiddee = dphidse * dseddee
        dphiddev = dphidp * dpddev
        temp = dsqrt(dphidse**2.d0 + 2.d0/9.d0*dphidp**2.d0)
        dtempddee = (2.d0*se*dseddee/phi/phi/Yy**4.d0 - se**2.d0/Yy**4.d0*2.d0/phi**3.d0*dphiddee &
                    + 2.d0/9.d0*dphidp**2.d0*(-2.d0)/phi*dphiddee)/2.d0/temp
        dtempddev = (q1**2.d0*q2**2.d0*fstar**2.d0/phi**2.d0/Yy**2.d0*sinh(1.5*q2/Yy*p)*cosh(1.5*q2/Yy*p) &
                    * 1.5*q2/Yy*dpddev + 2.d0/9.d0*dphidp**2.d0*(-2.d0)/phi*dphiddev-se**2.d0/yy**4.d0*2.d0&
                    /phi**3.d0*dphiddev)/2.d0/temp
        F1 = temp * dee /dt/e0dot - dphidse * phi**mm
        F2 = temp * dev /dt/e0dot - dphidp * phi**mm

        dF1ddee = temp/dt/e0dot - dseddee / Yy**2.d0 * phi**(mm-1.d0) - &
                    se/Yy**2.d0*(mm-1.d0)*phi**(mm-2.d0)*dphiddee + dee/dt/e0dot*dtempddee
        dF1ddev = - se/Yy**2.d0*(mm-1.d0)*phi**(mm-2.d0)*dphiddev + dee/dt/e0dot*dtempddev
        dF2ddee = dev/dt/e0dot * dtempddee - dphidp*(phi**mm)*(mm-1.d0)/(phi)*dphiddee
        dF2ddev = temp/dt/e0dot - dphidp*(phi**mm)*(mm-1.d0)/(phi)*dphiddev - 1.5*q1*q2*fstar/Yy*phi**(mm-1.d0)* &
                    cosh(1.5d0*q2*p/Yy) * dpddev*1.5d0*q2/Yy + dev/dt/e0dot * dtempddev

            matrix(1,1:2)=[dF1ddee,df1ddev]
            matrix(2,1:2)=[df2ddee,df2ddev]
            call invert_small(matrix,matrix_inv,matrix_det)
            line=[-f1,-f2]
            ddeevv=matmul(matrix_inv,line)
            !write(*,*)matrix_inv
           ! write(*,*) ddeevv
            error=dsqrt(ddeevv(1)**2.d0+ddeevv(2)**2.d0)
            dee=dee+ddeevv(1)
            dev=dev+ddeevv(2)


        end do

     taomat1=0.d0
        do i=1,3
            do j=1,3
                if(i==j) then
                    if(sigmastar/=0.d0) then
                    taomat1(i,j)=sstar(i,j)-dee*e/(1.d0+xnu)*1.5d0*sstar(i,j)/sigmastar+(pstar-e/(3.d0*(1.d0-2.d0*xnu))*dev)
                    else
                     taomat1(i,j)=sstar(i,j)+(pstar-e/(3.d0*(1.d0-2.d0*xnu))*dev)
                     end if
                else
                     if(sigmastar/=0.d0) then
                    taomat1(i,j)=sstar(i,j)-dee*e/(1.d0+xnu)*1.5d0*sstar(i,j)/sigmastar
                    else
                     taomat1(i,j)=sstar(i,j)
                     end if
                end if
            end do
        end do

        dematrix=e0dot*dtime/(1.d0-vf)*phi**mm*(dphidse**2.d0+2.d0/9.d0*dphidp**2.d0)**(-1.d0/2.d0)&
        *(dphidse*se+1.d0/3.d0*dphidp*p)

        ematrix1=ematrix+dematrix

end if

vf1=1.d0+(vf-1)*exp(-dev)+fn*dematrix/(ssn*dsqrt(2.d0*3.1415d0))*exp(-0.5d0*((ematrix-een)/ssn)**2.d0)

stress1(1)=taomat1(1,1)
stress1(2)=taomat1(2,2)
stress1(3)=taomat1(3,3)
stress1(4)=taomat1(1,2)
stress1(5)=taomat1(1,3)
stress1(6)=taomat1(2,3)
 !stressDn=stress0
updated_state_variables(1:6)=stress1(1:6)
updated_state_variables(7)=ematrix1
updated_state_variables(8)=vf1

 end subroutine gurson

!==========================SUBROUTINE fieldvars_linelast_3dbasic ==============================
subroutine fieldvars_gurson_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only: dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only: dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only: dNdy => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only : dNbardy => vol_avg_shape_function_derivatives_3D
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
    integer, intent( in )         :: n_field_variables                                      ! # field variables

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    character (len=100), intent(in) :: field_variable_names(n_field_variables)

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( inout )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step

    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables

    ! Local Variables
    logical      :: strcmp

    integer      :: n_points,kint,k,i,j,a,c

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
   ! real (prec)  :: E, xnu, D44, D11, D12,volume              ! Material properties
    real (prec)  :: p, smises                          ! Pressure and Mises stress
     real (prec)  ::  E, xnu,yy,e0,mm,q1,q2,q3,fn,en,sn,fc,ff              ! Material properties
    real (prec)  ::  F_incre(3,3),F_mid(3,3),F_mid_inv(3,3),JJ,L_incre(3,3),Lbar_incre(3,3),LL_incre(3,3)
    real (prec)  ::  eta,eta_incre,volume,e_incre(3,3),w_incre(3,3),II(3,3),R_incre(3,3),I_w_inv(3,3)
    real (prec)  :: deter,det,temp(3,3),Vff
    !
    !     Subroutine to compute element contribution to project element integration point data to nodes

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    nodal_fieldvariables = 0.d0


    E = element_properties(1)
    xnu = element_properties(2)
    yy = element_properties(3)
    e0 = element_properties(4)
    mm = element_properties(5)
    q1 = element_properties(6)
    q2 = element_properties(7)
    q3 = element_properties(8)
    fn = element_properties(9)
    en = element_properties(10)
    sn = element_properties(11)
    fc = element_properties(12)
    ff = element_properties(13)


    volume = 0.d0
    eta = 0.d0
    eta_incre = 0.d0
    dNbardy=0.d0
    F_mid = 0.d0
    F_incre = 0.d0

do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        F_incre = 0.d0
        do i=1,3
            do j=1,3
                do a=1,n_nodes
                    F_incre(i,j)=F_incre(i,j)+dNdx(a,j)*dof_increment(3*(a-1)+i)
                end do
            end do
        end do
        F_mid = 0.d0
        F_mid(1,1) = 1.d0
        F_mid(2,2) = 1.d0
        F_mid(3,3) = 1.d0
        do i=1,3
            do j=1,3
                do a=1,n_nodes
                    F_mid(i,j)=F_mid(i,j)+dNdx(a,j)*(dof_total(3*(a-1)+i)+0.5d0*(dof_increment(3*(a-1)+i)))
                end do
            end do
        end do

        call invert_small(F_mid,F_mid_inv,JJ)
        !write (*,*) 1,jj
        volume = volume+w(kint)*determinant
        !JJ = JJ+JJ*w(kint)*determinant
        L_incre = matmul(F_incre,F_mid_inv)
        eta = eta + JJ*w(kint)*determinant
        eta_incre = eta_incre + JJ*(L_incre(1,1)+L_incre(2,2)+L_incre(3,3))*w(kint)*determinant
        dNdy = matmul(dNdx,F_mid_inv)
        dNbardy = dNbardy + JJ*dNdy*w(kint)*determinant
     end do
        eta = eta/volume
        eta_incre = eta_incre/(volume*eta)
        dNbardy = dNbardy/(volume*eta)

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        F_incre = 0.d0
         do i=1,3
             do j=1,3
                 do a=1,n_nodes
                     F_incre(i,j)=F_incre(i,j)+dNdx(a,j)*dof_increment(3*(a-1)+i)
                 end do
             end do
         end do
         !write (*,*) 1,F_incre
        F_mid = 0.d0
        F_mid(1,1) = 1.d0
        F_mid(2,2) = 1.d0
        F_mid(3,3) = 1.d0
        do i=1,3
            do j=1,3
                do a=1,n_nodes
                    F_mid(i,j)=F_mid(i,j)+dNdx(a,j)*(dof_total(3*(a-1)+i)+0.5d0*(dof_increment(3*(a-1)+i)))
                end do
            end do
        end do
        call invert_small(F_mid,F_mid_inv,JJ)
        !write(*,*) 2,F_mid_inv
        Lbar_incre=0.d0
        Lbar_incre = matmul(F_incre,F_mid_inv)
        !write (*,*) 'Lbar_incre',Lbar_incre
        do i=1,3
            do j=1,3
                if (i==j) then
                    Lbar_incre(i,j) = Lbar_incre(i,j)+(eta_incre-&
                    (Lbar_incre(1,1)+Lbar_incre(2,2)+Lbar_incre(3,3)))/3.d0
                end if
            end do
        end do

!
        e_incre= (Lbar_incre+transpose(Lbar_incre))/2.d0
        w_incre= (Lbar_incre-transpose(Lbar_incre))/2.d0
        II=0.d0
        II(1,1)=1.d0
        II(2,2)=1.d0
        II(3,3)=1.d0
        temp=II-w_incre/2.d0
       ! write (*,*) 'I=',II
        !write (*,*) 'temp=',temp
        call invert_small(temp,I_w_inv,det)
        !write (*,*) 2, det
        R_incre = matmul(I_w_inv,(II+w_incre/2.d0))

        call gurson(element_properties,n_properties,8,initial_state_variables&
        (((kint-1)*8+1):((kint-1)*8+8)),updated_state_variables&
        (((kint-1)*8+1):((kint-1)*8+8)),e_incre,R_incre,stress)

        vff=updated_state_variables(((kint-1)*8+8))


        p = sum(stress(1:3))/3.d0
        sdev = stress
        sdev(1:3) = sdev(1:3)-p
        smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
        ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
        do k = 1,n_field_variables
            if (strcmp(field_variable_names(k),'S11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(5)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(6)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'SMISES',6) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + smises*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'Vf',2) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + Vff*N(1:n_nodes)*determinant*w(kint)
            endif
        end do

    end do

    return
end subroutine fieldvars_gurson_3dbasic

!
!
!
!
!
!
