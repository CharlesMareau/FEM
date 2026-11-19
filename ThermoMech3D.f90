program fem_main
 use utils
 use tensor
 use iso_fortran_env
 implicit none
 integer :: freq
! Variables maillage
 integer :: nnode,nnset,nelset,nelem3d,nelem2D,ndof,nnode_per_element3d,nnode_per_element2D,nntg,nntg_per_element3d,nntg_per_element2D
 real(8), allocatable :: coords(:,:),location(:,:)
 integer, allocatable :: connect3d(:,:),connect2d(:,:),idelem2d(:),idelem3d(:)
! Données matériaux
 character(5) :: cname
 character(len=64), allocatable :: nset_names(:),elset_names(:)
 integer, allocatable :: nset_nodes(:),nset_start(:),nset_end(:)
 integer, allocatable :: elset_elem(:),elset_start(:),elset_end(:)
! Variables pour assemblage
 integer :: i,j,k,l,ij,kl,idof,jdof,kdof,ldof,iddl,jddl,kddl,ielem,inode,jnode,knode,iter,istp,unit,istat,intg,jntg,igss,jgss,iset,nstat,nprop
 integer, allocatable :: bc_mech_type(:),bc_therm_type(:)
 real(8), allocatable :: bc_mech_val(:,:),bc_therm_val(:,:)
 real(8), allocatable :: interpolation(:,:),gradient(:,:,:),volume(:),wght(:)
 real(8), allocatable :: stiffness(:),residual_mech(:),increment_mech(:)
 real(8), allocatable :: conductance(:),residual_therm(:),increment_therm(:)
 real(8), allocatable :: position_new(:),position_old(:)
 real(8), allocatable :: temperature_new(:),temperature_old(:)
 real(8), allocatable :: temp_new(:),temp_old(:)
 real(8), allocatable :: capacity_new(:)
 real(8), allocatable :: source_new(:)
 type(Tensor2), allocatable :: conductivity_new(:)
 type(Tensor1), allocatable :: heatflow_new(:),heatflow_old(:)
 type(Tensor1), allocatable :: tempgrad_new(:),tempgrad_old(:)
 type(Tensor2), allocatable :: defgrad_new(:),defgrad_old(:)
 type(Tensor2), allocatable :: stress_new(:),stress_old(:)
 type(Tensor4), allocatable :: modulus_new(:)
 real(8), allocatable :: state_new(:,:),state_old(:,:)
 real(8), allocatable :: internalenergy_new(:),internalenergy_old(:)
 real(8), allocatable :: dissipatedenergy_new(:),dissipatedenergy_old(:)
 real(8), allocatable :: x(:),y(:),z(:),matrix(:,:),vector(:),auxiliary(:,:),props(:)
 real(8) :: dtime_min,dtime_max
 real(8) :: time_old,time_new,dtime
 real(8) :: duration,progress,delta
 integer :: nnz,mmz,nrhs,phase,msglvl
 integer, allocatable :: ia_mech(:),ja_mech(:),correspond_mech(:,:),permutation_mech(:)
 integer :: mtype_mech,iparm_mech(64),pt_mech(64)
 integer :: maxfct_mech,mnum_mech,error_mech
 integer, allocatable :: ia_therm(:),ja_therm(:),correspond_therm(:,:),permutation_therm(:)
 integer :: mtype_therm,iparm_therm(64),pt_therm(64)
 integer :: maxfct_therm,mnum_therm,error_therm
! *** variables for data output
 character(40) :: filevtu
 character(8) :: fil
 integer :: ifile
 real(8), parameter :: small=1.d-5
 real(8) :: res_mech_old,res_mech,damp_mech,toler_mech=1.d-2
 real(8) :: res_therm_old,res_therm,damp_therm,toler_therm=1.d-2
 integer :: idof_max
 character(80) :: filemsh,fileproc,filemat,prosa
 integer, parameter :: iter_max=10000

! *** read main input file
 open(unit=1,file='implicit3D.in',status='old')
  read(1,'(a)') prosa
  read(1,'(a)') filemsh
  read(1,'(a)') prosa
  read(1,'(a)') fileproc
  read(1,'(a)') prosa
  read(1,'(a)') filemat
  read(1,'(a)') prosa
  read(1,*) freq
 close(unit=1)

! *** read material input file
 call read_material(filemat,cname,nstat,nprop,props)

! *** read mesh file
 call read_inp(filemsh)
!  call read_msh("hole.msh", coords, conn, nnode, nelem3d)
 write(*,'("Number of 2D elements",i6)') nelem2d
 write(*,'("Number of 3D elements",i6)') nelem3d
 write(*,'("Number of nodes",i6)') nnode
 ndof=3*nnode
 write(*,'("Number of DOFs",i6)') ndof
 nntg=nntg_per_element3d*nelem3d
 write(*,'("Number of integration points",i6)') nntg
 write(*,'("Number of nodes per 2D element",i6)') nnode_per_element2d
 write(*,'("Number of nodes per 3D element",i6)') nnode_per_element3d

 allocate(matrix(3*nnode_per_element3d,3*nnode_per_element3d),vector(3*nnode_per_element3d),auxiliary(9,3*nnode_per_element3d))
! *** initializes variables
 allocate(position_new(ndof),position_old(ndof))
 allocate(temperature_new(nnode),temperature_old(nnode))
 allocate(temp_new(nntg),temp_old(nntg))
 allocate(tempgrad_new(nntg),tempgrad_old(nntg))
 allocate(heatflow_new(nntg),heatflow_old(nntg))
 allocate(defgrad_new(nntg),defgrad_old(nntg))
 allocate(stress_new(nntg),stress_old(nntg))
 allocate(modulus_new(nntg))
 allocate(source_new(nntg))
 allocate(capacity_new(nntg))
 allocate(conductivity_new(nntg))
 allocate(state_new(nstat,nntg),state_old(nstat,nntg))
 allocate(dissipatedenergy_new(nntg),dissipatedenergy_old(nntg))
 allocate(internalenergy_new(nntg),internalenergy_old(nntg))
 do ielem=1,nelem3d
  do inode=1,nnode_per_element3d
   idof=3*connect3d(inode,ielem)-2
   jdof=3*connect3d(inode,ielem)-1
   kdof=3*connect3d(inode,ielem)
   position_old(idof)=coords(1,connect3d(inode,ielem))
   position_old(jdof)=coords(2,connect3d(inode,ielem))
   position_old(kdof)=coords(3,connect3d(inode,ielem))
   position_new(idof)=position_old(idof)
   position_new(jdof)=position_old(jdof)
   position_new(kdof)=position_old(kdof)
   ldof=connect3d(inode,ielem)
   temperature_old(ldof)=0.d0
   temperature_new(ldof)=temperature_old(ldof)
  enddo
  do intg=1,nntg_per_element3d
   igss=(ielem-1)*nntg_per_element3d+intg
   internalenergy_old(igss)=0.d0
   internalenergy_new(igss)=0.d0
   dissipatedenergy_old(igss)=0.d0
   dissipatedenergy_new(igss)=0.d0
   stress_old(igss)=null2
   stress_new(igss)=null2
   defgrad_old(igss)=id2
   defgrad_new(igss)=id2
   temp_old(igss)=0.d0
   temp_new(igss)=0.d0
   tempgrad_old(igss)=null1
   tempgrad_new(igss)=null1
   heatflow_old(igss)=null1
   heatflow_new(igss)=null1
   do istat=1,nstat
    state_old(istat,igss)=0.d0
    state_new(istat,igss)=0.d0
   enddo
  enddo
 enddo
! *** calculates the shape function derivative
 allocate(interpolation(nnode_per_element3d,nntg),gradient(3,nnode_per_element3d,nntg),location(3,nntg),wght(nntg),volume(nelem3d))
 allocate(x(nnode_per_element3d),y(nnode_per_element3d),z(nnode_per_element3d))
 do ielem=1,nelem3d
  do inode=1,nnode_per_element3d
   x(inode)=coords(1,connect3d(inode,ielem))
   y(inode)=coords(2,connect3d(inode,ielem))
   z(inode)=coords(3,connect3d(inode,ielem))
  enddo
  igss=(ielem-1)*nntg_per_element3d+1
  jgss=(ielem-1)*nntg_per_element3d+nntg_per_element3d
  select case(nntg_per_element3d)
   case(1)
    call tetrahedronT4(x,y,z,interpolation(:,igss:jgss),gradient(:,:,igss:jgss),location(:,igss:jgss),wght(igss:jgss),volume(ielem))
   case(4)
    call tetrahedronT10(x,y,z,interpolation(:,igss:jgss),gradient(:,:,igss:jgss),location(:,igss:jgss),wght(igss:jgss),volume(ielem))
   case(8)
    call hexahedronH8(x,y,z,interpolation(:,igss:jgss),gradient(:,:,igss:jgss),location(:,igss:jgss),wght(igss:jgss),volume(ielem))
   case default
    write(*,'("Unknown type of element. Program will stop !")')
    stop
  end select
 enddo
 deallocate(x,y,z)
! *** generates the csr structure for storage
 call build_csr_structure
! *** sets the parameters for the PARDISO solver
 nrhs  = 1
 msglvl = 0
 mtype_mech = 11                ! matrice réelle symétrique indéfinie : 11 / non symétrique : 11
 error_mech = 0
 maxfct_mech = 1
 mnum_mech = 1
 pt_mech = 0
 iparm_mech = 0
 iparm_mech(1) = 1              ! utiliser iparm personnalisé
 iparm_mech(3) = 0              ! pas de message
 iparm_mech(7) = 1              ! nombre de threads
 iparm_mech(10) = 13            ! pivot perturbation
 iparm_mech(17) = -1            ! pas de réordonnancement METIS
 iparm_mech(34) = 0             ! utiliser zero-based indexing → =1 pour C, =0 pour Fortran
 mtype_therm = 11                ! matrice réelle symétrique indéfinie : 11 / non symétrique : 11
 error_therm = 0
 maxfct_therm = 1
 mnum_therm = 1
 pt_therm = 0
 iparm_therm = 0
 iparm_therm(1) = 1              ! utiliser iparm personnalisé
 iparm_therm(3) = 0              ! pas de message
 iparm_therm(7) = 1              ! nombre de threads
 iparm_therm(10) = 13            ! pivot perturbation
 iparm_therm(17) = -1            ! pas de réordonnancement METIS
 iparm_therm(34) = 0             ! utiliser zero-based indexing → =1 pour C, =0 pour Fortran
! *** finds nodes with prescribed BC
 allocate(bc_mech_type(ndof),bc_mech_val(2,ndof),bc_therm_type(nnode),bc_therm_val(2,nnode))
 call read_step(fileproc)
! *** saves initial state
 if (freq.ne.0) then
  ifile=0
  write(fil,'(i8)') ifile
  fil=adjustl(fil)
  filevtu="result"//trim(fil)//".vtu"
  call export_vtu(filevtu)
  call export_rct(0)
  call export_rct(1)
 endif
 time_old=0.d0
 istp=0
 progress=0.d0
 do while (time_old.lt.duration)
  istp=istp+1
!   dtime=min(duration-time_old,dtime)
  time_new=time_old+dtime
  if (time_new.gt.duration) then
   dtime=duration-time_old
   time_new=duration
  endif
  if ((time_new+dtime_min).gt.duration) then
   dtime=duration-time_old
   time_new=duration
  endif
  write(*,'("ISTEP:",i6," TIME:",2es12.4)') istp,time_old,dtime
  progress=progress+dtime/duration
! *** applies boundary condition
  !$omp parallel default(shared)
  !$omp do private(idof,i,inode) schedule(static)
  do idof=1,ndof
   if (bc_mech_type(idof).eq.1) then
    inode=(idof+2)/3
    i=idof-3*inode+3
    position_new(idof)=coords(i,inode)+bc_mech_val(1,idof)+(bc_mech_val(2,idof)-bc_mech_val(1,idof))*progress
   endif
  enddo
  !$omp enddo
  !$omp do private(inode) schedule(static)
  do inode=1,nnode
   if (bc_therm_type(inode).eq.1) then
    temperature_new(inode)=bc_therm_val(1,inode)+(bc_therm_val(2,inode)-bc_therm_val(1,inode))*progress
   endif
  enddo
  !$omp enddo
  !$omp end parallel
  res_mech=toler_mech+1.d0
  damp_mech=1.d0
  res_therm=toler_therm+1.d0
  damp_therm=1.d0
  iter=0
! *** solves mechanical problem
  do while ((res_mech.gt.toler_mech).or.(res_therm.gt.toler_therm))
   iter=iter+1
   if (iter.gt.iter_max) then
    write(*,'("MAX NUMBER OF ITERATIONS HAS BEEN REACHED.")')
    call export_vtu("fail.vtu")
    stop
   endif
! *** calculates the deformation gradient tensor
   !$omp parallel default(shared)
   !$omp do private(ielem,intg,igss,inode,idof,i,j) schedule(static)
   do ielem=1,nelem3d
    do intg=1,nntg_per_element3d
     igss=(ielem-1)*nntg_per_element3d+intg
     defgrad_new(igss)=null2
     do inode=1,nnode_per_element3d
      do i=1,3
       idof=3*(connect3d(inode,ielem)-1)+i
       do j=1,3
        defgrad_new(igss)%val(i,j)=defgrad_new(igss)%val(i,j)+position_new(idof)*gradient(j,inode,igss)
       enddo
      enddo
     enddo
    enddo
   enddo
   !$omp enddo
! *** calculates the temperature and its gradient
   !$omp do private(ielem,intg,igss,inode,idof,i) schedule(static)
   do ielem=1,nelem3d
    do intg=1,nntg_per_element3d
     igss=(ielem-1)*nntg_per_element3d+intg
     temp_new(igss)=0.d0
     tempgrad_new(igss)=null1
     do inode=1,nnode_per_element3d
      idof=connect3d(inode,ielem)
      temp_new(igss)=temp_new(igss)+temperature_new(idof)*interpolation(inode,igss)
      do i=1,3
       tempgrad_new(igss)%val(i)=tempgrad_new(igss)%val(i)+temperature_new(idof)*gradient(i,inode,igss)
      enddo
     enddo
    enddo
   enddo
   !$omp enddo
   !$omp end parallel
! *** integrates constitutive relations
   select case(cname)
    case("ELAST","elast","Elast")
     !$omp parallel default(shared)
     !$omp do private(intg) schedule(static)
     do intg=1,nntg
      call thermoelasticity(nprop,nstat,props,time_old,time_new,defgrad_old(intg),defgrad_new(intg),stress_old(intg),stress_new(intg),state_old(:,intg),state_new(:,intg) &
                           ,temp_old(intg),temp_new(intg),tempgrad_old(intg),tempgrad_new(intg),heatflow_old(intg),heatflow_new(intg) &
                           ,internalenergy_old(intg),internalenergy_new(intg),dissipatedenergy_old(intg),dissipatedenergy_new(intg) &
                           ,modulus_new(intg),capacity_new(intg),source_new(intg),conductivity_new(intg))
     enddo
     !$omp enddo
     !$omp end parallel
    case("PLAST","plast","Plast")
     !$omp parallel default(shared)
     !$omp do private(intg) schedule(static)
     do intg=1,nntg
      call thermoplasticity(nprop,nstat,props,time_old,time_new,defgrad_old(intg),defgrad_new(intg),stress_old(intg),stress_new(intg),state_old(:,intg),state_new(:,intg) &
                           ,temp_old(intg),temp_new(intg),tempgrad_old(intg),tempgrad_new(intg),heatflow_old(intg),heatflow_new(intg) &
                           ,internalenergy_old(intg),internalenergy_new(intg),dissipatedenergy_old(intg),dissipatedenergy_new(intg) &
                           ,modulus_new(intg),capacity_new(intg),source_new(intg),conductivity_new(intg))
     enddo
     !$omp enddo
     !$omp end parallel
    case default
     write(*,'("ERROR: UNKNOWN MATERIAL TYPE")')
     stop
   end select
! ! ! !    !$omp parallel default(shared)
! ! ! !    !$omp do private(intg) schedule(static)
! ! ! !    do intg=1,nntg
! ! ! !     call thermal(nprop,nstat,props,time_old,time_new,temp_old(intg),temp_new(intg),tempgrad_old(intg),tempgrad_new(intg),state_old(:,intg),state_new(:,intg) &
! ! ! !                    ,heatflow_old(intg),heatflow_new(intg),capacity_new(intg),conductivity_new(intg))
! ! ! !    enddo
! ! ! !    !$omp enddo
! ! ! !    !$omp end parallel
   stiffness=0.d0
   residual_mech=0.d0
   conductance=0.d0
   residual_therm=0.d0
! *** evaluates stiffness matrix
   !$omp parallel default(shared)
   !$omp do schedule(static) private(ielem,intg,igss,matrix,inode,iddl,idof,knode,kddl,kdof,jdof,i,j,k,l)
   do ielem=1,nelem3d
    do intg=1,nntg_per_element3d
     igss=(ielem-1)*nntg_per_element3d+intg
     do inode=1,nnode_per_element3d
      do i=1,3
       idof=3*(connect3d(inode,ielem)-1)+i
       iddl=3*(inode-1)+i
       do knode=1,nnode_per_element3d
        do k=1,3
         kdof=3*(connect3d(knode,ielem)-1)+k
         kddl=3*(knode-1)+k
         matrix(iddl,kddl)=0.d0
         do j=1,3
          do l=1,3
           matrix(iddl,kddl)=matrix(iddl,kddl)+gradient(j,inode,igss)*modulus_new(igss)%val(i,j,k,l)*gradient(l,knode,igss)
          enddo
         enddo
         matrix(iddl,kddl)=matrix(iddl,kddl)*wght(igss)
         jdof=correspond_mech(idof,kdof)
         !$omp atomic update
         stiffness(jdof)=stiffness(jdof)+matrix(iddl,kddl)
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo
   !$omp end do
! *** evaluates conductance matrix
   !$omp do schedule(static) private(ielem,intg,igss,matrix,inode,idof,knode,kdof,jdof,i,k)
   do ielem=1,nelem3d
    do intg=1,nntg_per_element3d
     igss=(ielem-1)*nntg_per_element3d+intg
     do inode=1,nnode_per_element3d
      idof=connect3d(inode,ielem)
      do knode=1,nnode_per_element3d
       kdof=connect3d(knode,ielem)
       matrix(inode,knode)=interpolation(inode,igss)*capacity_new(igss)*interpolation(knode,igss)/(time_new-time_old)
       do i=1,3
        do k=1,3
         matrix(inode,knode)=matrix(inode,knode)+0.5d0*gradient(i,inode,igss)*conductivity_new(igss)%val(i,k)*gradient(k,knode,igss)
        enddo
       enddo
       matrix(inode,knode)=matrix(inode,knode)*wght(igss)
       jdof=correspond_therm(idof,kdof)
       !$omp atomic update
       conductance(jdof)=conductance(jdof)+matrix(inode,knode)
      enddo
     enddo
    enddo
   enddo
   !$omp end do
! *** evaluates mechanical residual
   !$omp do schedule(static) private(ielem,intg,igss,vector,inode,iddl,idof,i,j)
   do ielem=1,nelem3d
    do intg=1,nntg_per_element3d
     igss=(ielem-1)*nntg_per_element3d+intg
     do inode=1,nnode_per_element3d
      do i=1,3
       idof=3*(connect3d(inode,ielem)-1)+i
       iddl=3*(inode-1)+i
       vector(iddl)=0.d0
       do j=1,3
        vector(iddl)=vector(iddl)+stress_new(igss)%val(i,j)*gradient(j,inode,igss)
       enddo
       vector(iddl)=vector(iddl)*wght(igss)
       !$omp atomic update
       residual_mech(idof)=residual_mech(idof)+vector(iddl)
      enddo
     enddo
    enddo
   enddo
   !$omp end do
! *** evaluates thermal residual
   !$omp do schedule(static) private(ielem,intg,igss,vector,inode,idof,i)
   do ielem=1,nelem3d
    do intg=1,nntg_per_element3d
     igss=(ielem-1)*nntg_per_element3d+intg
     do inode=1,nnode_per_element3d
      idof=connect3d(inode,ielem)
      vector(inode)=interpolation(inode,igss)*(capacity_new(igss)*(temp_new(igss)-temp_old(igss))/(time_new-time_old)-source_new(igss))
!       vector(inode)=interpolation(inode,igss)*capacity_new(igss)*(temp_new(igss)-temp_old(igss))/(time_new-time_old)
      do i=1,3
       vector(inode)=vector(inode)-0.5d0*gradient(i,inode,igss)*(heatflow_new(igss)%val(i)+heatflow_old(igss)%val(i))
      enddo
      vector(inode)=vector(inode)*wght(igss)
      !$omp atomic update
      residual_therm(idof)=residual_therm(idof)+vector(inode)
     enddo
    enddo
   enddo
   !$omp end do
! *** prescribes mechanical boundary conditions
   !$omp do schedule(static) private(idof,jdof,kdof)
   do idof=1,ndof
    if (bc_mech_type(idof).eq.1) then
     residual_mech(idof)=0.d0
     do jdof=1,ndof
      kdof=correspond_mech(idof,jdof)
      stiffness(kdof)=0.d0
     enddo
     kdof=correspond_mech(idof,idof)
     stiffness(kdof)=1.d0
    endif
    if (bc_mech_type(idof).eq.2) then
     !$omp atomic update
     residual_mech(idof)=residual_mech(idof)-bc_mech_val(1,idof)-(bc_mech_val(2,idof)-bc_mech_val(1,idof))*progress
    endif
   enddo
   !$omp end do
! *** prescribes thermal boundary conditions
   !$omp do schedule(static) private(inode,jnode,knode)
   do inode=1,nnode
    if (bc_therm_type(inode).eq.1) then
     residual_therm(inode)=0.d0
     do jnode=1,nnode
      knode=correspond_therm(inode,jnode)
      conductance(knode)=0.d0
     enddo
     knode=correspond_therm(inode,inode)
     conductance(knode)=1.d0
    endif
    if (bc_therm_type(inode).eq.2) then
     !$omp atomic update
     residual_therm(inode)=residual_therm(inode)-bc_therm_val(1,inode)-(bc_therm_val(2,inode)-bc_therm_val(1,inode))*progress
    endif
   enddo
   !$omp end do
   !$omp end parallel
! *** calculates the residuals
   res_mech_old=res_mech
   res_mech=0.d0
   do idof=1,ndof
    if (bc_mech_type(idof).ne.1) then
     res_mech=res_mech+residual_mech(idof)**2.d0
    endif
   enddo
   res_mech=sqrt(res_mech)
   if (res_mech.gt.res_mech_old) then
    damp_mech=0.9d0*damp_mech
   endif
   write(*,'("  ITER:",i4," MECHANICAL RESIDUAL:",es12.4)') iter,res_mech
   res_therm_old=res_therm
   res_therm=0.d0
   do inode=1,nnode
    if (bc_therm_type(inode).ne.1) then
     res_therm=res_therm+residual_therm(inode)**2.d0
    endif
   enddo
   res_therm=sqrt(res_therm)
   if (res_therm.gt.res_therm_old) then
    damp_therm=0.9d0*damp_therm
   endif
   write(*,'("            THERMAL RESIDUAL:",es12.4)') res_therm
! *** uses Newton-Raphson algorithm to update solution with PARDISO solver
   if (istp.eq.1.and.iter.eq.1) then
    phase = 11 ! *** analysis
    call pardiso(pt_mech,maxfct_mech,mnum_mech,mtype_mech,phase,ndof,stiffness,ia_mech,ja_mech,permutation_mech,nrhs,iparm_mech,msglvl,residual_mech,increment_mech,error_mech)
    call pardiso(pt_therm,maxfct_therm,mnum_therm,mtype_therm,phase,nnode,conductance,ia_therm,ja_therm,permutation_therm,nrhs,iparm_therm,msglvl,residual_therm,increment_therm,error_therm)
   endif
   phase=23 ! *** factorization + solution
   call pardiso(pt_mech,maxfct_mech,mnum_mech,mtype_mech,phase,ndof,stiffness,ia_mech,ja_mech,permutation_mech,nrhs,iparm_mech,msglvl,residual_mech,increment_mech,error_mech)
   call pardiso(pt_therm,maxfct_therm,mnum_therm,mtype_therm,phase,nnode,conductance,ia_therm,ja_therm,permutation_therm,nrhs,iparm_therm,msglvl,residual_therm,increment_therm,error_therm)
   !$omp parallel default(shared)
   !$omp do schedule(static) private(idof)
   do idof=1,ndof
    position_new(idof)=position_new(idof)-increment_mech(idof)*damp_mech
   enddo
   !$omp end do
   !$omp do schedule(static) private(inode)
   do inode=1,nnode
    temperature_new(inode)=temperature_new(inode)-increment_therm(inode)*damp_therm
   enddo
   !$omp end do
   !$omp end parallel
  enddo

!   open(unit=20, file="temperature.out", status="replace", action="write")
!   do inode=1,nnode
!    write(20,'(i6,2es12.4)') inode,position_new(3*inode),temperature_new(inode)
!   enddo
!   close(unit=20)
!
!    stop
! *** updates variables
  time_old=time_new
  !$omp parallel default(shared)
  !$omp do schedule(static) private(idof,delta)
  do idof=1,ndof
   delta=position_new(idof)-position_old(idof)
   position_old(idof)=position_new(idof)
   position_new(idof)=position_old(idof)+delta
  enddo
  !$omp end do
  !$omp do schedule(static) private(inode,delta)
  do inode=1,nnode
   delta=temperature_new(inode)-temperature_old(inode)
   temperature_old(inode)=temperature_new(inode)
   temperature_new(inode)=temperature_old(inode)+delta
  enddo
  !$omp end do
  !$omp do schedule(static) private(ielem,intg,igss,istat,delta)
  do ielem=1,nelem3d
   do intg=1,nntg_per_element3d
    igss=(ielem-1)*nntg_per_element3d+intg
    internalenergy_old(igss)=internalenergy_new(igss)
    dissipatedenergy_old(igss)=dissipatedenergy_new(igss)
    stress_old(igss)=stress_new(igss)
    defgrad_old(igss)=defgrad_new(igss)
    temp_old(igss)=temp_new(igss)
    tempgrad_old(igss)=tempgrad_new(igss)
    heatflow_old(igss)=heatflow_new(igss)
    do istat=1,nstat
     delta=state_new(istat,igss)-state_old(istat,igss)
     state_old(istat,igss)=state_new(istat,igss)
     state_new(istat,igss)=state_old(istat,igss)+delta
    enddo
   enddo
  enddo
  !$omp end do
  !$omp end parallel
  if (freq.ne.0) then
   if (modulo(istp,freq).eq.0) then
    ifile=ifile+1
    write(fil,'(i8)') ifile
    fil=adjustl(fil)
    filevtu="result"//trim(fil)//".vtu"
    call export_vtu(filevtu)
    call export_rct(1)
!     stop
   endif
  endif
 enddo
! *** export final results
 call export_vtu("result.vtu")
 call export_rct(1)
 call export_rct(2)
! *** cleans memory for PARDISO solver
 phase = 0
 call pardiso(pt_mech,maxfct_mech,mnum_mech,mtype_mech,phase,ndof,stiffness,ia_mech,ja_mech,permutation_mech,nrhs,iparm_mech,msglvl,residual_mech,increment_mech,error_mech)
 call pardiso(pt_therm,maxfct_therm,mnum_therm,mtype_therm,phase,nnode,conductance,ia_therm,ja_therm,permutation_therm,nrhs,iparm_therm,msglvl,residual_therm,increment_therm,error_therm)

contains
! *******************************************************************
! *** subroutine build_csr_structure
! *******************************************************************
subroutine build_csr_structure
 implicit none
 integer, allocatable :: tag(:,:)
 integer :: cnt
! *** mechanical problem
 allocate(tag(ndof,ndof))
 tag=0
 do ielem=1,nelem3d
  do inode=1,nnode_per_element3d
   do iddl=-2,0
    idof=3*connect3d(inode,ielem)+iddl
    do jnode=1,nnode_per_element3d
     do jddl=-2,0
      jdof=3*connect3d(jnode,ielem)+jddl
      tag(idof,jdof)=1
     enddo
    enddo
   enddo
  enddo
 enddo
! *** finds the number of non-zero coefficients
 nnz=0
 do idof=1,ndof
  do jdof=1,ndof
   if (tag(idof,jdof).ne.0) then
    nnz=nnz+1
   endif
  enddo
 enddo
 write(*,'("NUMBER OF NON-ZERO VALUES IN STIFFNESS MATRIX:",i6)') nnz
 allocate(stiffness(nnz),residual_mech(ndof),increment_mech(ndof))
 allocate(ja_mech(nnz),ia_mech(ndof+1),correspond_mech(ndof,ndof),permutation_mech(ndof))
 cnt=1
 ia_mech(1)=1
 do idof=1,ndof
  do jdof=1,ndof
   if (tag(idof,jdof).ne.0) then
    correspond_mech(idof,jdof)=cnt
    ja_mech(cnt)=jdof
    cnt=cnt+1
   endif
  enddo
  ia_mech(idof+1)=cnt
 enddo
 deallocate(tag)
! *** thermal problem
 allocate(tag(nnode,nnode))
 tag=0
 do ielem=1,nelem3d
  do inode=1,nnode_per_element3d
   idof=connect3d(inode,ielem)
   do jnode=1,nnode_per_element3d
    jdof=connect3d(jnode,ielem)
    tag(idof,jdof)=1
   enddo
  enddo
 enddo
! *** finds the number of non-zero coefficients
 mmz=0
 do inode=1,nnode
  do jnode=1,nnode
   if (tag(inode,jnode).ne.0) then
    mmz=mmz+1
   endif
  enddo
 enddo
 write(*,'("NUMBER OF NON-ZERO VALUES IN CONDUCTANCE MATRIX:",i6)') mmz
 allocate(conductance(mmz),residual_therm(nnode),increment_therm(nnode))
 allocate(ja_therm(mmz),ia_therm(nnode+1),correspond_therm(nnode,nnode),permutation_therm(nnode))
 cnt=1
 ia_therm(1)=1
 do inode=1,nnode
  do jnode=1,nnode
   if (tag(inode,jnode).ne.0) then
    correspond_therm(inode,jnode)=cnt
    ja_therm(cnt)=jnode
    cnt=cnt+1
   endif
  enddo
  ia_therm(inode+1)=cnt
 enddo
 deallocate(tag)

end subroutine
! *******************************************************************
! *** subroutine sort
! *******************************************************************
subroutine sort(arr, n)
 integer, intent(inout) :: arr(:)
 integer, intent(in) :: n
 integer :: i, j, tmp
 do i = 1, n-1
  do j = i+1, n
   if (arr(j) < arr(i)) then
    tmp = arr(i); arr(i) = arr(j); arr(j) = tmp
   end if
  end do
 end do
end subroutine
! *******************************************************************
! *** subroutine export_rct
! *******************************************************************
subroutine export_rct(iopt)
 implicit none
 integer :: iset,iopt
 real(8) :: reaction(3),moment(3),disp(3)
 select case(iopt)
 case(0)
  open(unit=20, file="reaction.out", status="replace", action="write")
  write(20,'(" TIME")',advance='no')
  do iset=1,nnset
   if (iset.eq.nnset) then
    write(20,'(" RF1_",a," RF2_",a," RF3_",a," U1_",a," U2_",a," U3_",a)',advance='yes') trim(nset_names(iset)),trim(nset_names(iset)),trim(nset_names(iset)),trim(nset_names(iset)),trim(nset_names(iset)),trim(nset_names(iset))
   else
    write(20,'(" RF1_",a," RF2_",a," RF3_",a," U1_",a," U2_",a," U3_",a)',advance='no') trim(nset_names(iset)),trim(nset_names(iset)),trim(nset_names(iset)),trim(nset_names(iset)),trim(nset_names(iset)),trim(nset_names(iset))
   endif
  enddo
 case(1)
  write(20,'(1es12.4)',advance='no') time_old
  residual_mech=0.d0
  !$omp parallel default(shared)
  !$omp do schedule(static) private(ielem,intg,igss,vector,inode,iddl,idof,i,j)
  do ielem=1,nelem3d
   do intg=1,nntg_per_element3d
    igss=(ielem-1)*nntg_per_element3d+intg
    do inode=1,nnode_per_element3d
     do i=1,3
      idof=3*(connect3d(inode,ielem)-1)+i
      iddl=3*(inode-1)+i
      vector(iddl)=0.d0
      do j=1,3
       vector(iddl)=vector(iddl)+stress_old(igss)%val(i,j)*gradient(j,inode,igss)*wght(igss)
      enddo
      !$omp atomic update
      residual_mech(idof)=residual_mech(idof)+vector(iddl)
     enddo
    enddo
   enddo
  enddo
  !$omp end do
  !$omp end parallel
  do idof=1,ndof
   if (bc_mech_type(idof).eq.2) then
    residual_mech(idof)=residual_mech(idof)-bc_mech_val(1,idof)-(bc_mech_val(2,idof)-bc_mech_val(1,idof))*progress
   endif
  enddo
  do iset=1,nnset
!    write(*,'(a)') nset_names(iset)
   reaction(1)=0.d0
   reaction(2)=0.d0
   reaction(3)=0.d0
   moment(1)=0.d0
   moment(2)=0.d0
   moment(3)=0.d0
   disp(1)=0.d0
   disp(2)=0.d0
   disp(3)=0.d0
   do inode=nset_start(iset),nset_end(iset)
    jnode=nset_nodes(inode)
    idof=3*jnode-2
    jdof=3*jnode-1
    kdof=3*jnode
    disp(1)=disp(1)+position_old(idof)-coords(1,jnode)
    disp(2)=disp(2)+position_old(jdof)-coords(2,jnode)
    disp(3)=disp(3)+position_old(kdof)-coords(3,jnode)
    if (bc_mech_type(idof).eq.1) then
!      write(*,'("Node",i6," DOF:",i6," Reaction residual_mech:",es12.4)') jnode,idof,residual_mech(idof)
     reaction(1)=reaction(1)+residual_mech(idof)
    endif
    if (bc_mech_type(jdof).eq.1) then
!      write(*,'("Node",i6," DOF:",i6," Reaction residual_mech:",es12.4)') jnode,idof,residual_mech(idof)
     reaction(2)=reaction(2)+residual_mech(jdof)
    endif
    if (bc_mech_type(jdof).eq.1) then
!      write(*,'("Node",i6," DOF:",i6," Reaction force:",es12.4)') jnode,idof,force(idof)
     reaction(3)=reaction(3)+residual_mech(kdof)
    endif
    moment(1)=moment(1)+residual_mech(kdof)*position_old(jdof)-residual_mech(jdof)*position_old(kdof)
    moment(2)=moment(2)+residual_mech(idof)*position_old(kdof)-residual_mech(kdof)*position_old(idof)
    moment(3)=moment(3)+residual_mech(jdof)*position_old(idof)-residual_mech(idof)*position_old(jdof)
   enddo
   disp(1)=disp(1)/real(nset_end(iset)-nset_start(iset)+1,8)
   disp(2)=disp(2)/real(nset_end(iset)-nset_start(iset)+1,8)
   disp(3)=disp(3)/real(nset_end(iset)-nset_start(iset)+1,8)
   if (iset.eq.nnset) then
    write(20,'(6es12.4)',advance='yes') reaction(:),disp(:)
   else
    write(20,'(6es12.4)',advance='no') reaction(:),disp(:)
   endif
  enddo
 case(2)
  close(unit=20)
 end select

end subroutine
! *******************************************************************
! *** subroutine export_vtu
! *******************************************************************
subroutine export_vtu(filename)
 implicit none
 character(*), intent(in) :: filename
 real(8) :: stress(3,3),defgrad(3,3),state(nstat),energy,temp,flux(3)
 integer :: vtk_type
! *** check element type: C3D4 or C3D8
 select case(nnode_per_element3d)
 case(4)
  if (nntg_per_element3d.ne.1) then
   write(*,*) 'Error export_vtu: C3D4 requires nntg_per_element3d=1'
   stop
  endif
  vtk_type=10
 case(8)
  if (nntg_per_element3d.ne.8) then
   write(*,*) 'Error export_vtu: C3D8 requires nntg_per_element3d=8'
   stop
  endif
  vtk_type=12
 case(10)
  if (nntg_per_element3d.ne.4) then
   write(*,*) 'Error export_vtu: C3D10 requires nntg_per_element3d=4'
   stop
  endif
  vtk_type=24
 case default
  write(*,*) 'Error export_vtu: unsupported element type'
  stop
 end select

 open(newunit=unit,file=filename,status="replace",action="write")
  write(unit,'(a)') '<?xml version="1.0"?>'
  write(unit,'(a)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
  write(unit,'(a)') '  <UnstructuredGrid>'
  write(unit,'(a,I0,a,I0,a)') '    <Piece NumberOfPoints="',nnode,'" NumberOfCells="',nelem3d,'">'
  write(unit,'(a)') '      <Points>'
  write(unit,'(a)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
  do inode=1,nnode
   idof=3*inode-2
   jdof=3*inode-1
   kdof=3*inode
   write(unit,'(3E16.6E3)') position_old(idof),position_old(jdof),position_old(kdof)
  enddo
  write(unit,'(a)') '        </DataArray>'
  write(unit,'(a)') '      </Points>'
  write(unit,'(a)') '      <Cells>'
  write(unit,'(a)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
  do ielem=1,nelem3d
   do inode=1,nnode_per_element3d
    write(unit,'(I0,1X)',advance='no') connect3d(inode,ielem)-1
   enddo
   write(unit,*)
  enddo
  write(unit,'(a)') '        </DataArray>'
  write(unit,'(a)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
  do ielem=1,nelem3d
   write(unit,'(I0)') ielem*nnode_per_element3d
  enddo
  write(unit,'(a)') '        </DataArray>'
  write(unit,'(a)') '        <DataArray type="UInt8" Name="types" format="ascii">'
  do ielem=1,nelem3d
   write(unit,'(I0)') vtk_type
  enddo
  write(unit,'(a)') '        </DataArray>'
  write(unit,'(a)') '      </Cells>'
  write(unit,'(a)') '      <CellData Scalars="scalars">'
! *** Compute average stress tensor
  write(unit,'(a)') '        <DataArray type="Float64" Name="Stress" NumberOfComponents="9" format="ascii">'
  do ielem=1,nelem3d
   do i=1,3
    do j=1,3
     stress(i,j)=0.d0
     do intg=1,nntg_per_element3d
      igss=(ielem-1)*nntg_per_element3d+intg
      stress(i,j)=stress(i,j)+wght(igss)*stress_old(igss)%val(i,j)/volume(ielem)
     enddo
    enddo
   enddo
   write(unit,'(9E16.6E3)') ((stress(i,j),j=1,3),i=1,3)
  enddo
  write(unit,'(a)') '        </DataArray>'
! *** Compute average deformation gradient
  write(unit,'(a)') '        <DataArray type="Float64" Name="DeformationGradient" NumberOfComponents="9" format="ascii">'
  do ielem=1,nelem3d
   do i=1,3
    do j=1,3
     defgrad(i,j)=0.d0
     do intg=1,nntg_per_element3d
      igss=(ielem-1)*nntg_per_element3d+intg
      defgrad(i,j)=defgrad(i,j)+wght(igss)*defgrad_old(igss)%val(i,j)/volume(ielem)
     enddo
    enddo
   enddo
   write(unit,'(9E16.6E3)') ((defgrad(i,j),j=1,3),i=1,3)
  enddo
  write(unit,'(a)') '        </DataArray>'
! *** Compute average state variables
  if (nstat.gt.0) then
   write(unit,'(a,I1,a)') '        <DataArray type="Float64" Name="StateVariable" NumberOfComponents="',nstat,'" format="ascii">'
   do ielem=1,nelem3d
    do istat=1,nstat
     state(istat)=0.d0
     do intg=1,nntg_per_element3d
      igss=(ielem-1)*nntg_per_element3d+intg
      state(istat)=state(istat)+wght(igss)*state_old(istat,igss)/volume(ielem)
     enddo
    enddo
    do istat=1,nstat-1
     write(unit,'(1E16.6E3)',advance='no') state(istat)
    enddo
    write(unit,'(1E16.6E3)',advance='yes') state(nstat)
   end do
   write(unit,'(a)') '        </DataArray>'
  endif
! *** Compute average internal energy
  write(unit,'(a)') '        <DataArray type="Float64" Name="InternalEnergy" NumberOfComponents="1" format="ascii">'
  do ielem=1,nelem3d
   energy=0.d0
   do intg=1,nntg_per_element3d
    igss=(ielem-1)*nntg_per_element3d+intg
    energy=energy+wght(igss)*internalenergy_old(igss)/volume(ielem)
   enddo
   write(unit,'(1E16.6E3)') energy
  enddo
  write(unit,'(a)') '        </DataArray>'
! *** Compute average dissipated energy
  write(unit,'(a)') '        <DataArray type="Float64" Name="DissipatedEnergy" NumberOfComponents="1" format="ascii">'
  do ielem=1,nelem3d
   energy=0.d0
   do intg=1,nntg_per_element3d
    igss=(ielem-1)*nntg_per_element3d+intg
    energy=energy+wght(igss)*dissipatedenergy_old(igss)/volume(ielem)
   enddo
   write(unit,'(1E16.6E3)') energy
  enddo
  write(unit,'(a)') '        </DataArray>'
  write(unit,'(a)') '        <DataArray type="Float64" Name="Temperature" NumberOfComponents="1" format="ascii">'
  do ielem=1,nelem3d
   temp=0.d0
   do intg=1,nntg_per_element3d
    igss=(ielem-1)*nntg_per_element3d+intg
    temp=temp+wght(igss)*temp_old(igss)/volume(ielem)
   enddo
   write(unit,'(1E16.6E3)') temp
  enddo
  write(unit,'(a)') '        </DataArray>'
  write(unit,'(a)') '        <DataArray type="Float64" Name="HeatFlux" NumberOfComponents="3" format="ascii">'
  do ielem=1,nelem3d
   do i=1,3
    flux(i)=0.d0
    do intg=1,nntg_per_element3d
     igss=(ielem-1)*nntg_per_element3d+intg
     flux(i)=flux(i)+wght(igss)*heatflow_old(igss)%val(i)/volume(ielem)
    enddo
   enddo
   write(unit,'(3E16.6E3)') (flux(i),i=1,3)
  enddo
  write(unit,'(a)') '        </DataArray>'
  write(unit,'(a)') '      </CellData>'
  write(unit,'(a)') '    </Piece>'
  write(unit,'(a)') '  </UnstructuredGrid>'
  write(unit,'(a)') '</VTKFile>'
 close(unit)

end subroutine export_vtu
! *******************************************************************
! *** subroutine tetrahedronT4
! *******************************************************************
subroutine tetrahedronT4(x,y,z,N,G,L,W,V)
 implicit none
 real(8), intent(in)  :: x(4),y(4),z(4)
 real(8), intent(out) :: W(1),V,L(3,1),N(4,1),G(3,4,1)
 real(8) :: B(9,12,1),dN(4,3)
 real(8) :: dNdx(4),dNdy(4),dNdz(4)
 real(8) :: J(3,3),detJ,invJ(3,3)
 integer :: i

 V=0.d0
 B=0.d0
 L=0.d0
 W=0.d0

 ! *** Shape function derivatives in natural coordinates
 !     N1=1-xi-eta-zeta, N2=xi, N3=eta, N4=zeta
 N(1,1)=0.25d0
 N(2,1)=0.25d0
 N(3,1)=0.25d0
 N(4,1)=0.25d0
 dN(1,:)=(/-1.d0,-1.d0,-1.d0/)
 dN(2,:)=(/ 1.d0, 0.d0, 0.d0/)
 dN(3,:)=(/ 0.d0, 1.d0, 0.d0/)
 dN(4,:)=(/ 0.d0, 0.d0, 1.d0/)

 ! *** Jacobian: J(a,b) = d(coordinate a)/d(natural coordinate b)
 J=0.d0
 do i=1,4
  J(1,1)=J(1,1)+dN(i,1)*x(i)
  J(1,2)=J(1,2)+dN(i,2)*x(i)
  J(1,3)=J(1,3)+dN(i,3)*x(i)
  J(2,1)=J(2,1)+dN(i,1)*y(i)
  J(2,2)=J(2,2)+dN(i,2)*y(i)
  J(2,3)=J(2,3)+dN(i,3)*y(i)
  J(3,1)=J(3,1)+dN(i,1)*z(i)
  J(3,2)=J(3,2)+dN(i,2)*z(i)
  J(3,3)=J(3,3)+dN(i,3)*z(i)
 enddo

 ! *** Jacobian determinant
 detJ = J(1,1)*(J(2,2)*J(3,3)-J(2,3)*J(3,2)) &
      - J(1,2)*(J(2,1)*J(3,3)-J(2,3)*J(3,1)) &
      + J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1))

 if (detJ.le.0.d0) then
  write(*,'("Non-positive Jacobian in tetrahedronT4")')
  stop
 endif

 ! *** Inverse Jacobian
 invJ(1,1)=(J(2,2)*J(3,3)-J(2,3)*J(3,2))/detJ
 invJ(1,2)=(J(1,3)*J(3,2)-J(1,2)*J(3,3))/detJ
 invJ(1,3)=(J(1,2)*J(2,3)-J(1,3)*J(2,2))/detJ
 invJ(2,1)=(J(2,3)*J(3,1)-J(2,1)*J(3,3))/detJ
 invJ(2,2)=(J(1,1)*J(3,3)-J(1,3)*J(3,1))/detJ
 invJ(2,3)=(J(1,3)*J(2,1)-J(1,1)*J(2,3))/detJ
 invJ(3,1)=(J(2,1)*J(3,2)-J(2,2)*J(3,1))/detJ
 invJ(3,2)=(J(1,2)*J(3,1)-J(1,1)*J(3,2))/detJ
 invJ(3,3)=(J(1,1)*J(2,2)-J(1,2)*J(2,1))/detJ

 ! *** Derivatives in physical coordinates
 !     grad_x(N) = transpose(invJ) * grad_nat(N)
 do i=1,4
  dNdx(i)=invJ(1,1)*dN(i,1) &
         +invJ(2,1)*dN(i,2) &
         +invJ(3,1)*dN(i,3)
  dNdy(i)=invJ(1,2)*dN(i,1) &
         +invJ(2,2)*dN(i,2) &
         +invJ(3,2)*dN(i,3)
  dNdz(i)=invJ(1,3)*dN(i,1) &
         +invJ(2,3)*dN(i,2) &
         +invJ(3,3)*dN(i,3)
  G(1,i,1)=dNdx(i)
  G(2,i,1)=dNdy(i)
  G(3,i,1)=dNdz(i)
 enddo

 ! *** Integration point: centroid
 L(1,1)=0.25d0*(x(1)+x(2)+x(3)+x(4))
 L(2,1)=0.25d0*(y(1)+y(2)+y(3)+y(4))
 L(3,1)=0.25d0*(z(1)+z(2)+z(3)+z(4))

 ! *** Displacement-gradient matrix B
 do i=1,4
  B(1,3*i-2,1)=dNdx(i)
  B(2,3*i-1,1)=dNdx(i)
  B(3,3*i  ,1)=dNdx(i)
  B(4,3*i-2,1)=dNdy(i)
  B(5,3*i-1,1)=dNdy(i)
  B(6,3*i  ,1)=dNdy(i)
  B(7,3*i-2,1)=dNdz(i)
  B(8,3*i-1,1)=dNdz(i)
  B(9,3*i  ,1)=dNdz(i)
 enddo

 ! *** Integration weight and volume
 W(1)=detJ/6.d0
 V=W(1)

end subroutine tetrahedronT4
! *******************************************************************
! *** subroutine tetrahedronT10
! *** Quadratic tetrahedron, 10 nodes (C3D10)
! *** 4 integration points
! *******************************************************************
subroutine tetrahedronT10(x,y,z,N,G,L,W,V)
 implicit none
 real(8), intent(in)  :: x(10),y(10),z(10)
 real(8), intent(out) :: W(4),V,L(3,4),N(10,4),G(3,10,4)
 real(8) :: B(9,30,4),dN(10,3)
 real(8) :: dNdx(10),dNdy(10),dNdz(10)
 real(8) :: J(3,3),detJ,invJ(3,3)
 real(8) :: r,s,t,lam(4),dlam(4,3)
 real(8) :: aa,bb
 integer :: i,p

 aa = 0.5854101966249685d0
 bb = 0.1381966011250105d0
 V=0.d0
 B=0.d0
 L=0.d0
 W=0.d0
 ! *** Derivatives of barycentric coordinates
 dlam(1,:)=(/-1.d0,-1.d0,-1.d0/)
 dlam(2,:)=(/ 1.d0, 0.d0, 0.d0/)
 dlam(3,:)=(/ 0.d0, 1.d0, 0.d0/)
 dlam(4,:)=(/ 0.d0, 0.d0, 1.d0/)
 do p=1,4
! *** Integration point in natural coordinates
  select case(p)
  case(1)
   r=bb
   s=bb
   t=bb
  case(2)
   r=aa
   s=bb
   t=bb
  case(3)
   r=bb
   s=aa
   t=bb
  case(4)
   r=bb
   s=bb
   t=aa
  end select
  lam(1)=1.d0-r-s-t
  lam(2)=r
  lam(3)=s
  lam(4)=t
! *** Shape functions
  do i=1,4
   N(i,p)=lam(i)*(2.d0*lam(i)-1.d0)
  enddo
  N(5,p)=4.d0*lam(1)*lam(2)
  N(6,p)=4.d0*lam(2)*lam(3)
  N(7,p)=4.d0*lam(3)*lam(1)
  N(8,p)=4.d0*lam(1)*lam(4)
  N(9,p)=4.d0*lam(2)*lam(4)
  N(10,p)=4.d0*lam(3)*lam(4)
! *** Physical coordinates of integration point
  L(:,p)=0.d0
  do i=1,10
   L(1,p)=L(1,p)+N(i,p)*x(i)
   L(2,p)=L(2,p)+N(i,p)*y(i)
   L(3,p)=L(3,p)+N(i,p)*z(i)
  enddo
! *** Derivatives of shape functions in natural coordinates
  do i=1,4
   dN(i,:)=(4.d0*lam(i)-1.d0)*dlam(i,:)
  enddo
  dN(5,:)=4.d0*(dlam(1,:)*lam(2)+lam(1)*dlam(2,:))
  dN(6,:)=4.d0*(dlam(2,:)*lam(3)+lam(2)*dlam(3,:))
  dN(7,:)=4.d0*(dlam(3,:)*lam(1)+lam(3)*dlam(1,:))
  dN(8,:)=4.d0*(dlam(1,:)*lam(4)+lam(1)*dlam(4,:))
  dN(9,:)=4.d0*(dlam(2,:)*lam(4)+lam(2)*dlam(4,:))
  dN(10,:)=4.d0*(dlam(3,:)*lam(4)+lam(3)*dlam(4,:))
! *** Jacobian: J(a,b) = d(coordinate a)/d(natural coordinate b)
  J=0.d0
  do i=1,10
   J(1,1)=J(1,1)+dN(i,1)*x(i)
   J(1,2)=J(1,2)+dN(i,2)*x(i)
   J(1,3)=J(1,3)+dN(i,3)*x(i)
   J(2,1)=J(2,1)+dN(i,1)*y(i)
   J(2,2)=J(2,2)+dN(i,2)*y(i)
   J(2,3)=J(2,3)+dN(i,3)*y(i)
   J(3,1)=J(3,1)+dN(i,1)*z(i)
   J(3,2)=J(3,2)+dN(i,2)*z(i)
   J(3,3)=J(3,3)+dN(i,3)*z(i)
  enddo
! *** Jacobian determinant
  detJ = J(1,1)*(J(2,2)*J(3,3)-J(2,3)*J(3,2)) &
       - J(1,2)*(J(2,1)*J(3,3)-J(2,3)*J(3,1)) &
       + J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1))
  if (detJ.le.0.d0) then
   write(*,'("Non-positive Jacobian in tetrahedronT10, point ",I2)') p
   write(*,*) 'detJ = ',detJ
   stop
  endif
! *** Inverse Jacobian
  invJ(1,1)=(J(2,2)*J(3,3)-J(2,3)*J(3,2))/detJ
  invJ(1,2)=(J(1,3)*J(3,2)-J(1,2)*J(3,3))/detJ
  invJ(1,3)=(J(1,2)*J(2,3)-J(1,3)*J(2,2))/detJ
  invJ(2,1)=(J(2,3)*J(3,1)-J(2,1)*J(3,3))/detJ
  invJ(2,2)=(J(1,1)*J(3,3)-J(1,3)*J(3,1))/detJ
  invJ(2,3)=(J(1,3)*J(2,1)-J(1,1)*J(2,3))/detJ
  invJ(3,1)=(J(2,1)*J(3,2)-J(2,2)*J(3,1))/detJ
  invJ(3,2)=(J(1,2)*J(3,1)-J(1,1)*J(3,2))/detJ
  invJ(3,3)=(J(1,1)*J(2,2)-J(1,2)*J(2,1))/detJ
! *** Derivatives in physical coordinates
  do i=1,10
   dNdx(i)=invJ(1,1)*dN(i,1) &
          +invJ(2,1)*dN(i,2) &
          +invJ(3,1)*dN(i,3)
   dNdy(i)=invJ(1,2)*dN(i,1) &
          +invJ(2,2)*dN(i,2) &
          +invJ(3,2)*dN(i,3)
   dNdz(i)=invJ(1,3)*dN(i,1) &
          +invJ(2,3)*dN(i,2) &
          +invJ(3,3)*dN(i,3)
   G(1,i,p)=dNdx(i)
   G(2,i,p)=dNdy(i)
   G(3,i,p)=dNdz(i)
  enddo
! *** Displacement-gradient matrix B
  B(:,:,p)=0.d0
  do i=1,10
   B(1,3*i-2,p)=dNdx(i)
   B(2,3*i-1,p)=dNdx(i)
   B(3,3*i  ,p)=dNdx(i)
   B(4,3*i-2,p)=dNdy(i)
   B(5,3*i-1,p)=dNdy(i)
   B(6,3*i  ,p)=dNdy(i)
   B(7,3*i-2,p)=dNdz(i)
   B(8,3*i-1,p)=dNdz(i)
   B(9,3*i  ,p)=dNdz(i)
  enddo
! *** Integration weight and volume
  W(p)=detJ/24.d0
  V=V+W(p)
 enddo

end subroutine tetrahedronT10
! *******************************************************************
! *** subroutine hexahedronH8
! *******************************************************************
subroutine hexahedronH8(x,y,z,N,G,L,W,V)
 implicit none
 real(8), intent(in)  :: x(8),y(8),z(8)
 real(8), intent(out) :: W(8),V,L(3,8),N(8,8),G(3,8,8)
 real(8) :: B(9,24,8),dN(8,3)
 real(8) :: dNdx(8),dNdy(8),dNdz(8)
 real(8) :: J(3,3),detJ,invJ(3,3)
 real(8) :: xi,eta,zeta
 integer :: i,p
 real(8), parameter :: q=1.d0/sqrt(3.d0)
 integer, parameter :: sx(8)=(/ -1, 1, 1,-1,-1, 1, 1,-1 /)
 integer, parameter :: sy(8)=(/ -1,-1, 1, 1,-1,-1, 1, 1 /)
 integer, parameter :: sz(8)=(/ -1,-1,-1,-1, 1, 1, 1, 1 /)

 V=0.d0
 ! *** 8 Gauss integration points
 do p=1,8
  xi  = q*dble(sx(p))
  eta = q*dble(sy(p))
  zeta= q*dble(sz(p))
 ! *** Evaluates shape functions
  do i=1,8
   N(i,p)=0.125d0*(1.d0+dble(sx(i))*xi)*(1.d0+dble(sy(i))*eta)*(1.d0+dble(sz(i))*zeta)
  enddo
 ! *** Physical coordinates of the Gauss point
  L(:,p)=0.d0
  do i=1,8
   L(1,p)=L(1,p)+N(i,p)*x(i)
   L(2,p)=L(2,p)+N(i,p)*y(i)
   L(3,p)=L(3,p)+N(i,p)*z(i)
  enddo
 ! *** Derivatives of shape functions in natural coordinates
  do i=1,8
   dN(i,1)=0.125d0*dble(sx(i))*(1.d0+dble(sy(i))*eta)*(1.d0+dble(sz(i))*zeta)
   dN(i,2)=0.125d0*dble(sy(i))*(1.d0+dble(sx(i))*xi)*(1.d0+dble(sz(i))*zeta)
   dN(i,3)=0.125d0*dble(sz(i))*(1.d0+dble(sx(i))*xi)*(1.d0+dble(sy(i))*eta)
  enddo
 ! *** Jacobian: J(a,b) = d(coordinate a)/d(natural coordinate b)
  J=0.d0
  do i=1,8
   J(1,1)=J(1,1)+dN(i,1)*x(i)
   J(1,2)=J(1,2)+dN(i,2)*x(i)
   J(1,3)=J(1,3)+dN(i,3)*x(i)
   J(2,1)=J(2,1)+dN(i,1)*y(i)
   J(2,2)=J(2,2)+dN(i,2)*y(i)
   J(2,3)=J(2,3)+dN(i,3)*y(i)
   J(3,1)=J(3,1)+dN(i,1)*z(i)
   J(3,2)=J(3,2)+dN(i,2)*z(i)
   J(3,3)=J(3,3)+dN(i,3)*z(i)
  enddo
 ! *** Jacobian determinant
  detJ = J(1,1)*(J(2,2)*J(3,3)-J(2,3)*J(3,2)) &
       - J(1,2)*(J(2,1)*J(3,3)-J(2,3)*J(3,1)) &
       + J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1))
  if (detJ.le.0.d0) then
   write(*,'("Non-positive Jacobian in hexahedronH8, point ",I2)') p
   stop
  endif
! *** Inverse Jacobian
  invJ(1,1)=(J(2,2)*J(3,3)-J(2,3)*J(3,2))/detJ
  invJ(1,2)=(J(1,3)*J(3,2)-J(1,2)*J(3,3))/detJ
  invJ(1,3)=(J(1,2)*J(2,3)-J(1,3)*J(2,2))/detJ
  invJ(2,1)=(J(2,3)*J(3,1)-J(2,1)*J(3,3))/detJ
  invJ(2,2)=(J(1,1)*J(3,3)-J(1,3)*J(3,1))/detJ
  invJ(2,3)=(J(1,3)*J(2,1)-J(1,1)*J(2,3))/detJ
  invJ(3,1)=(J(2,1)*J(3,2)-J(2,2)*J(3,1))/detJ
  invJ(3,2)=(J(1,2)*J(3,1)-J(1,1)*J(3,2))/detJ
  invJ(3,3)=(J(1,1)*J(2,2)-J(1,2)*J(2,1))/detJ
! *** Derivatives in physical coordinates
 !     grad_x(N) = transpose(invJ) * grad_nat(N)
  do i=1,8
   dNdx(i)=invJ(1,1)*dN(i,1) &
          +invJ(2,1)*dN(i,2) &
          +invJ(3,1)*dN(i,3)
   dNdy(i)=invJ(1,2)*dN(i,1) &
          +invJ(2,2)*dN(i,2) &
          +invJ(3,2)*dN(i,3)
   dNdz(i)=invJ(1,3)*dN(i,1) &
          +invJ(2,3)*dN(i,2) &
          +invJ(3,3)*dN(i,3)
   G(1,i,p)=dNdx(i)
   G(2,i,p)=dNdy(i)
   G(3,i,p)=dNdz(i)
  enddo
 ! *** Displacement-gradient matrix B
  B(:,:,p)=0.d0
  do i=1,8
   B(1,3*i-2,p)=dNdx(i)
   B(2,3*i-1,p)=dNdx(i)
   B(3,3*i  ,p)=dNdx(i)
   B(4,3*i-2,p)=dNdy(i)
   B(5,3*i-1,p)=dNdy(i)
   B(6,3*i  ,p)=dNdy(i)
   B(7,3*i-2,p)=dNdz(i)
   B(8,3*i-1,p)=dNdz(i)
   B(9,3*i  ,p)=dNdz(i)
  enddo
! *** Integration wght
  W(p)=detJ
  V=V+W(p)
 enddo

end subroutine hexahedronH8
! *******************************************************************
! *** subroutine quadrangleQ4
! *******************************************************************
subroutine quadrangleQ4(x,y,z,N,W,S)
 implicit none
 real(8), intent(in)  :: x(4), y(4), z(4)
 real(8), intent(out) :: N(4,4), W(4), S
 real(8) :: xi(4), eta(4)
 real(8) :: a
 real(8) :: dNdxi(4), dNdeta(4)
 real(8) :: xr,xs,yr,ys,zr,zs
 real(8) :: cx,cy,cz,Js
 integer :: g,i

! *** local coordinates and weights of integration points
 a = 1.d0/sqrt(3.d0)
 xi(1) = -a
 xi(2) =  a
 xi(3) =  a
 xi(4) = -a
 eta(1) = -a
 eta(2) = -a
 eta(3) =  a
 eta(4) =  a
 S = 0.d0
! *** loop over integration points
 do g = 1,4
! *** interpolation functions
  N(1,g) = 0.25d0*(1.d0-xi(g))*(1.d0-eta(g))
  N(2,g) = 0.25d0*(1.d0+xi(g))*(1.d0-eta(g))
  N(3,g) = 0.25d0*(1.d0+xi(g))*(1.d0+eta(g))
  N(4,g) = 0.25d0*(1.d0-xi(g))*(1.d0+eta(g))
! *** derivatives wrt local coordinates
  dNdxi(1) = -0.25d0*(1.d0-eta(g))
  dNdxi(2) =  0.25d0*(1.d0-eta(g))
  dNdxi(3) =  0.25d0*(1.d0+eta(g))
  dNdxi(4) = -0.25d0*(1.d0+eta(g))
  dNdeta(1) = -0.25d0*(1.d0-xi(g))
  dNdeta(2) = -0.25d0*(1.d0+xi(g))
  dNdeta(3) =  0.25d0*(1.d0+xi(g))
  dNdeta(4) =  0.25d0*(1.d0-xi(g))
! *** jacobian
  xr = 0.d0
  xs = 0.d0
  yr = 0.d0
  ys = 0.d0
  zr = 0.d0
  zs = 0.d0
  do i = 1,4
   xr = xr+dNdxi(i)*x(i)
   xs = xs+dNdeta(i)*x(i)
   yr = yr+dNdxi(i)*y(i)
   ys = ys+dNdeta(i)*y(i)
   zr = zr+dNdxi(i)*z(i)
   zs = zs+dNdeta(i)*z(i)
  enddo
! *** cross product
  cx = yr*zs-zr*ys
  cy = zr*xs-xr*zs
  cz = xr*ys-yr*xs
  Js = sqrt(cx*cx+cy*cy+cz*cz)
! *** physical weights and total area
  W(g) = Js
  S = S+W(g)
 enddo

end subroutine quadrangleQ4
! *******************************************************************
! *** subroutine triangleT6
! *******************************************************************
subroutine triangleT6(x,y,z,N,W,S)
 implicit none
 real(8), intent(in)  :: x(6), y(6), z(6)
 real(8), intent(out) :: N(6,3), W(3), S
 real(8) :: xi(3), eta(3), wg(3)
 real(8) :: L1,L2,L3
 real(8) :: dNdxi(6), dNdeta(6)
 real(8) :: xr,xs,yr,ys,zr,zs
 real(8) :: cx,cy,cz,Js
 integer :: g,a
! *** local coordinates and weights of integration points
 xi(1) = 1.d0/6.d0
 xi(2) = 2.d0/3.d0
 xi(3) = 1.d0/6.d0
 eta(1) = 1.d0/6.d0
 eta(2) = 1.d0/6.d0
 eta(3) = 2.d0/3.d0
 wg(:) = 1.d0/6.d0
 S = 0.d0
! *** loop over integration points
 do g = 1,3
  L1 = 1.d0-xi(g)-eta(g)
  L2 = xi(g)
  L3 = eta(g)
! *** interpolation functions
  N(1,g) = L1*(2.d0*L1-1.d0)
  N(2,g) = L2*(2.d0*L2-1.d0)
  N(3,g) = L3*(2.d0*L3-1.d0)
  N(4,g) = 4.d0*L1*L2
  N(5,g) = 4.d0*L2*L3
  N(6,g) = 4.d0*L3*L1
! *** derivatives wrt local coordinates
  dNdxi(1) = -(4.d0*L1-1.d0)
  dNdxi(2) =  4.d0*L2-1.d0
  dNdxi(3) =  0.d0
  dNdxi(4) =  4.d0*(L1-L2)
  dNdxi(5) =  4.d0*L3
  dNdxi(6) = -4.d0*L3
  dNdeta(1) = -(4.d0*L1-1.d0)
  dNdeta(2) =  0.d0
  dNdeta(3) =  4.d0*L3-1.d0
  dNdeta(4) = -4.d0*L2
  dNdeta(5) =  4.d0*L2
  dNdeta(6) =  4.d0*(L1-L3)
! *** jacobian
  xr = 0.d0
  xs = 0.d0
  yr = 0.d0
  ys = 0.d0
  zr = 0.d0
  zs = 0.d0
  do a = 1,6
   xr = xr + dNdxi(a)*x(a)
   xs = xs + dNdeta(a)*x(a)
   yr = yr + dNdxi(a)*y(a)
   ys = ys + dNdeta(a)*y(a)
   zr = zr + dNdxi(a)*z(a)
   zs = zs + dNdeta(a)*z(a)
  enddo
! *** cross product
  cx = yr*zs-zr*ys
  cy = zr*xs-xr*zs
  cz = xr*ys-yr*xs
  Js = sqrt(cx*cx+cy*cy+cz*cz)
! *** physical weights and total surface
  W(g) = wg(g)*Js
  S = S+W(g)
 enddo

end subroutine triangleT6
! *******************************************************************
! *** subroutine triangleT3
! *******************************************************************
subroutine triangleT3(x,y,z,N,W,S)
 implicit none
 real(8), intent(in)  :: x(3), y(3), z(3)
 real(8), intent(out) :: N(3,1), W(1), S
 real(8) :: xr,yr,zr
 real(8) :: xs,ys,zs
 real(8) :: cx,cy,cz,Js
! *** interpolation functions
  N(1,1) = 1.d0/3.d0
  N(2,1) = 1.d0/3.d0
  N(3,1) = 1.d0/3.d0
! *** derivatives wrt local coordinates
  xr = x(2)-x(1)
  yr = y(2)-y(1)
  zr = z(2)-z(1)
  xs = x(3)-x(1)
  ys = y(3)-y(1)
  zs = z(3)-z(1)
! *** cross product
  cx = yr*zs-zr*ys
  cy = zr*xs-xr*zs
  cz = xr*ys-yr*xs
  Js = sqrt(cx*cx+cy*cy+cz*cz)
! *** physical weights and total surface
  W(1) = Js/2.d0
  S = W(1)

end subroutine triangleT3
! *******************************************************************
! *** subroutine read_material
! *******************************************************************
subroutine read_material(filename,cname,nstat,nprop,props)
 character(5), intent(out) :: cname
 character(80) :: filename, line
 integer, intent(out) :: nstat,nprop
 real(8), allocatable, intent(out) :: props(:)
 integer :: unit,iprop
 open(newunit=unit, file=filename, status='old')
  read(unit,'(a)') line
  read(unit,'(a)') cname
  read(unit,'(a)') line
  read(unit,*) nstat
  read(unit,'(a)') line
  read(unit,*) nprop
  allocate(props(nprop))
  read(unit,'(a)') line
  do iprop=1,nprop
   read(unit,*) props(iprop)
  enddo
 close(unit)
end subroutine
! *******************************************************************
! *** subroutine read_msh
! *******************************************************************
subroutine read_msh(filename, coords, connect3d, nnode, nelem3d)
 implicit none
 character(len=*), intent(in) :: filename
 integer, intent(out) :: nnode, nelem3d
 real(8), allocatable, intent(out) :: coords(:,:)
 integer, allocatable, intent(out) :: connect3d(:,:)
 integer :: i, id, etype, ntags, tag1, tag2
 integer :: node1, node2, node3
 character(len=256) :: line
 integer :: unit, num_elem, dummy
 open(newunit=unit, file=filename, status='old')
! Aller à la section $Nodes
 do
  read(unit,'(A)',end=100) line
  if (trim(line) == '$Nodes') exit
 end do
 read(unit,*) nnode
 allocate(coords(2, nnode))
 do i = 1, nnode
  read(unit,*) id, coords(1, id), coords(2, id), dummy
 end do
! Aller à la section $Elements
  do
   read(unit,'(A)',end=100) line
   if (trim(line) == '$Elements') exit
  end do
  read(unit,*) num_elem
  allocate(connect3d(3, num_elem)) ! sur-allocation temporaire
  nelem3d = 0
  do i = 1, num_elem
   read(unit,*) id, etype, ntags, tag1, tag2, node1, node2, node3
   if (etype == 2) then ! triangle linéaire
    nelem3d = nelem3d + 1
    connect3d(1, nelem3d) = node1
    connect3d(2, nelem3d) = node2
    connect3d(3, nelem3d) = node3
   end if
  end do
  if (nelem3d < num_elem) then
   call move_alloc(connect3d, connect3d)
   allocate(connect3d(3, nelem3d))
   connect3d = connect3d(:, 1:nelem3d)
  end if

 close(unit)
 return
100 print *, "Erreur lors de la lecture du fichier .msh"
 stop
end subroutine read_msh
! *******************************************************************
! *** subroutine solve_system
! *******************************************************************
subroutine solve_system(k, f, u)
 implicit none
 real(8), intent(in) :: k(:,:), f(:)
 real(8), intent(out) :: u(:)
 integer :: n, info
 integer, allocatable :: ipiv(:)
 n = size(f)
 allocate(ipiv(n))
 u=f
 call dgesv(n, 1, k, n, ipiv, u, n, info)
 if (info.ne.0) then
  print *, "Erreur lors de la résolution : info =", info
  stop
 end if
end subroutine solve_system
! *******************************************************************
! *** subroutine read_step
! *******************************************************************
subroutine read_step(filename)
 implicit none
 character(len=*), intent(in) :: filename
 character(len=256) :: line
 logical :: reading_step,reading_disp,reading_force,reading_pressure,reading_temp,reading_flux
 integer :: i,j,k,nbc,iset,inode,jnode,bc_axis,jelem
 character(len=64) :: bc_name
 real(8) :: bc_value1,bc_value2
 real(8) :: area
 real(8), allocatable :: interpolation2D(:,:),wght2D(:)

 allocate(interpolation2D(nnode_per_element2d,nntg_per_element2d),wght2D(nntg_per_element2d))
 allocate(x(nnode_per_element2d),y(nnode_per_element2d),z(nnode_per_element2d))
 reading_step=.false.
 reading_disp=.false.
 reading_force=.false.
 reading_pressure=.false.
 reading_temp=.false.
 reading_flux=.false.
 bc_mech_type(:)=0
 bc_mech_val(:,:)=0.d0
 bc_therm_type(:)=0
 bc_therm_val(:,:)=0.d0
 open(unit=10, file=filename, status='old')
 do
  read(10,'(A)', end=100) line
! Détection des sections
  if (index(line, "*STEP") > 0) then
   reading_step=.true.
   cycle
  endif
  if (index(line, "*DISPLACEMENT") > 0) then
   reading_disp=.true.
   reading_force=.false.
   reading_pressure=.false.
   reading_temp=.false.
   reading_flux=.false.
   cycle
  endif
  if (index(line, "*FORCE") > 0) then
   reading_disp=.false.
   reading_force=.true.
   reading_pressure=.false.
   reading_temp=.false.
   reading_flux=.false.
   cycle
  endif
  if (index(line, "*PRESSURE") > 0) then
   reading_disp=.false.
   reading_force=.false.
   reading_pressure=.true.
   reading_temp=.false.
   reading_flux=.false.
   cycle
  endif
  if (index(line, "*TEMPERATURE") > 0) then
   reading_disp=.false.
   reading_force=.false.
   reading_pressure=.false.
   reading_temp=.true.
   reading_flux=.false.
   cycle
  endif
  if (index(line, "*FLUX") > 0) then
   reading_disp=.false.
   reading_force=.false.
   reading_pressure=.false.
   reading_temp=.false.
   reading_flux=.true.
   cycle
  endif
  if (reading_step) then
   read(line,*) duration,dtime,dtime_min,dtime_max
   write(*,'("Duration",es12.4)') duration
   write(*,'("Initial time increment",es12.4)') dtime
   write(*,'("Minimal time increment",es12.4)') dtime_min
   write(*,'("Maximal time increment",es12.4)') dtime_max
   reading_step=.false.
  endif
  if (reading_disp) then
   i=index(line,",")
   j=index(line(i+1:len(line)),",")
   k=index(line(i+j+1:len(line)),",")
   bc_name=line(1:i-1)
   write(*,'("  Nset for prescribed displacement: ",a)') trim(bc_name)
   read(line(i+1:i+j),*) bc_axis
   write(*,'("  Axis for prescribed displacement: ",i2)') bc_axis
   read(line(i+j+1:i+j+k),*) bc_value1
   write(*,'("  Initial value for prescribed displacement: ",f12.4)') bc_value1
   read(line(i+j+k+1:len(line)),*) bc_value2
   write(*,'("  Final value for prescribed displacement: ",f12.4)') bc_value2
   do iset=1,nnset
    if (bc_name.eq.nset_names(iset)) then
     write(*,*) "found set",iset
     do inode=nset_start(iset),nset_end(iset)
      jnode=nset_nodes(inode)
!       write(*,*) jnode
      idof=3*jnode+bc_axis-3
      bc_mech_type(idof)=1
      bc_mech_val(1,idof)=bc_value1
      bc_mech_val(2,idof)=bc_value2
     enddo
    endif
   enddo
   reading_disp=.false.
  endif
  if (reading_force) then
   i=index(line,",")
   j=index(line(i+1:len(line)),",")
   k=index(line(i+j+1:len(line)),",")
   bc_name=line(1:i-1)
   write(*,'("  Nset for prescribed force: ",a)') trim(bc_name)
   j=index(line(i+1:len(line)),",")
   read(line(i+1:i+j),*) bc_axis
   write(*,'("  Axis for prescribed force: ",i2)') bc_axis
   read(line(i+j+1:i+j+k),*) bc_value1
   write(*,'("  Initial value for prescribed force: ",f12.4)') bc_value1
   read(line(i+j+k+1:len(line)),*) bc_value2
   write(*,'("  Final value for prescribed force: ",f12.4)') bc_value2
   do iset=1,nnset
    if (bc_name.eq.nset_names(iset)) then
     write(*,*) "found set",iset
     do inode=nset_start(iset),nset_end(iset)
      jnode=nset_nodes(inode)
!       write(*,*) jnode
      idof=3*jnode+bc_axis-3
      bc_mech_type(idof)=2
      bc_mech_val(1,idof)=bc_value1
      bc_mech_val(2,idof)=bc_value2
     enddo
    endif
   enddo
   reading_force=.false.
  endif
  if (reading_pressure) then
   i=index(line,",")
   j=index(line(i+1:len(line)),",")
   k=index(line(i+j+1:len(line)),",")
   bc_name=line(1:i-1)
   write(*,'("  Elset for prescribed pressure: ",a)') trim(bc_name)
   j=index(line(i+1:len(line)),",")
   read(line(i+1:i+j),*) bc_axis
   write(*,'("  Axis for prescribed pressure: ",i2)') bc_axis
   read(line(i+j+1:i+j+k),*) bc_value1
   write(*,'("  Initial value for prescribed pressure: ",f12.4)') bc_value1
   read(line(i+j+k+1:len(line)),*) bc_value2
   write(*,'("  Final value for prescribed pressure: ",f12.4)') bc_value2
   do iset=1,nelset
    if (bc_name.eq.elset_names(iset)) then
     write(*,*) "found set",iset
     do ielem=elset_start(iset),elset_end(iset)
      jelem=elset_elem(ielem)
      do inode=1,nnode_per_element2d
       x(inode)=coords(1,connect2d(inode,jelem))
       y(inode)=coords(2,connect2d(inode,jelem))
       z(inode)=coords(3,connect2d(inode,jelem))
      enddo
      select case(nntg_per_element2d)
       case(1)
        call triangleT3(x,y,z,interpolation2D(:,:),wght2D(:),area)
       case(3)
        call triangleT6(x,y,z,interpolation2D(:,:),wght2D(:),area)
       case(4)
        call quadrangleQ4(x,y,z,interpolation2D(:,:),wght2D(:),area)
       case default
        write(*,'("Unknown type of 2D element. Program will stop !")')
        stop
      end select
      do inode=1,nnode_per_element2d
       jnode=connect2d(inode,jelem)
       idof=3*jnode+bc_axis-3
       bc_mech_type(idof)=2
       do intg=1,nntg_per_element2d
        bc_mech_val(1,idof)=bc_mech_val(1,idof)-interpolation2D(inode,intg)*wght2D(intg)*bc_value1
        bc_mech_val(2,idof)=bc_mech_val(2,idof)-interpolation2D(inode,intg)*wght2D(intg)*bc_value2
       enddo
      enddo
     enddo
    endif
   enddo
   reading_pressure=.false.
  endif
  if (reading_temp) then
   i=index(line,",")
   j=index(line(i+1:len(line)),",")
   bc_name=line(1:i-1)
   write(*,'("  Nset for prescribed temperature: ",a)') trim(bc_name)
   read(line(i+1:i+j),*) bc_value1
   write(*,'("  Initial value for prescribed temperature: ",f12.4)') bc_value1
   read(line(i+j+1:len(line)),*) bc_value2
   write(*,'("  Final value for prescribed temperature: ",f12.4)') bc_value2
   do iset=1,nnset
    if (bc_name.eq.nset_names(iset)) then
      write(*,*) "found set",iset
     do inode=nset_start(iset),nset_end(iset)
      jnode=nset_nodes(inode)
      bc_therm_type(jnode)=1
      bc_therm_val(1,jnode)=bc_value1
      bc_therm_val(2,jnode)=bc_value2
     enddo
    endif
   enddo
   reading_temp=.false.
  endif
  if (reading_flux) then
   i=index(line,",")
   j=index(line(i+1:len(line)),",")
   bc_name=line(1:i-1)
   write(*,'("  Nset for prescribed flux: ",a)') trim(bc_name)
   read(line(i+1:i+j),*) bc_value1
   write(*,'("  Initial value for prescribed flux: ",f12.4)') bc_value1
   read(line(i+j+1:len(line)),*) bc_value2
   write(*,'("  Final value for prescribed flux: ",f12.4)') bc_value2
   do iset=1,nnset
    if (bc_name.eq.nset_names(iset)) then
      write(*,*) "found set",iset
     do inode=nset_start(iset),nset_end(iset)
      jnode=nset_nodes(inode)
      bc_therm_type(jnode)=2
      bc_therm_val(1,jnode)=bc_value1
      bc_therm_val(2,jnode)=bc_value2
     enddo
    endif
   enddo
   reading_flux=.false.
  endif
 enddo
100 continue
 close(unit=10)
 deallocate(interpolation2D,wght2D,x,y,z)
end subroutine
! *******************************************************************
! *** subroutine read_inp
! *******************************************************************
subroutine read_inp(filename)
 implicit none
 character(len=*), intent(in) :: filename
 integer :: max_nodes,max_elem
 real(8) :: x,y,z
 character(len=256) :: line
 integer :: i, j, k, id, n(10)
 integer :: node_id, ierr
 logical :: reading_node,reading_elem2d,reading_elem3d,reading_nset,reading_elset,found
 character(len=64) :: temp_name

 nnode=0
 nelem2d=0
 nelem3d=0
 nnset=0
 max_nodes=0
 nelset=0
 max_elem=0
 reading_node=.false.
 reading_elem2d=.false.
 reading_elem3d=.false.
 reading_nset=.false.
 reading_elset=.false.
 open(unit=10, file=filename, status='old')
 do
  read(10,'(A)', end=100) line
! Détection des sections
  if (index(line, "*NODE") > 0) then
   reading_node = .true.
   reading_elem2d = .false.
   reading_elem3d = .false.
   reading_nset  = .false.
   reading_elset  = .false.
   cycle
  else if (index(line, "*ELEMENT") > 0 .and. index(line, "CPS3") > 0) then
   reading_node = .false.
   reading_elem2d = .true.
   reading_elem3d = .false.
   reading_nset  = .false.
   reading_elset  = .false.
   nnode_per_element2d=3
   cycle
  else if (index(line, "*ELEMENT") > 0 .and. index(line, "CPS4") > 0) then
   reading_node = .false.
   reading_elem2d = .true.
   reading_elem3d = .false.
   reading_nset  = .false.
   reading_elset  = .false.
   nnode_per_element2d=4
   cycle
  else if (index(line, "*ELEMENT") > 0 .and. index(line, "CPS6") > 0) then
   reading_node = .false.
   reading_elem2d = .true.
   reading_elem3d = .false.
   reading_nset  = .false.
   reading_elset  = .false.
   nnode_per_element2d=6
   cycle
  else if (index(line, "*ELEMENT") > 0 .and. index(line, "C3D8") > 0) then
   reading_node = .false.
   reading_elem2d = .false.
   reading_elem3d = .true.
   reading_nset  = .false.
   reading_elset  = .false.
   nnode_per_element3d=8
   cycle
  else if (index(line, "*ELEMENT") > 0 .and. index(line, "C3D4") > 0) then
   reading_node = .false.
   reading_elem2d = .false.
   reading_elem3d = .true.
   reading_nset  = .false.
   reading_elset  = .false.
   nnode_per_element3d=4
   cycle
  else if (index(line, "*ELEMENT") > 0 .and. index(line, "C3D10") > 0) then
   reading_node = .false.
   reading_elem2d = .false.
   reading_elem3d = .true.
   reading_nset  = .false.
   reading_elset  = .false.
   nnode_per_element3d=10
   cycle
  else if (index(line, "*NSET") > 0) then
   reading_node = .false.
   reading_elem2d = .false.
   reading_elem3d = .false.
   reading_nset  = .true.
   reading_elset  = .false.
! Extraire le nom du NSET
   call parse_nset_name(line,"NSET=",temp_name)
   nnset = nnset + 1
   cycle
  else if (index(line, "*ELSET") > 0) then
   reading_node = .false.
   reading_elem2d = .false.
   reading_elem3d = .false.
   reading_nset  = .false.
   reading_elset  = .true.
! Extraire le nom du NSET
   call parse_nset_name(line,"ELSET=",temp_name)
   nelset = nelset + 1
   cycle
  else if (index(line, "*") == 1) then
   reading_node = .false.
   reading_elem2d = .false.
   reading_elem3d = .false.
   reading_nset  = .false.
   reading_elset  = .false.
   cycle
  end if
! Lecture des données
  if (reading_node) then
   write(*,'("READING NODE")')
   read(line,*,iostat=ierr) id,x,y,z
   if (ierr.eq.0) then
    nnode=nnode+1
   endif
   write(*,'("READING NODE",i4)') id
  else if (reading_elem2d) then
   write(*,'("READING ELEMS")')
   read(line,*,iostat=ierr) id, (n(inode),inode=1,nnode_per_element2d)
   if (ierr.eq.0) then
    nelem2d=nelem2d+1
   endif
  else if (reading_elem3d) then
   write(*,'("READING ELEMS")')
   read(line,*,iostat=ierr) id, (n(inode),inode=1,nnode_per_element3d)
   if (ierr.eq.0) then
    nelem3d=nelem3d+1
   endif
  else if (reading_nset) then
   write(*,'("READING NSET")')
   write(*,'(i4)') count_commas(line)
   write(*,'("READING LINE in NSET")')
   do i = 1, count_commas(line)
    read(line,*) node_id
    max_nodes = max_nodes + 1
   end do
   write(*,*) nnset,max_nodes
  else if (reading_elset) then
   write(*,'("READING ELSET")')
   write(*,'(i4)') count_commas(line)
   write(*,'("READING LINE in ELSET")')
   do i = 1, count_commas(line)
    read(line,*) id
    max_elem = max_elem + 1
   end do
   write(*,*) nelset,max_elem
  end if
 end do
100 continue
 close(10)
 print *, "Number of nodes:", nnode
 print *, "Number of 2D elements:", nelem2d
 print *, "Number of 3D elements:", nelem3d
 print *, "Number of nodes per element:", nnode_per_element3d
 print *, "Number of NSETs:", nnset
 print *, "Total number of nodes in NSETs:", max_nodes
 print *, "Number of ELSETs:", nelset
 print *, "Total number of elements in ELSETs:", max_elem
 allocate(coords(3, nnode))
 allocate(connect3d(nnode_per_element3d,nelem3d),connect2d(nnode_per_element2d,nelem2d))
 allocate(idelem3d(nelem3d),idelem2d(nelem2d))
 allocate(nset_start(nnset),nset_end(nnset),nset_nodes(max_nodes),nset_names(nnset))
 allocate(elset_start(nelset),elset_end(nelset),elset_elem(max_elem),elset_names(nelset))
 nnode=0
 nelem2d=0
 nelem3d=0
 nnset=0
 max_nodes=0
 nelset=0
 max_elem=0
 reading_node=.false.
 reading_elem2d = .false.
 reading_elem3d=.false.
 reading_nset=.false.
 reading_elset=.false.
 open(unit=10, file=filename, status='old')
 do
  read(10,'(A)', end=200) line
! Détection des sections
  if (index(line, "*NODE") > 0) then
   write(*,'("NODES")')
   reading_node = .true.
   reading_elem2d = .false.
   reading_elem3d = .false.
   reading_nset  = .false.
   reading_elset  = .false.
   cycle
  else if (index(line, "*ELEMENT") > 0 .and. index(line, "CPS3") > 0) then
   reading_node = .false.
   reading_elem2d = .true.
   reading_elem3d = .false.
   reading_nset  = .false.
   reading_elset  = .false.
   nnode_per_element2d=3
   nntg_per_element2d=1
   cycle
  else if (index(line, "*ELEMENT") > 0 .and. index(line, "CPS4") > 0) then
   reading_node = .false.
   reading_elem2d = .true.
   reading_elem3d = .false.
   reading_nset  = .false.
   reading_elset  = .false.
   nnode_per_element2d=4
   nntg_per_element2d=4
   cycle
  else if (index(line, "*ELEMENT") > 0 .and. index(line, "CPS6") > 0) then
   reading_node = .false.
   reading_elem2d = .true.
   reading_elem3d = .false.
   reading_nset  = .false.
   reading_elset  = .false.
   nnode_per_element2d=6
   nntg_per_element2d=3
   cycle
  else if (index(line, "*ELEMENT") > 0 .and. index(line, "C3D8") > 0) then
   reading_node = .false.
   reading_elem2d = .false.
   reading_elem3d = .true.
   reading_nset  = .false.
   nnode_per_element3d=8
   nntg_per_element3d=8
   cycle
  else if (index(line, "*ELEMENT") > 0 .and. index(line, "C3D4") > 0) then
   reading_node = .false.
   reading_elem2d = .false.
   reading_elem3d = .true.
   reading_nset  = .false.
   reading_elset  = .false.
   nnode_per_element3d=4
   nntg_per_element3d=1
   cycle
  else if (index(line, "*ELEMENT") > 0 .and. index(line, "C3D10") > 0) then
   reading_node = .false.
   reading_elem2d = .false.
   reading_elem3d = .true.
   reading_nset  = .false.
   reading_elset  = .false.
   nnode_per_element3d=10
   nntg_per_element3d=4
   cycle
  else if (index(line, "*NSET") > 0) then
   write(*,'("NSET")')
   reading_node = .false.
   reading_elem2d = .false.
   reading_elem3d = .false.
   reading_nset  = .true.
   reading_elset  = .false.
! Extraire le nom du NSET
   call parse_nset_name(line,"NSET=",temp_name)
   nnset = nnset + 1
   nset_names(nnset) = temp_name
   nset_start(nnset)=max_nodes+1
   cycle
  else if (index(line, "*ELSET") > 0) then
   reading_node = .false.
   reading_elem2d = .false.
   reading_elem3d = .false.
   reading_nset  = .false.
   reading_elset  = .true.
! Extraire le nom du ELSET
   call parse_nset_name(line,"ELSET=",temp_name)
   nelset = nelset + 1
   elset_names(nelset) = temp_name
   elset_start(nelset) = max_elem+1
   cycle
  else if (index(line, "*") == 1) then
   reading_node = .false.
   reading_elem2d = .false.
   reading_elem3d = .false.
   reading_nset  = .false.
   reading_elset  = .false.
   cycle
  end if
! Lecture des données
  if (reading_node) then
   write(*,'("READING NODE")')
   read(line,*,iostat=ierr) id,x,y,z
   if (ierr.eq.0.and.id.eq.(nnode+1)) then
    nnode=nnode+1
    coords(1,nnode)=x
    coords(2,nnode)=y
    coords(3,nnode)=z
   else
    write(*,'("ERROR WHEN READING NODE ",i6," POSITION")') id
!     stop
   endif
   write(*,'("READING NODE",i4)') id
  else if (reading_elem2d) then
   write(*,'("READING ELEMS")')
   read(line,*,iostat=ierr) id, (n(inode),inode=1,nnode_per_element2d)
   if (ierr.eq.0) then
    nelem2d=nelem2d+1
    idelem2d(nelem2d)=id
    do inode=1,nnode_per_element2d
     connect2d(inode,nelem2d)=n(inode)
    enddo
   else
    write(*,'("ERROR WHEN READING ELEMENT ",i6," NODES")') id
   endif
  else if (reading_elem3d) then
   write(*,'("READING ELEMS")')
   read(line,*,iostat=ierr) id,(n(inode),inode=1,nnode_per_element3d)
   if (ierr.eq.0) then
    nelem3d=nelem3d+1
    idelem3d(nelem3d)=id
    do inode=1,nnode_per_element3d
     connect3d(inode,nelem3d)=n(inode)
    enddo
   else
    write(*,'("ERROR WHEN READING ELEMENT ",i6," NODES")') id
   endif
  else if (reading_nset) then
   write(*,'("READING NSET")')
   write(*,'(i4)') count_commas(line)
   write(*,'("READING LINE in NSET")')
   k=1
   do i = 1, count_commas(line)
    max_nodes = max_nodes + 1
    j=index(line, ",")
    read(line(k:j-1),*) node_id
    line(j:j)=" "
    k=j+1
    nset_nodes(max_nodes)=node_id
   end do
  else if (reading_elset) then
   write(*,'("READING ELSET")')
   write(*,'(i4)') count_commas(line)
   write(*,'("READING LINE in ELSET")')
   k=1
   do i = 1, count_commas(line)
    max_elem = max_elem + 1
    j=index(line, ",")
    read(line(k:j-1),*) id
    line(j:j)=" "
    k=j+1
    elset_elem(max_elem)=id
   end do
  end if
 end do

200 continue
 close(10)

 do id=1,nnode
  write(*,'("Node",i4," Coords",3es12.4)') id,coords(:,id)
 enddo

 do id=1,nelem2d
  write(*,'("Element",i6," Nodes",6i6)') idelem2d(id),connect2d(:,id)
 enddo

 do id=1,nelem3d
  write(*,'("Element",i6," Nodes",10i6)') idelem3d(id),connect3d(:,id)
 enddo


 do id=1,nnset-1
  nset_end(id)=nset_start(id+1)-1
 enddo
 nset_end(nnset)=max_nodes

 do id=1,nelset-1
  elset_end(id)=elset_start(id+1)-1
 enddo
 elset_end(nelset)=max_elem

 do id=1,nnset
  write(*,'(a)') nset_names(id)
  do i=nset_start(id),nset_end(id)
   write(*,'("    Node:",i6)') nset_nodes(i)
  enddo
 enddo

 do id=1,nelset
  write(*,'(a)') elset_names(id)
  do i=elset_start(id),elset_end(id)
   write(*,'("    Element:",i6)') elset_elem(i)
   found=.true.
   ielem=0
   do while (found.and.ielem.lt.nelem2d)
    ielem=ielem+1
    if (elset_elem(i).eq.idelem2d(ielem)) then
     found=.false.
     write(*,'("    ...corresponding 2d element!",i6)') ielem
     elset_elem(i)=ielem
    endif
   enddo
   ielem=0
   do while (found.and.ielem.lt.nelem3d)
    ielem=ielem+1
    if (elset_elem(i).eq.idelem3d(ielem)) then
     found=.false.
     write(*,'("    ...corresponding 3d element!",i6)') ielem
     elset_elem(i)=ielem
    endif
   enddo
   if (found) then
    write(*,'("Cannot find corresponding element in elset ",a)') elset_names(id)
    stop
   endif
  enddo
 enddo

end subroutine
! *******************************************************************
! *** subroutine parse_nset_name
! *******************************************************************
subroutine parse_nset_name(line, tag, name)
 character(len=*), intent(in) :: line,tag
 character(len=*), intent(out) :: name
 integer :: p1
 name = ""
 p1 = index(line, tag)
 if (p1 > 0) then
  read(line(p1+len(tag):), '(A)') name
  p1 = index(name, ",")
  if (p1 > 0) name = name(:p1-1)
 end if
end subroutine
! *******************************************************************
! *** function count_commas
! *******************************************************************
integer function count_commas(line)
 character(len=*), intent(in) :: line
 integer :: i
 count_commas = 0
 do i = 1, len_trim(line)
  if (line(i:i) == ",") then
   count_commas = count_commas + 1
  endif
 end do
end function
! *******************************************************************
! *** subroutine modulus
! *******************************************************************
subroutine modulus(intg)
 use utils
 use tensor
 implicit none
 integer, intent(in) :: intg
 type(Tensor2) :: defgrad_trl,stress_trl
 real(8) :: state_trl(nstat)
 type(Tensor4) :: modulus_trl
 integer :: i,j,k,l
 real(8) :: small=1.d-5

 call plasticity(nprop,nstat,props,time_old,time_new,defgrad_old(intg),defgrad_new(intg),stress_old(intg),stress_new(intg),state_old(:,intg),state_new(:,intg) &
                ,internalenergy_old(intg),internalenergy_new(intg),dissipatedenergy_old(intg),dissipatedenergy_new(intg),modulus_new(intg))
 write(*,'("MODULUS from subroutine",9es12.4)') modulus_new(intg)%val(:,:,:,:)
 do k=1,3
  do l=1,3
   defgrad_trl=defgrad_new(intg)
   defgrad_trl%val(k,l)=defgrad_trl%val(k,l)+small
   call plasticity(nprop,nstat,props,time_old,time_new,defgrad_old(intg),defgrad_trl,stress_old(intg),stress_trl,state_old(:,intg),state_trl(:) &
                  ,internalenergy_old(intg),internalenergy_new(intg),dissipatedenergy_old(intg),dissipatedenergy_new(intg),modulus_new(intg))
   do i=1,3
    do j=1,3
     modulus_trl%val(i,j,k,l)=(stress_trl%val(i,j)-stress_new(intg)%val(i,j))/(defgrad_trl%val(k,l)-defgrad_new(intg)%val(k,l))
    enddo
   enddo
  enddo
 enddo
 write(*,'("MODULUS from numerical procedure",9es12.4)') modulus_trl%val(:,:,:,:)
end subroutine
end program
! *******************************************************************
! *** subroutine average
! *******************************************************************
subroutine average(wght,npair,pair,radius,nntg,field,field_average)
 implicit none
 integer, intent(in) :: nntg,npair
 integer, intent(in) :: pair(2,npair)
 real(8), intent(in) :: field(nntg),wght(nntg),radius
 real(8), intent(out) :: field_average(nntg)
 integer :: intg,jntg,i,ipair
 real(8), parameter :: pi=3.14159265359d0
 do intg=1,nntg
  field_average(intg)=field(intg)
 enddo
 do ipair=1,npair
  intg=pair(1,ipair)
  jntg=pair(2,ipair)
  if (intg.ne.jntg) then
   field_average(intg)=field_average(intg)+(field(jntg)-field(intg))*wght(jntg)/(pi*radius**2.d0)
   field_average(jntg)=field_average(jntg)+(field(intg)-field(jntg))*wght(intg)/(pi*radius**2.d0)
  endif
 enddo
end subroutine
! *******************************************************************
! *** subroutine variance
! *******************************************************************
subroutine variance(wght,npair,pair,radius,nntg,field,field_average,field_variance)
 implicit none
 integer, intent(in) :: nntg,npair
 integer, intent(in) :: pair(2,npair)
 real(8), intent(in) :: field(nntg),field_average(nntg),wght(nntg),radius
 real(8), intent(out) :: field_variance(nntg)
 integer :: intg,jntg,i,ipair
 real(8), parameter :: pi=3.14159265359d0
 do intg=1,nntg
  field_variance(intg)=field(intg)**2.d0-field_average(intg)**2.d0
 enddo
 do ipair=1,npair
  intg=pair(1,ipair)
  jntg=pair(2,ipair)
  if (intg.ne.jntg) then
   field_variance(intg)=field_variance(intg)+(field(jntg)**2.d0-field(intg)**2.d0)*wght(jntg)/(pi*radius**2.d0)
   field_variance(jntg)=field_variance(jntg)+(field(intg)**2.d0-field(jntg)**2.d0)*wght(intg)/(pi*radius**2.d0)
  endif
 enddo
end subroutine
! *******************************************************************
! *** function kronecker
! *******************************************************************
function kronecker(i,j)
 implicit none
 real(8) :: kronecker
 integer :: i,j

 if (i.eq.j) then
  kronecker=1.d0
 else
  kronecker=0.d0
 endif
end function
! *******************************************************************
! *** subroutine thermal
! *******************************************************************
subroutine thermal(nprop,nstat,props,time_old,time_new,temp_old,temp_new,tempgrad_old,tempgrad_new,state_old,state_new,heatflow_old,heatflow_new,capacity_new,conductivity_new)
 use utils
 use tensor
 implicit none
! *** interface
 integer, intent(in) :: nprop,nstat
 real(8), intent(in) :: props(nprop),state_old(nstat),state_new(nstat),time_old,time_new
 real(8), intent(in) :: temp_old,temp_new
 type(Tensor1), intent(in) :: tempgrad_old,tempgrad_new,heatflow_old
 type(Tensor1), intent(out) :: heatflow_new
 type(Tensor2), intent(out) :: conductivity_new
 real(8), intent(out) :: capacity_new

 if (nprop.lt.3) then
  write(*,'("ERROR: INSUFFICENT NUMBER OF STATE VARIABLES AND/OR MATERIAL PROPERTIES")')
  stop
 endif
! *** recover material properties
 capacity_new=props(nprop-2)*props(nprop-1)
 conductivity_new=props(nprop)*id2
 heatflow_new=(-props(nprop))*tempgrad_new

end subroutine
! *******************************************************************
! *** subroutine elasticity
! *******************************************************************
subroutine elasticity(nprop,nstat,props,time_old,time_new,defgrad_old,defgrad_new,piola_old,piola_new,state_old,state_new,internalenergy_old,internalenergy_new,dissipatedenergy_old,dissipatedenergy_new,modulus_new)
 use utils
 use tensor
 implicit none
! *** interface
 integer, intent(in) :: nprop,nstat
 real(8), intent(in) :: props(nprop),state_old(nstat),time_old,time_new
 real(8), intent(out) :: state_new(nstat)
 real(8), intent(in) :: internalenergy_old,dissipatedenergy_old
 real(8), intent(out) :: internalenergy_new,dissipatedenergy_new
 type(Tensor2), intent(in) :: defgrad_old,defgrad_new,piola_old
 type(Tensor2), intent(out) :: piola_new
 type(Tensor4), intent(out) :: modulus_new
! *** internal
 type(Tensor2s) :: cauchy_new,strain_new,stress_new,piola2_new
 type(Tensor4ss) :: derivative_new,tangent2_new,auxiliary_new,tangent_new
 real(8) :: bulk,shear

! *** calculate final Cauchy-Green tensor
 cauchy_new=sym(trsp(defgrad_new).dot.defgrad_new)
! *** calculate final strain tensor
 call isotropic(func="log",TA=cauchy_new,TC=strain_new)
 strain_new=0.5d0*strain_new
! *** calculate final stress tensor
 bulk=props(1)/(3.d0-6.d0*props(2))
 shear=0.5d0*props(1)/(1.d0+props(2))
 stress_new=((3.d0*bulk)*sph(strain_new))+((2.d0*shear)*dev(strain_new))
! *** calculate final dissipated energy density
 dissipatedenergy_new=dissipatedenergy_old
! *** calculate final internal energy density
 internalenergy_new=0.5d0*(stress_new.ddot.strain_new)
! *** calculate PK2 stress tensor
 call isotropic(func="log",TA=cauchy_new,TB=stress_new,dTCdTA=derivative_new,TBddTCddTA=auxiliary_new)
 piola2_new=stress_new.ddot.derivative_new
! *** calculate PK1 stress tensor
 piola_new=defgrad_new.dot.piola2_new
! *** calculate tangent stress tensor
 tangent_new=((3.d0*bulk)*ps4ss)+((2.d0*shear)*pd4ss)
 tangent2_new=(trsp(derivative_new).ddot.(tangent_new.ddot.derivative_new))+(2.d0*auxiliary_new)
 modulus_new=rtrsp(defgrad_new.dot.(tangent2_new.dot.trsp(defgrad_new)))+(id2s.bartimes.piola2_new)
 state_new=state_old

end subroutine
! *******************************************************************
! *** subroutine plasticity
! *******************************************************************
subroutine plasticity(nprop,nstat,props,time_old,time_new,defgrad_old,defgrad_new,piola_old,piola_new,state_old,state_new,internalenergy_old,internalenergy_new,dissipatedenergy_old,dissipatedenergy_new,modulus_new)
 use utils
 use tensor
 implicit none
! *** interface
 integer, intent(in) :: nprop,nstat
 real(8), intent(in) :: props(nprop),state_old(nstat),time_old,time_new
 real(8), intent(out) :: state_new(nstat)
 real(8), intent(in) :: internalenergy_old,dissipatedenergy_old
 real(8), intent(out) :: internalenergy_new,dissipatedenergy_new
 type(Tensor2), intent(in) :: defgrad_old,defgrad_new,piola_old
 type(Tensor2), intent(out) :: piola_new
 type(Tensor4), intent(out) :: modulus_new
! *** internal
 type(Tensor2s) :: cauchy_old,cauchy_new
 type(Tensor2s) :: strain_old,strain_new
 type(Tensor2s) :: elastic_strain_old,elastic_strain_new
 type(Tensor2s) :: plastic_strain_old,plastic_strain_new
 type(Tensor2s) :: stress_old,stress_new,piola2_new
 type(Tensor2s) :: flow_old,flow_new
 real(8) :: cumulated_old,cumulated_new
 real(8) :: yield_old,yield_new
 real(8) :: equivalent_old,equivalent_new
 real(8) :: rate_old,rate_new
 type(Tensor4ss) :: derivative_new,tangent2_new,auxiliary_new,tangent_new,stiffness
 real(8) :: bulk,shear,threshold,hardening,coefficient
 integer :: ij,iter
 real(8), parameter :: toler=1.d-2,small=1.d-8
 integer, parameter :: iter_max=100

 if (nstat.lt.7.or.nprop.lt.5) then
  write(*,'("ERROR: INSUFFICENT NUMBER OF STATE VARIABLES AND/OR MATERIAL PROPERTIES")')
  stop
 endif
! *** recover material properties
 bulk=props(1)/(3.d0-6.d0*props(2))
 shear=0.5d0*props(1)/(1.d0+props(2))
 threshold=props(3)
 hardening=props(4)
 coefficient=props(5)
 stiffness=((3.d0*bulk)*ps4ss)+((2.d0*shear)*pd4ss)
! *** calculate initial Cauchy-Green tensor
 cauchy_old=sym(trsp(defgrad_old).dot.defgrad_old)
! *** calculate final Cauchy-Green tensor
 cauchy_new=sym(trsp(defgrad_new).dot.defgrad_new)
! *** calculate initial strain tensor
 call isotropic(func="log",TA=cauchy_old,TC=strain_old)
 strain_old=0.5d0*strain_old
! *** calculate final strain tensor
 call isotropic(func="log",TA=cauchy_new,TC=strain_new)
 strain_new=0.5d0*strain_new
! *** recover initial plastic strain tensor
 do ij=1,6
  plastic_strain_old%val(ij)=state_old(ij)
 enddo
! *** recover initial cumulated plastic strain tensor
 cumulated_old=max(state_old(7),small)
! *** calculate initial stress tensor
 elastic_strain_old=strain_old-plastic_strain_old
 stress_old=((3.d0*bulk)*sph(elastic_strain_old))+((2.d0*shear)*dev(elastic_strain_old))
! *** calculate trial stress tensor
 plastic_strain_new=plastic_strain_old
 elastic_strain_new=strain_new-plastic_strain_new
 stress_new=((3.d0*bulk)*sph(elastic_strain_new))+((2.d0*shear)*dev(elastic_strain_new))
! *** calculate trial equivalent stress
 equivalent_new=max(sqrt(1.5d0*(dev(stress_new).ddot.dev(stress_new))),small)
! *** calculate initial yield stress
 yield_old=threshold+hardening*cumulated_old**coefficient
! *** calculate trial yield stress
 cumulated_new=cumulated_old
 yield_new=threshold+hardening*cumulated_new**coefficient
! *** calculate trial flow direction
 flow_new=(1.5d0/equivalent_new)*dev(stress_new)
! *** checks whether plastic flow will occur or not
 if ((equivalent_new-yield_new).gt.toler) then
  iter=0
  do while ((equivalent_new-yield_new).gt.toler)
   iter=iter+1
   rate_new=3.d0*shear+hardening*coefficient*cumulated_new**(coefficient-1.d0)
   cumulated_new=cumulated_new+(equivalent_new-yield_new)/rate_new
   plastic_strain_new=plastic_strain_old+((cumulated_new-cumulated_old)*flow_new)
   elastic_strain_new=strain_new-plastic_strain_new
   stress_new=((3.d0*bulk)*sph(elastic_strain_new))+((2.d0*shear)*dev(elastic_strain_new))
   equivalent_new=max(sqrt(1.5d0*(dev(stress_new).ddot.dev(stress_new))),small)
   yield_new=threshold+hardening*cumulated_new**coefficient
   flow_new=(1.5d0/equivalent_new)*dev(stress_new)
!    write(*,'("iter",i4,2es12.4)') iter,(equivalent_new-yield_new)
   if (iter.gt.iter_max) then
    write(*,'("ERROR: EXCESSIVE NUMBER OF ITERATIONS IN plasticity")')
    stop
   endif
  enddo
!   stop
  tangent_new=stiffness-((((4.d0*shear*shear)/rate_new)*flow_new)*flow_new)
 else
  tangent_new=stiffness
 endif
!  write(*,'("strain_new",6es12.4)') strain_new%val(:)
!  write(*,'("plastic_strain_new",6es12.4)') plastic_strain_new%val(:)
!  write(*,'("stress_new",6es12.4)') stress_new%val(:)
!  write(*,'("cumulated_new",6es12.4)') cumulated_new
!  write(*,'(20es12.4)') strain_new%val(:),
! *** calculate final dissipated energy density
!  dissipatedenergy_new=dissipatedenergy_old+0.5d0*((stress_old+stress_new).ddot.(plastic_strain_new-plastic_strain_old))-(0.5d0*(yield_old+yield_new)-threshold)*(cumulated_new-cumulated_old)
 dissipatedenergy_new=dissipatedenergy_old+0.5d0*((stress_old+stress_new).ddot.(plastic_strain_new-plastic_strain_old))
!                                           -0.5d0*((yield_old-threshold)*cumulated_old**(coefficient-1.d0)+(yield_new-threshold)*cumulated_new**(coefficient-1.d0))*(cumulated_new-cumulated_old)
! *** calculate final internal energy density
!  internalenergy_new=0.5d0*(stress_new.ddot.elastic_strain_new)+hardening/(coefficient+1.d0)*cumulated_new**(coefficient+1.d0)
 internalenergy_new=0.5d0*(stress_new.ddot.elastic_strain_new)!+0.5d0*hardening*cumulated_new**(2.d0*coefficient)
!  write(*,'("energy_new",6es12.4)') dissipatedenergy_new,internalenergy_new
! *** calculate PK2 stress tensor
 call isotropic(func="log",TA=cauchy_new,TB=stress_new,dTCdTA=derivative_new,TBddTCddTA=auxiliary_new)
 piola2_new=stress_new.ddot.derivative_new
! *** calculate PK1 stress tensor
 piola_new=defgrad_new.dot.piola2_new
! *** calculate tangent stress tensor
 tangent2_new=(trsp(derivative_new).ddot.(tangent_new.ddot.derivative_new))+(2.d0*auxiliary_new)
 modulus_new=rtrsp(defgrad_new.dot.(tangent2_new.dot.trsp(defgrad_new)))+(id2s.bartimes.piola2_new)
! *** assign final plastic strain tensor
 do ij=1,6
  state_new(ij)=plastic_strain_new%val(ij)
 enddo
! *** assign final cumulated plastic strain
 state_new(7)=cumulated_new

end subroutine
! *******************************************************************
! *** subroutine thermoplasticity
! *******************************************************************
subroutine thermoplasticity(nprop,nstat,props,time_old,time_new,defgrad_old,defgrad_new,piola_old,piola_new,state_old,state_new &
                           ,temp_old,temp_new,tempgrad_old,tempgrad_new,heatflow_old,heatflow_new &
                           ,internalenergy_old,internalenergy_new,dissipatedenergy_old,dissipatedenergy_new &
                           ,modulus_new,capacity_new,source_new,conductivity_new)
  use utils
 use tensor
 implicit none
! *** interface
 integer, intent(in) :: nprop,nstat
 real(8), intent(in) :: props(nprop),state_old(nstat),time_old,time_new
 real(8), intent(out) :: state_new(nstat)
 real(8), intent(in) :: internalenergy_old,dissipatedenergy_old
 real(8), intent(out) :: internalenergy_new,dissipatedenergy_new
 type(Tensor2), intent(in) :: defgrad_old,defgrad_new,piola_old
 type(Tensor2), intent(out) :: piola_new
 type(Tensor4), intent(out) :: modulus_new
 real(8), intent(in) :: temp_old,temp_new
 type(Tensor1), intent(in) :: tempgrad_old,tempgrad_new,heatflow_old
 type(Tensor1), intent(out) :: heatflow_new
 type(Tensor2), intent(out) :: conductivity_new
 real(8), intent(out) :: capacity_new,source_new
! *** internal
 type(Tensor2s) :: cauchy_old,cauchy_new
 type(Tensor2s) :: strain_old,strain_new
 type(Tensor2s) :: elastic_strain_old,elastic_strain_new
 type(Tensor2s) :: plastic_strain_old,plastic_strain_new
 type(Tensor2s) :: stress_old,stress_new,piola2_new
 type(Tensor2s) :: flow_old,flow_new
 real(8) :: cumulated_old,cumulated_new
 real(8) :: yield_old,yield_new
 real(8) :: equivalent_old,equivalent_new
 real(8) :: rate_old,rate_new
 type(Tensor4ss) :: derivative_new,tangent2_new,auxiliary_new,tangent_new,stiffness
 real(8) :: bulk,shear,threshold,hardening,coefficient,crit,density,capac,conduc
 integer :: ij,iter
 real(8), parameter :: toler=1.d-2,small=1.d-8
 integer, parameter :: iter_max=100

 if (nstat.lt.7.or.nprop.lt.9) then
  write(*,'("ERROR: INSUFFICENT NUMBER OF STATE VARIABLES AND/OR MATERIAL PROPERTIES")')
  stop
 endif
! *** recover material properties
 bulk=props(1)/(3.d0-6.d0*props(2))
 shear=0.5d0*props(1)/(1.d0+props(2))
 threshold=props(3)
 hardening=props(4)
 coefficient=props(5)
 crit=props(6)
 density=props(7)
 capac=props(8)
 conduc=props(9)
 stiffness=((3.d0*bulk)*ps4ss)+((2.d0*shear)*pd4ss)
! *** calculate initial Cauchy-Green tensor
 cauchy_old=sym(trsp(defgrad_old).dot.defgrad_old)
! *** calculate final Cauchy-Green tensor
 cauchy_new=sym(trsp(defgrad_new).dot.defgrad_new)
! *** calculate initial strain tensor
 call isotropic(func="log",TA=cauchy_old,TC=strain_old)
 strain_old=0.5d0*strain_old
! *** calculate final strain tensor
 call isotropic(func="log",TA=cauchy_new,TC=strain_new)
 strain_new=0.5d0*strain_new
! *** recover initial plastic strain tensor
 do ij=1,6
  plastic_strain_old%val(ij)=state_old(ij)
 enddo
! *** recover initial cumulated plastic strain tensor
 cumulated_old=max(state_old(7),small)
! *** calculate initial stress tensor
 elastic_strain_old=strain_old-plastic_strain_old
 stress_old=((3.d0*bulk)*sph(elastic_strain_old))+((2.d0*shear)*dev(elastic_strain_old))
! *** calculate trial stress tensor
 plastic_strain_new=plastic_strain_old
 elastic_strain_new=strain_new-plastic_strain_new
 stress_new=((3.d0*bulk)*sph(elastic_strain_new))+((2.d0*shear)*dev(elastic_strain_new))
! *** calculate trial equivalent stress
 equivalent_new=max(sqrt(1.5d0*(dev(stress_new).ddot.dev(stress_new))),small)
! *** calculate initial yield stress
 yield_old=threshold+hardening*cumulated_old**coefficient
! *** calculate trial yield stress
 cumulated_new=cumulated_old
 yield_new=threshold+hardening*cumulated_new**coefficient
! *** calculate trial flow direction
 flow_new=(1.5d0/equivalent_new)*dev(stress_new)
! *** checks whether plastic flow will occur or not
 if ((equivalent_new-yield_new).gt.toler) then
  iter=0
  do while ((equivalent_new-yield_new).gt.toler)
   iter=iter+1
   rate_new=3.d0*shear+hardening*coefficient*cumulated_new**(coefficient-1.d0)
   cumulated_new=cumulated_new+(equivalent_new-yield_new)/rate_new
   plastic_strain_new=plastic_strain_old+((cumulated_new-cumulated_old)*flow_new)
   elastic_strain_new=strain_new-plastic_strain_new
   stress_new=((3.d0*bulk)*sph(elastic_strain_new))+((2.d0*shear)*dev(elastic_strain_new))
   equivalent_new=max(sqrt(1.5d0*(dev(stress_new).ddot.dev(stress_new))),small)
   yield_new=threshold+hardening*cumulated_new**coefficient
   flow_new=(1.5d0/equivalent_new)*dev(stress_new)
!    write(*,'("iter",i4,2es12.4)') iter,(equivalent_new-yield_new)
   if (iter.gt.iter_max) then
    write(*,'("ERROR: EXCESSIVE NUMBER OF ITERATIONS IN plasticity")')
    stop
   endif
  enddo
!   stop
  tangent_new=stiffness-((((4.d0*shear*shear)/rate_new)*flow_new)*flow_new)
 else
  tangent_new=stiffness
 endif
!  write(*,'("strain_new",6es12.4)') strain_new%val(:)
!  write(*,'("plastic_strain_new",6es12.4)') plastic_strain_new%val(:)
!  write(*,'("stress_new",6es12.4)') stress_new%val(:)
!  write(*,'("cumulated_new",6es12.4)') cumulated_new
!  write(*,'(20es12.4)') strain_new%val(:),
! *** calculate final dissipated energy density
!  dissipatedenergy_new=dissipatedenergy_old+0.5d0*((stress_old+stress_new).ddot.(plastic_strain_new-plastic_strain_old))-(0.5d0*(yield_old+yield_new)-threshold)*(cumulated_new-cumulated_old)
 dissipatedenergy_new=dissipatedenergy_old+0.5d0*((stress_old+stress_new).ddot.(plastic_strain_new-plastic_strain_old))-hardening*((cumulated_new*crit)/(cumulated_new+crit))**coefficient**(crit/(cumulated_new+crit))**2.d0*(cumulated_new-cumulated_old)
! *** calculate final internal energy density
!  internalenergy_new=0.5d0*(stress_new.ddot.elastic_strain_new)+hardening/(coefficient+1.d0)*cumulated_new**(coefficient+1.d0)
 internalenergy_new=0.5d0*(stress_new.ddot.elastic_strain_new)+hardening/(coefficient+1.d0)*((cumulated_new*crit)/(cumulated_new+crit))**(coefficient+1.d0)
!  write(*,'("energy_new",6es12.4)') dissipatedenergy_new,internalenergy_new
! *** calculate PK2 stress tensor
 call isotropic(func="log",TA=cauchy_new,TB=stress_new,dTCdTA=derivative_new,TBddTCddTA=auxiliary_new)
 piola2_new=stress_new.ddot.derivative_new
! *** calculate PK1 stress tensor
 piola_new=defgrad_new.dot.piola2_new
! *** calculate tangent stress tensor
 tangent2_new=(trsp(derivative_new).ddot.(tangent_new.ddot.derivative_new))+(2.d0*auxiliary_new)
 modulus_new=rtrsp(defgrad_new.dot.(tangent2_new.dot.trsp(defgrad_new)))+(id2s.bartimes.piola2_new)
! *** assign final plastic strain tensor
 do ij=1,6
  state_new(ij)=plastic_strain_new%val(ij)
 enddo
! *** assign final cumulated plastic strain
 state_new(7)=cumulated_new
! *** calculate heat capacity per volume unit
 capacity_new=density*capac
! *** calculate thermal conductivity tensor and heat flow surface density
 conductivity_new=conduc*id2
 heatflow_new=(-conduc)*tempgrad_new
! *** calculate heat source density
 source_new=(dissipatedenergy_new-dissipatedenergy_old)/(time_new-time_old)

end subroutine
! *******************************************************************
! *** subroutine thermoelasticity
! *******************************************************************
subroutine thermoelasticity(nprop,nstat,props,time_old,time_new,defgrad_old,defgrad_new,piola_old,piola_new,state_old,state_new &
                           ,temp_old,temp_new,tempgrad_old,tempgrad_new,heatflow_old,heatflow_new &
                           ,internalenergy_old,internalenergy_new,dissipatedenergy_old,dissipatedenergy_new &
                           ,modulus_new,capacity_new,source_new,conductivity_new)
  use utils
 use tensor
 implicit none
! *** interface
 integer, intent(in) :: nprop,nstat
 real(8), intent(in) :: props(nprop),state_old(nstat),time_old,time_new
 real(8), intent(out) :: state_new(nstat)
 real(8), intent(in) :: internalenergy_old,dissipatedenergy_old
 real(8), intent(out) :: internalenergy_new,dissipatedenergy_new
 type(Tensor2), intent(in) :: defgrad_old,defgrad_new,piola_old
 type(Tensor2), intent(out) :: piola_new
 type(Tensor4), intent(out) :: modulus_new
 real(8), intent(in) :: temp_old,temp_new
 type(Tensor1), intent(in) :: tempgrad_old,tempgrad_new,heatflow_old
 type(Tensor1), intent(out) :: heatflow_new
 type(Tensor2), intent(out) :: conductivity_new
 real(8), intent(out) :: capacity_new,source_new
! *** internal
 type(Tensor2s) :: cauchy_old,cauchy_new
 type(Tensor2s) :: strain_old,strain_new
 type(Tensor2s) :: thermal_strain_old,thermal_strain_new
 type(Tensor2s) :: stress_old,stress_new
 type(Tensor2s) :: piola2_old,piola2_new
 type(Tensor4ss) :: derivative_new,tangent2_new,auxiliary_new,tangent_new
 real(8) :: bulk,shear,density,capac,conduc,expan,ref

 if (nprop.lt.7) then
  write(*,'("ERROR: INSUFFICENT NUMBER OF MATERIAL PROPERTIES")')
  stop
 endif
 bulk=props(1)/(3.d0-6.d0*props(2))
 shear=0.5d0*props(1)/(1.d0+props(2))
 density=props(3)
 capac=props(4)
 conduc=props(5)
 expan=props(6)
 ref=props(7)
! *** calculate initial Cauchy-Green tensor
 cauchy_old=sym(trsp(defgrad_old).dot.defgrad_old)
! *** calculate final Cauchy-Green tensor
 cauchy_new=sym(trsp(defgrad_new).dot.defgrad_new)
! *** calculate initial strain tensor
 call isotropic(func="log",TA=cauchy_old,TC=strain_old)
 strain_old=0.5d0*strain_old
! *** calculate final strain tensor
 call isotropic(func="log",TA=cauchy_new,TC=strain_new)
 strain_new=0.5d0*strain_new
! *** calculate final thermal strain tensor
 thermal_strain_new=(expan*(temp_new-ref))*id2s
! *** calculate final stress tensor
 stress_new=((3.d0*bulk)*sph(strain_new-thermal_strain_new))+((2.d0*shear)*dev(strain_new))
! *** calculate final dissipated energy density
 dissipatedenergy_new=dissipatedenergy_old
! *** calculate final internal energy density
 internalenergy_new=0.5d0*bulk*trc(strain_new-thermal_strain_new)**2.d0+shear*(dev(strain_new).ddot.dev(strain_new)) &
                   +3.d0*bulk*trc(strain_new-thermal_strain_new)*expan*temp_new+density*capac*(temp_new-ref)
! *** calculate PK2 stress tensor
 call isotropic(func="log",TA=cauchy_new,TB=stress_new,dTCdTA=derivative_new,TBddTCddTA=auxiliary_new)
 piola2_new=stress_new.ddot.derivative_new
! *** calculate PK1 stress tensor
 piola_new=defgrad_new.dot.piola2_new
! *** calculate tangent stress tensor
 tangent_new=((3.d0*bulk)*ps4ss)+((2.d0*shear)*pd4ss)
 tangent2_new=(trsp(derivative_new).ddot.(tangent_new.ddot.derivative_new))+(2.d0*auxiliary_new)
 modulus_new=rtrsp(defgrad_new.dot.(tangent2_new.dot.trsp(defgrad_new)))+(id2s.bartimes.piola2_new)
 state_new=state_old
! *** calculate heat capacity per volume unit
 capacity_new=density*capac-9.d0*bulk*expan**2.d0*temp_new
! *** calculate thermal conductivity tensor and heat flow surface density
 conductivity_new=conduc*id2
 heatflow_new=(-conduc)*tempgrad_new
! *** calculate heat source density
 source_new=-3.d0*bulk*temp_new*expan*(trc(strain_new)-trc(strain_old))/(time_new-time_old)

end subroutine
