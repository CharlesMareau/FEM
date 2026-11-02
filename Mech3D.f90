program fem_main
 use utils
 use tensor
 use iso_fortran_env
 implicit none
 integer :: freq
! Variables maillage
 integer :: nnode,nelem,ndof,nnode_per_element,nntg,nntg_per_element,npair
 real(8), allocatable :: coords(:,:),location(:,:)
 integer, allocatable :: conn(:,:)
! Données matériaux
 character(5) :: cname
 character(len=64), allocatable :: nset_names(:)
 integer, allocatable :: nset_sizes(:),nset_nodes(:,:)
! Variables pour assemblage
 integer :: i,j,k,l,ij,kl,idof,jdof,kdof,iddl,jddl,ielem,inode,jnode,iter,istp,unit,istat,intg,jntg,igss,jgss,iset,nstat,nprop
 integer, allocatable :: bc_type(:)
 real(8), allocatable :: bc_val(:)
 real(8), allocatable :: stiffness(:),force(:),increment(:),derivative(:,:,:),volume(:),wght(:)
 real(8), allocatable :: position_new(:),position_old(:)
 type(Tensor2), allocatable :: defgrad_new(:),defgrad_old(:)
 type(Tensor2), allocatable :: stress_new(:),stress_old(:)
 type(Tensor4), allocatable :: modulus_new(:),modulus_old(:)
 real(8), allocatable :: state_new(:,:),state_old(:,:)
 real(8), allocatable :: internalenergy_new(:),internalenergy_old(:)
 real(8), allocatable :: dissipatedenergy_new(:),dissipatedenergy_old(:)
 real(8), allocatable :: x(:),y(:),z(:),matrix(:,:),vector(:),auxiliary(:,:),props(:)
 real(8) :: dtime_min,dtime_max
 real(8) :: time_old,time_new,dtime
 real(8) :: duration,progress,toler=1.d-2,delta
 integer :: nnz,mtype,nrhs
 integer, allocatable :: ia(:), ja(:),correspond(:,:),permutation(:)
 integer :: iparm(64), pt(64)
 integer :: maxfct, mnum, phase, error, msglvl
! *** variables for data output
 character(40) :: filevtu
 character(8) :: fil
 integer :: ifile
 real(8), parameter :: small=1.d-5
 real(8) :: residual_old,residual,damping,res_max
 integer :: idof_max
 character(80) :: filemsh,fileproc,filemat,prosa

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
 call read_inp(filemsh,coords,conn,nnode,nelem,nset_names,nset_sizes,nset_nodes)
!  call read_msh("hole.msh", coords, conn, nnode, nelem)
 write(*,'("Number of elements",i6)') nelem
 write(*,'("Number of nodes",i6)') nnode
 ndof=3*nnode
 write(*,'("Number of DOFs",i6)') ndof
 nntg=nntg_per_element*nelem
 write(*,'("Number of integration points",i6)') nntg

 allocate(x(nnode_per_element),y(nnode_per_element),z(nnode_per_element))
 allocate(matrix(3*nnode_per_element,3*nnode_per_element),vector(3*nnode_per_element),auxiliary(9,3*nnode_per_element))
! *** initializes variables
 allocate(position_new(ndof),position_old(ndof))
 allocate(defgrad_new(nntg),defgrad_old(nntg))
 allocate(stress_new(nntg),stress_old(nntg))
 allocate(modulus_new(nntg),modulus_old(nntg))
 allocate(state_new(nstat,nntg),state_old(nstat,nntg))
 allocate(dissipatedenergy_new(nntg),dissipatedenergy_old(nntg))
 allocate(internalenergy_new(nntg),internalenergy_old(nntg))
 do ielem=1,nelem
  do inode=1,nnode_per_element
   idof=3*conn(inode,ielem)-2
   jdof=3*conn(inode,ielem)-1
   kdof=3*conn(inode,ielem)
   position_old(idof)=coords(1,conn(inode,ielem))
   position_old(jdof)=coords(2,conn(inode,ielem))
   position_old(kdof)=coords(3,conn(inode,ielem))
   position_new(idof)=position_old(idof)
   position_new(jdof)=position_old(jdof)
   position_new(kdof)=position_old(kdof)
  enddo
  do intg=1,nntg_per_element
   igss=(ielem-1)*nntg_per_element+intg
   internalenergy_old(igss)=0.d0
   internalenergy_new(igss)=0.d0
   dissipatedenergy_old(igss)=0.d0
   dissipatedenergy_new(igss)=0.d0
   stress_old(igss)=null2
   stress_new(igss)=null2
   defgrad_old(igss)=id2
   defgrad_new(igss)=id2
   do istat=1,nstat
    state_old(istat,igss)=0.d0
    state_new(istat,igss)=0.d0
   enddo
  enddo
 enddo
! *** calculates the shape function derivative
 allocate(derivative(9,3*nnode_per_element,nntg),location(3,nntg),wght(nntg),volume(nelem))
 do ielem=1,nelem
  do inode=1,nnode_per_element
   x(inode)=coords(1,conn(inode,ielem))
   y(inode)=coords(2,conn(inode,ielem))
   z(inode)=coords(3,conn(inode,ielem))
  enddo
  igss=(ielem-1)*nntg_per_element+1
  jgss=(ielem-1)*nntg_per_element+nntg_per_element
  select case(nntg_per_element)
   case(1)
    call tetrahedronT4(x,y,z,derivative(:,:,igss:jgss),location(:,igss:jgss),wght(igss:jgss),volume(ielem))
   case(4)
    call tetrahedronT10(x,y,z,derivative(:,:,igss:jgss),location(:,igss:jgss),wght(igss:jgss),volume(ielem))
   case(8)
    call hexahedronH8(x,y,z,derivative(:,:,igss:jgss),location(:,igss:jgss),wght(igss:jgss),volume(ielem))
   case default
    write(*,'("Unknown type of element. Program will stop !")')
    stop
  end select
 enddo
! *** generates the csr structure for storage
 call build_csr_structure
 stiffness=0.0_8
 force=0.0_8
 increment=0.0_8
! *** sets the parameters for the PARDISO solver
 mtype = 11                ! matrice réelle symétrique indéfinie : 11 / non symétrique : 11
 nrhs  = 1
 error = 0
 msglvl = 0
 maxfct = 1
 mnum = 1
 pt = 0
 iparm = 0
 iparm(1) = 1              ! utiliser iparm personnalisé
 iparm(3) = 0              ! pas de message
 iparm(7) = 1              ! nombre de threads
 iparm(10) = 13            ! pivot perturbation
 iparm(17) = -1            ! pas de réordonnancement METIS
 iparm(34) = 0             ! utiliser zero-based indexing → =1 pour C, =0 pour Fortran
! *** finds nodes with prescribed BC
 allocate(bc_type(ndof),bc_val(ndof))
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
  write(*,'("ISTEP:",i6," TIME:",es12.4)') istp,time_old
  dtime=min(duration-time_old,dtime)
  time_new=time_old+dtime
  progress=progress+dtime/duration
! *** applies boundary condition
  do idof=1,ndof
   if (bc_type(idof).eq.1) then
    position_new(idof)=position_old(idof)+bc_val(idof)*dtime/duration
   endif
  enddo
  residual=toler+1.d0
  damping=1.d0
  iter=0
  do while (residual.gt.toler)
   iter=iter+1
   if (iter.gt.2000) then
    write(*,'("MAX NUMBER OF ITERATIONS HAS BEEN REACHED.")')
    call export_vtu("fail.vtu")
    stop
   endif
! *** calculates the deformation gradient tensor
   do ielem=1,nelem
    do intg=1,nntg_per_element
     igss=(ielem-1)*nntg_per_element+intg
     defgrad_new(igss)=null2
     do inode=1,nnode_per_element
      do iddl=-2,0
       idof=3*conn(inode,ielem)+iddl
       defgrad_new(igss)%val(:,:)=defgrad_new(igss)%val(:,:)+reshape(derivative(:,3*inode+iddl,igss),(/3,3/))*position_new(idof)
      enddo
     enddo
    enddo
   enddo
! *** integrates constitutive relations
   select case(cname)
    case("ELAST","elast","Elast")
     !$omp parallel default(shared)
     !$omp do private(intg) schedule(static)
     do intg=1,nntg
      call elasticity(nprop,nstat,props,time_old,time_new,defgrad_old(intg),defgrad_new(intg),stress_old(intg),stress_new(intg),state_old(:,intg),state_new(:,intg) &
                   ,internalenergy_old(intg),internalenergy_new(intg),dissipatedenergy_old(intg),dissipatedenergy_new(intg),modulus_new(intg))
     enddo
     !$omp enddo
     !$omp end parallel
    case("PLAST","plast","Plast")
     !$omp parallel default(shared)
     !$omp do private(intg) schedule(static)
     do intg=1,nntg
      call plasticity(nprop,nstat,props,time_old,time_new,defgrad_old(intg),defgrad_new(intg),stress_old(intg),stress_new(intg),state_old(:,intg),state_new(:,intg) &
                   ,internalenergy_old(intg),internalenergy_new(intg),dissipatedenergy_old(intg),dissipatedenergy_new(intg),modulus_new(intg))
     enddo
     !$omp enddo
     !$omp end parallel
    case default
     write(*,'("ERROR: UNKNOWN MATERIAL TYPE")')
     stop
   end select
! *** evaluates stiffness matrix
  stiffness=0.d0
  !$omp parallel default(shared)
  !$omp do schedule(static) private(ielem,intg,igss,matrix,inode,iddl,idof,jnode,jddl,jdof,kdof)
  do ielem=1,nelem
   do intg=1,nntg_per_element
    igss=(ielem-1)*nntg_per_element+intg
     matrix=matmul(transpose(derivative(:,:,igss)),matmul(reshape(modulus_new(igss)%val(:,:,:,:),(/9,9/)),derivative(:,:,igss)))*wght(igss)
     do inode=1,nnode_per_element
      do iddl=-2,0
       idof=3*conn(inode,ielem)+iddl
       do jnode=1,nnode_per_element
        do jddl=-2,0
         jdof=3*conn(jnode,ielem)+jddl
         kdof=correspond(idof,jdof)
         !$omp atomic update
         stiffness(kdof)=stiffness(kdof)+matrix(3*inode+iddl,3*jnode+jddl)
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo
   !$omp end do
   !$omp end parallel
! *** evaluates force vector
   force=0.d0
  !$omp parallel default(shared)
  !$omp do schedule(static) private(ielem,intg,igss,vector,inode,iddl,idof)
   do ielem=1,nelem
    do intg=1,nntg_per_element
     igss=(ielem-1)*nntg_per_element+intg
     vector=matmul(transpose(derivative(:,:,igss)),reshape(stress_new(igss)%val(:,:),(/9/)))*wght(igss)
     do inode=1,nnode_per_element
      do iddl=-2,0
       idof=3*conn(inode,ielem)+iddl
       !$omp atomic update
       force(idof)=force(idof)+vector(3*inode+iddl)
      enddo
     enddo
    enddo
   enddo
   !$omp end do
   !$omp end parallel
! *** prescribes boundary conditions
   do idof=1,ndof
    if (bc_type(idof).eq.1) then
     force(idof)=0.d0
     do jdof=1,ndof
      stiffness(correspond(idof,jdof))=0.d0
     enddo
     stiffness(correspond(idof,idof))=1.d0
    endif
    if (bc_type(idof).eq.2) then
     force(idof)=force(idof)-bc_val(idof)*progress
    endif
   enddo
! *** calculates the residual
   residual_old=residual
   residual=0.d0
   res_max=0.d0
   do idof=1,ndof
    if (bc_type(idof).ne.1) then
     residual=residual+force(idof)**2.d0
    endif
   enddo
   residual=sqrt(residual)
   if (residual.gt.residual_old) then
    damping=0.9d0*damping
   endif
   write(*,'("  ITER:",i4," RESIDUAL:",es12.4)') iter,residual
!    pause
! *** uses Newton-Raphson algorithm to update solution with PARDISO solver
   if (istp.eq.1.and.iter.eq.1) then
    phase = 11 ! *** analysis
    call pardiso(pt,maxfct,mnum,mtype,phase,ndof,stiffness,ia,ja,permutation,nrhs,iparm,msglvl,force,increment,error)
   endif
   phase=23 ! *** factorization + solution
   call pardiso(pt,maxfct,mnum,mtype,phase,ndof,stiffness,ia,ja,permutation,nrhs,iparm,msglvl,force,increment,error)
   do idof=1,ndof
    position_new(idof)=position_new(idof)-increment(idof)*damping
   enddo
  enddo
! *** updates variables
  time_old=time_new
  do idof=1,ndof
   delta=position_new(idof)-position_old(idof)
   position_old(idof)=position_new(idof)
   position_new(idof)=position_old(idof)+delta
  enddo
  do ielem=1,nelem
   do intg=1,nntg_per_element
    igss=(ielem-1)*nntg_per_element+intg
    internalenergy_old(igss)=internalenergy_new(igss)
    dissipatedenergy_old(igss)=dissipatedenergy_new(igss)
    stress_old(igss)=stress_new(igss)
    defgrad_old(igss)=defgrad_new(igss)
    do istat=1,nstat
     delta=state_new(istat,igss)-state_old(istat,igss)
     state_old(istat,igss)=state_new(istat,igss)
     state_new(istat,igss)=state_old(istat,igss)+delta
    enddo
   enddo
  enddo
  if (freq.ne.0) then
   if (modulo(istp,freq).eq.0) then
    ifile=ifile+1
    write(fil,'(i8)') ifile
    fil=adjustl(fil)
    filevtu="result"//trim(fil)//".vtu"
    call export_vtu(filevtu)
    call export_rct(1)
   endif
  endif
 enddo
! *** export final results
 call export_vtu("result.vtu")
 call export_rct(1)
 call export_rct(2)
! *** cleans memory for PARDISO solver
 phase = 0
 call pardiso(pt,maxfct,mnum,mtype,phase,ndof,stiffness,ia,ja,permutation,nrhs,iparm,msglvl,force,increment,error)

contains
! *******************************************************************
! *** subroutine build_csr_structure
! *******************************************************************
subroutine build_csr_structure
 implicit none
 integer, allocatable :: tag(:,:)
 integer :: cnt
 allocate(tag(ndof,ndof))
 tag=0
 do ielem=1,nelem
  do inode=1,nnode_per_element
   do iddl=-2,0
    idof=3*conn(inode,ielem)+iddl
    do jnode=1,nnode_per_element
     do jddl=-2,0
      jdof=3*conn(jnode,ielem)+jddl
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
 allocate(stiffness(nnz),force(ndof),increment(ndof))
 allocate(ja(nnz),ia(ndof+1),correspond(ndof,ndof),permutation(ndof))
 cnt=1
 ia(1)=1
 do idof=1,ndof
  do jdof=1,ndof
   if (tag(idof,jdof).ne.0) then
    correspond(idof,jdof)=cnt
    ja(cnt)=jdof
    cnt=cnt+1
   endif
  enddo
  ia(idof+1)=cnt
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
! *** export_rct
! *******************************************************************
subroutine export_rct(iopt)
 implicit none
 integer :: iset,iopt
 real(8) :: reaction(3),moment(3)
 select case(iopt)
 case(0)
  open(unit=20, file="reaction.out", status="replace", action="write")
  write(20,'(" TIME")',advance='no')
  do iset=1,size(nset_names)
   if (iset.eq.size(nset_names)) then
    write(20,'(" RF1_",a," RF2_",a," RF3_",a," M1_",a," M2_",a," M3_",a)',advance='yes') trim(nset_names(iset)),trim(nset_names(iset)),trim(nset_names(iset)),trim(nset_names(iset)),trim(nset_names(iset)),trim(nset_names(iset))
   else
    write(20,'(" RF1_",a," RF2_",a," RF3_",a," M1_",a," M2_",a," M3_",a)',advance='no') trim(nset_names(iset)),trim(nset_names(iset)),trim(nset_names(iset)),trim(nset_names(iset)),trim(nset_names(iset)),trim(nset_names(iset))
   endif
  enddo
 case(1)
  write(20,'(1es12.4)',advance='no') time_old
  force=0.d0
  do ielem=1,nelem
   do intg=1,nntg_per_element
    igss=(ielem-1)*nntg_per_element+intg
    vector=matmul(transpose(derivative(:,:,igss)),reshape(stress_old(igss)%val(:,:),(/9/)))*wght(igss)
    do inode=1,nnode_per_element
     do iddl=-2,0
      idof=3*conn(inode,ielem)+iddl
      force(idof)=force(idof)+vector(3*inode+iddl)
     enddo
    enddo
   enddo
  enddo
  do idof=1,ndof
   if (bc_type(idof).eq.2) then
    force(idof)=force(idof)-bc_val(idof)*progress
   endif
  enddo
  do iset=1,size(nset_names)
!    write(*,'(a)') nset_names(iset)
   reaction(1)=0.d0
   reaction(2)=0.d0
   reaction(3)=0.d0
   moment(1)=0.d0
   moment(2)=0.d0
   moment(3)=0.d0
   do inode=1,nset_sizes(iset)
    jnode=nset_nodes(iset,inode)
    idof=3*jnode-2
    jdof=3*jnode-1
    kdof=3*jnode
    if (bc_type(idof).eq.1) then
!      write(*,'("Node",i6," DOF:",i6," Reaction force:",es12.4)') jnode,idof,force(idof)
     reaction(1)=reaction(1)+force(idof)
    endif
    if (bc_type(jdof).eq.1) then
!      write(*,'("Node",i6," DOF:",i6," Reaction force:",es12.4)') jnode,idof,force(idof)
     reaction(2)=reaction(2)+force(jdof)
    endif
    if (bc_type(jdof).eq.1) then
!      write(*,'("Node",i6," DOF:",i6," Reaction force:",es12.4)') jnode,idof,force(idof)
     reaction(3)=reaction(3)+force(kdof)
    endif
    moment(1)=moment(1)+force(kdof)*position_new(jdof)-force(jdof)*position_new(kdof)
    moment(2)=moment(2)+force(idof)*position_new(kdof)-force(kdof)*position_new(idof)
    moment(3)=moment(3)+force(jdof)*position_new(idof)-force(idof)*position_new(jdof)
   enddo
   if (iset.eq.size(nset_names)) then
    write(20,'(6es12.4)',advance='yes') reaction(:),moment(:)
   else
    write(20,'(6es12.4)',advance='no') reaction(:),moment(:)
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
 real(8) :: stress(3,3),defgrad(3,3),state(nstat),energy
 integer :: vtk_type
! *** check element type: C3D4 or C3D8
 select case(nnode_per_element)
 case(4)
  if (nntg_per_element.ne.1) then
   write(*,*) 'Error export_vtu: C3D4 requires nntg_per_element=1'
   stop
  endif
  vtk_type=10
 case(8)
  if (nntg_per_element.ne.8) then
   write(*,*) 'Error export_vtu: C3D8 requires nntg_per_element=8'
   stop
  endif
  vtk_type=12
 case(10)
  if (nntg_per_element.ne.4) then
   write(*,*) 'Error export_vtu: C3D10 requires nntg_per_element=4'
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
  write(unit,'(a,I0,a,I0,a)') '    <Piece NumberOfPoints="',nnode,'" NumberOfCells="',nelem,'">'
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
  do ielem=1,nelem
   do inode=1,nnode_per_element
    write(unit,'(I0,1X)',advance='no') conn(inode,ielem)-1
   enddo
   write(unit,*)
  enddo
  write(unit,'(a)') '        </DataArray>'
  write(unit,'(a)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
  do ielem=1,nelem
   write(unit,'(I0)') ielem*nnode_per_element
  enddo
  write(unit,'(a)') '        </DataArray>'
  write(unit,'(a)') '        <DataArray type="UInt8" Name="types" format="ascii">'
  do ielem=1,nelem
   write(unit,'(I0)') vtk_type
  enddo
  write(unit,'(a)') '        </DataArray>'
  write(unit,'(a)') '      </Cells>'
  write(unit,'(a)') '      <CellData Scalars="scalars">'
! *** Compute average stress tensor
  write(unit,'(a)') '        <DataArray type="Float64" Name="Stress" NumberOfComponents="9" format="ascii">'
  do ielem=1,nelem
   do i=1,3
    do j=1,3
     stress(i,j)=0.d0
     do intg=1,nntg_per_element
      igss=(ielem-1)*nntg_per_element+intg
      stress(i,j)=stress(i,j)+wght(igss)*stress_old(igss)%val(i,j)/volume(ielem)
     enddo
    enddo
   enddo
   write(unit,'(9E16.6E3)') ((stress(i,j),j=1,3),i=1,3)
  enddo
  write(unit,'(a)') '        </DataArray>'
! *** Compute average deformation gradient
  write(unit,'(a)') '        <DataArray type="Float64" Name="DeformationGradient" NumberOfComponents="9" format="ascii">'
  do ielem=1,nelem
   do i=1,3
    do j=1,3
     defgrad(i,j)=0.d0
     do intg=1,nntg_per_element
      igss=(ielem-1)*nntg_per_element+intg
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
   do ielem=1,nelem
    do istat=1,nstat
     state(istat)=0.d0
     do intg=1,nntg_per_element
      igss=(ielem-1)*nntg_per_element+intg
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
  do ielem=1,nelem
   energy=0.d0
   do intg=1,nntg_per_element
    igss=(ielem-1)*nntg_per_element+intg
    energy=energy+wght(igss)*internalenergy_old(igss)/volume(ielem)
   enddo
   write(unit,'(1E16.6E3)') energy
  enddo
  write(unit,'(a)') '        </DataArray>'
! *** Compute average dissipated energy
  write(unit,'(a)') '        <DataArray type="Float64" Name="DissipatedEnergy" NumberOfComponents="1" format="ascii">'
  do ielem=1,nelem
   energy=0.d0
   do intg=1,nntg_per_element
    igss=(ielem-1)*nntg_per_element+intg
    energy=energy+wght(igss)*dissipatedenergy_old(igss)/volume(ielem)
   enddo
   write(unit,'(1E16.6E3)') energy
  enddo
  write(unit,'(a)') '        </DataArray>'
  write(unit,'(a)') '      </CellData>'
  write(unit,'(a)') '    </Piece>'
  write(unit,'(a)') '  </UnstructuredGrid>'
  write(unit,'(a)') '</VTKFile>'
 close(unit)

end subroutine export_vtu
! *******************************************************************
! *** subroutine export_vtu
! *******************************************************************
subroutine export_vtu_old(filename)
 implicit none
 character(*), intent(in) :: filename
 real(8) :: stress(3,3),defgrad(3,3),state(nstat),energy

 open(newunit=unit, file=filename, status="replace", action="write")
  write(unit,'(a)') '<?xml version="1.0"?>'
  write(unit,'(a)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
  write(unit,'(a)') '  <UnstructuredGrid>'
  write(unit,'(a,I6,a,I6,a)') '    <Piece NumberOfPoints="', nnode, '" NumberOfCells="', nelem, '">'
  write(unit,'(a)') '      <Points>'
  write(unit,'(a)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
  do inode=1,nnode
   idof=3*inode-2
   jdof=3*inode-1
   kdof=3*inode
   write(unit,'(3E16.6E3)') position_old(idof),position_old(jdof),position_old(kdof)
!    write(unit,'(3E14.6E3)') coords(1,inode),coords(2,inode),0.0d0
  end do
  write(unit,'(a)') '        </DataArray>'
  write(unit,'(a)') '      </Points>'
  ! Connectivité
  write(unit,'(a)') '      <Cells>'
  write(unit,'(a)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
  do ielem=1,nelem
   do inode=1,(nnode_per_element-1)
     write(unit,'(I8)',advance='no') conn(inode,ielem)-1
   enddo
   write(unit,'(I8)',advance='yes') conn(inode,ielem)-1
!     write(unit,'(3I8)') conn(1,ielem)-1,conn(2,ielem)-1,conn(3,ielem)-1  ! VTK uses 0-based indexing
  end do
  write(unit,'(a)') '        </DataArray>'
  ! Offsets
  write(unit,'(a)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
  do ielem=1,nelem
    write(unit,'(I8)') ielem * nnode_per_element
  end do
  write(unit,'(a)') '        </DataArray>'
  write(unit,'(a)') '        <DataArray type="UInt8" Name="types" format="ascii">'
  do ielem=1,nelem
    write(unit,'(I8)') 12
  end do
  write(unit,'(a)') '        </DataArray>'
  write(unit,'(a)') '      </Cells>'
  write(unit,'(a)') '      <CellData Scalars="scalars">'
  write(unit,'(a)') '        <DataArray type="Float64" Name="Stress" NumberOfComponents="9" format="ascii">'
  do ielem=1,nelem
   do i=1,3
    do j=1,3
     stress(i,j)=0.d0
     do intg=1,nntg_per_element
      igss=(ielem-1)*nntg_per_element+intg
      stress(i,j)=stress(i,j)+wght(igss)*stress_old(igss)%val(i,j)/volume(ielem)
     enddo
    enddo
   enddo
   write(unit,'(9E16.6E3)') ((stress(i,j),j=1,3),i=1,3)
  end do
  write(unit,*) '        </DataArray>'
  write(unit,'(a)') '        <DataArray type="Float64" Name="DeformationGradient" NumberOfComponents="9" format="ascii">'
  do ielem=1,nelem
   do i=1,3
    do j=1,3
     defgrad(i,j)=0.d0
     do intg=1,nntg_per_element
      igss=(ielem-1)*nntg_per_element+intg
      defgrad(i,j)=defgrad(i,j)+wght(igss)*defgrad_old(igss)%val(i,j)/volume(ielem)
     enddo
    enddo
   enddo
   write(unit,'(9E16.6E3)') ((defgrad(i,j),j=1,3),i=1,3)
  end do
  write(unit,'(a)') '        </DataArray>'
  if (nstat.gt.0) then
   write(unit,'(a,I1,a)') '        <DataArray type="Float64" Name="StateVariable" NumberOfComponents="',nstat,'" format="ascii">'
   do ielem=1,nelem
    do istat=1,nstat
     state(istat)=0.d0
     do intg=1,nntg_per_element
      igss=(ielem-1)*nntg_per_element+intg
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
  write(unit,'(a)') '        <DataArray type="Float64" Name="InternalEnergy" NumberOfComponents="1" format="ascii">'
  do ielem=1,nelem
   energy=0.d0
   do intg=1,nntg_per_element
    igss=(ielem-1)*nntg_per_element+intg
    energy=energy+wght(igss)*internalenergy_old(igss)/volume(ielem)
   enddo
   write(unit,'(1E16.6E3)') energy
  enddo
  write(unit,'(a)') '        </DataArray>'
  write(unit,'(a)') '        <DataArray type="Float64" Name="DissipatedEnergy" NumberOfComponents="1" format="ascii">'
  do ielem=1,nelem
   energy=0.d0
   do intg=1,nntg_per_element
    igss=(ielem-1)*nntg_per_element+intg
    energy=energy+wght(igss)*dissipatedenergy_old(igss)/volume(ielem)
   enddo
   write(unit,'(1E16.6E3)') energy
  enddo
  write(unit,'(a)') '        </DataArray>'
  write(unit,'(a)') '      </CellData>'
  write(unit,'(a)') '    </Piece>'
  write(unit,'(a)') '  </UnstructuredGrid>'
  write(unit,'(a)') '</VTKFile>'
 close(unit)

end subroutine

! *******************************************************************
! *** subroutine tetrahedronT4
! *******************************************************************
subroutine tetrahedronT4(x,y,z,B,L,W,V)
 implicit none
 real(8), intent(in)  :: x(4),y(4),z(4)
 real(8), intent(out) :: B(9,12,1),W(1),V,L(3,1)
 real(8) :: dN(4,3)
 real(8) :: dNdx(4),dNdy(4),dNdz(4)
 real(8) :: J(3,3),detJ,invJ(3,3)
 integer :: i

 V=0.d0
 B=0.d0
 L=0.d0
 W=0.d0

 ! *** Shape function derivatives in natural coordinates
 !     N1=1-xi-eta-zeta, N2=xi, N3=eta, N4=zeta
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
subroutine tetrahedronT10(x,y,z,B,L,W,V)
 implicit none
 real(8), intent(in)  :: x(10),y(10),z(10)
 real(8), intent(out) :: B(9,30,4),W(4),V,L(3,4)
 real(8) :: N(10),dN(10,3)
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
   N(i)=lam(i)*(2.d0*lam(i)-1.d0)
  enddo
  N(5)=4.d0*lam(1)*lam(2)
  N(6)=4.d0*lam(2)*lam(3)
  N(7)=4.d0*lam(3)*lam(1)
  N(8)=4.d0*lam(1)*lam(4)
  N(9)=4.d0*lam(2)*lam(4)
  N(10)=4.d0*lam(3)*lam(4)
! *** Physical coordinates of integration point
  L(:,p)=0.d0
  do i=1,10
   L(1,p)=L(1,p)+N(i)*x(i)
   L(2,p)=L(2,p)+N(i)*y(i)
   L(3,p)=L(3,p)+N(i)*z(i)
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
subroutine hexahedronH8(x,y,z,B,L,W,V)
 implicit none
 real(8), intent(in)  :: x(8),y(8),z(8)
 real(8), intent(out) :: B(9,24,8),W(8),V,L(3,8)
 real(8) :: N(8),dN(8,3)
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
   N(i)=0.125d0*(1.d0+dble(sx(i))*xi)*(1.d0+dble(sy(i))*eta)*(1.d0+dble(sz(i))*zeta)
  enddo
 ! *** Physical coordinates of the Gauss point
  L(:,p)=0.d0
  do i=1,8
   L(1,p)=L(1,p)+N(i)*x(i)
   L(2,p)=L(2,p)+N(i)*y(i)
   L(3,p)=L(3,p)+N(i)*z(i)
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
! *** subroutine triangleT6
! *******************************************************************
subroutine triangleT6(x,y,B,L,W,A)
 implicit none
 real(8), intent(in) :: x(6),y(6)      ! Coordonnées (x,y) des 6 nœuds
 real(8), intent(out):: B(4,12,3),W(3),A,L(2,3)
 real(8) :: dNdx(6),dNdy(6),N(6),dN(6,2)
 real(8) :: J(2,2),detJ,invJ(2,2)
 real(8) :: L1,L2,L3
 integer :: i,p
 real(8), parameter :: gauss(2,3)=(/1.d0/6.d0,1.d0/6.d0,2.d0/3.d0,1.d0/6.d0,1.d0/6.d0,2.d0/3.d0/)
 real(8), parameter :: wght(3) = (/1.d0/3.d0,1d0/3d0,1d0/3.d0/)
 real(8) :: xg,yg

 A=0.d0
 do p=1,3
! *** evaluates element coordinates
  L1=1.d0-gauss(1,p)-gauss(2,p)
  L2=gauss(1,p)
  L3=gauss(2,p)
! *** evaluates shape functions
  N(1)=L1*(2d0*L1-1d0)
  N(2)=L2*(2d0*L2-1d0)
  N(3)=L3*(2d0*L3-1d0)
  N(4)=4d0*L1*L2
  N(5)=4d0*L2*L3
  N(6)=4d0*L3*L1
  L(:,p)=0.d0
  do i=1,6
   L(1,p)=L(1,p)+N(i)*x(i)
   L(2,p)=L(2,p)+N(i)*y(i)
  enddo
!   write(*,'("GAUSS POINT",2es12.4)') xg,yg
! *** evaluates shape function derivatives wrt element coordinates
  dN(1,1)=4.d0*gauss(1,p)+4.d0*gauss(2,p)-3.d0
  dN(1,2)=4.d0*gauss(2,p)+4.d0*gauss(1,p)-3.d0
  dN(2,1)=4.d0*gauss(1,p)-1.d0
  dN(2,2)=0.d0
  dN(3,1)=0.d0
  dN(3,2)=4.d0*gauss(2,p)-1.d0
  dN(4,1)=4.d0-8.d0*gauss(1,p)-4.d0*gauss(2,p)
  dN(4,2)=-4.d0*gauss(1,p)
  dN(5,1)=4.d0*gauss(2,p)
  dN(5,2)=4.d0*gauss(1,p)
  dN(6,1)=-4.d0*gauss(2,p)
  dN(6,2)=4.d0-8.d0*gauss(2,p)-4.d0*gauss(1,p)
! *** evaluates the Jacobian and its determinant
  J(:,:)=0.d0
  do i=1,6
   J(1,1)=J(1,1)+dN(i,1)*x(i)
   J(1,2)=J(1,2)+dN(i,2)*x(i)
   J(2,1)=J(2,1)+dN(i,1)*y(i)
   J(2,2)=J(2,2)+dN(i,2)*y(i)
  enddo
!   write(*,'("Jacobian",2es12.4)') J
  detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1)
  if (detJ.le.0.d0) then
   write(*,'("Negative jacobian in TriangleT6")')
   stop
  endif
! *** evaluates the Jacobian inverse
  invJ(1,1)=J(2,2)/detJ
  invJ(1,2)=-J(1,2)/detJ
  invJ(2,1)=-J(2,1)/detJ
  invJ(2,2)=J(1,1)/detJ
! *** evaluates the derivative of the wght functions
  do i = 1,6
   dNdx(i)=invJ(1,1)*dN(i,1)+invJ(2,1)*dN(i,2) ! *** to be verified
   dNdy(i)=invJ(1,2)*dN(i,1)+invJ(2,2)*dN(i,2) ! *** to be verified
  end do
  do i = 1,6
   B(1,2*i-1,p)=dNdx(i)
   B(2,2*i,p)=dNdy(i)
   B(3,2*i-1,p)=dNdy(i)
   B(4,2*i,p)=dNdx(i)
  end do
  W(p)=detJ*wght(p)*.5d0
  A=A+W(p)
 enddo

end subroutine triangleT6
! *******************************************************************
! *** subroutine triangleT3
! *******************************************************************
subroutine triangleT3(x, y, B, L, W, A)
 implicit none
 real(8), intent(in) :: x(3), y(3)  ! material coordinates of three nodes
 real(8) :: B(4,6,1), W(1), A, L(2,1)
 integer :: i
! Calcul de l'aire du triangle
 A = 0.5d0 * ((x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)))
 B(1,1,1)=(y(2)-y(3))/(2.0d0*A) ! *** dN1/dx
 B(2,1,1)=0.d0
 B(3,1,1)=(x(3)-x(2))/(2.0d0*A) ! *** dN1/dy
 B(4,1,1)=0.d0
 B(1,2,1)=0.d0
 B(2,2,1)=(x(3)-x(2))/(2.0d0*A) ! *** dN1/dy
 B(3,2,1)=0.d0
 B(4,2,1)=(y(2)-y(3))/(2.0d0*A) ! *** dN1/dx
 B(1,3,1)=(y(3)-y(1))/(2.0d0*A) ! *** dN2/dx
 B(2,3,1)=0.d0
 B(3,3,1)=(x(1)-x(3))/(2.0d0*A) ! *** dN2/dy
 B(4,3,1)=0.d0
 B(1,4,1)=0.d0
 B(2,4,1)=(x(1)-x(3))/(2.0d0*A) ! *** dN2/dy
 B(3,4,1)=0.d0
 B(4,4,1)=(y(3)-y(1))/(2.0d0*A) ! *** dN2/dx
 B(1,5,1)=(y(1)-y(2))/(2.0d0*A) ! *** dN3/dx
 B(2,5,1)=0.d0
 B(3,5,1)=(x(2)-x(1))/(2.0d0*A) ! *** dN3/dy
 B(4,5,1)=0.d0
 B(1,6,1)=0.d0
 B(2,6,1)=(x(2)-x(1))/(2.0d0*A) ! *** dN3/dy
 B(3,6,1)=0.d0
 B(4,6,1)=(y(1)-y(2))/(2.0d0*A) ! *** dN3/dx
 W(1)=A

 L(1,1)=(x(1)+x(2)+x(3))/3.d0
 L(2,1)=(y(1)+y(2)+y(3))/3.d0

end subroutine triangleT3
! *******************************************************************
! *** subroutine read_material
! *******************************************************************
subroutine read_material(filename,cname,nstat,nprop,props)
 character(5), intent(out) :: cname
 character(80) :: filename, line
 integer, intent(out) :: nstat,nprop
 real(8), allocatable, intent(out) :: props(:)
 integer :: unit
 open(newunit=unit, file=filename, status='old')
  read(unit,'(a)') line
  read(unit,'(a)') cname
  read(unit,'(a)') line
  read(unit,*) nstat
  read(unit,'(a)') line
  read(unit,*) nprop
  allocate(props(nprop))
  read(unit,'(a)') line
  read(unit,*) props(:)
 close(unit)
end subroutine
! *******************************************************************
! *** subroutine read_msh
! *******************************************************************
subroutine read_msh(filename, coords, conn, nnode, nelem)
 implicit none
 character(len=*), intent(in) :: filename
 integer, intent(out) :: nnode, nelem
 real(8), allocatable, intent(out) :: coords(:,:)
 integer, allocatable, intent(out) :: conn(:,:)
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
  allocate(conn(3, num_elem)) ! sur-allocation temporaire
  nelem = 0
  do i = 1, num_elem
   read(unit,*) id, etype, ntags, tag1, tag2, node1, node2, node3
   if (etype == 2) then ! triangle linéaire
    nelem = nelem + 1
    conn(1, nelem) = node1
    conn(2, nelem) = node2
    conn(3, nelem) = node3
   end if
  end do
  if (nelem < num_elem) then
   call move_alloc(conn, conn)
   allocate(conn(3, nelem))
   conn = conn(:, 1:nelem)
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
 logical :: reading_step,reading_disp,reading_force
 integer :: i,j,nbc,iset,inode,jnode,bc_axis
 character(len=64) :: bc_name
 real(8) :: bc_value

 reading_step=.false.
 reading_disp=.false.
 reading_force=.false.
 bc_type(:)=0
 open(unit=10, file=filename, status='old')
 do
  read(10,'(A)', end=100) line
! Détection des sections
  if (index(line, "*STEP") > 0) then
   reading_step=.true.
   cycle
  endif
  if (reading_step) then
   if (index(line, "*DISPLACEMENT") > 0) then
    reading_disp=.true.
    reading_force=.false.
    cycle
   endif
   if (index(line, "*FORCE") > 0) then
    reading_disp=.false.
    reading_force=.true.
    cycle
   endif
   if (not(reading_disp.or.reading_force)) then
    read(line,*) duration,dtime,dtime_min,dtime_max
    write(*,'("Duration",es12.4)') duration
    write(*,'("Initial time increment",es12.4)') dtime
    write(*,'("Minimal time increment",es12.4)') dtime_min
    write(*,'("Maximal time increment",es12.4)') dtime_max
   endif
  endif
  if (reading_disp) then
   i=index(line,",")
   bc_name=line(1:i-1)
   write(*,'("  Nset for prescribed displacement: ",a)') trim(bc_name)
   j=index(line(i+1:len(line)),",")
   read(line(i+1:i+j),*) bc_axis
   write(*,'("  Axis for prescribed displacement: ",i2)') bc_axis
   read(line(i+j+1:len(line)),*) bc_value
   write(*,'("  Value for prescribed displacement: ",f12.4)') bc_value
   do iset=1,size(nset_names)
    if (bc_name.eq.nset_names(iset)) then
      write(*,*) "found set",iset
     do inode=1,nset_sizes(iset)
      jnode=nset_nodes(iset,inode)
!       write(*,*) jnode
      idof=3*jnode+bc_axis-3
      bc_type(idof)=1
      bc_val(idof)=bc_value
     enddo
    endif
   enddo
   reading_disp=.false.
  endif
  if (reading_force) then
   i=index(line,",")
   bc_name=line(1:i-1)
   write(*,'("  Nset for prescribed force: ",a)') trim(bc_name)
   j=index(line(i+1:len(line)),",")
   read(line(i+1:i+j),*) bc_axis
   write(*,'("  Axis for prescribed force: ",i2)') bc_axis
   read(line(i+j+1:len(line)),*) bc_value
   write(*,'("  Value for prescribed force: ",f12.4)') bc_value
   do iset=1,size(nset_names)
    if (bc_name.eq.nset_names(iset)) then
      write(*,*) "found set",iset
     do inode=1,nset_sizes(iset)
      jnode=nset_nodes(iset,inode)
!       write(*,*) jnode
      idof=3*jnode+bc_axis-3
      bc_type(idof)=2
      bc_val(idof)=bc_value
     enddo
    endif
   enddo
   reading_force=.false.
  endif
 enddo
100 continue
 close(unit=10)

end subroutine
! *******************************************************************
! *** subroutine read_inp
! *******************************************************************
subroutine read_inp(filename,coords,conn,nnode,nelem,nset_names,nset_sizes,nset_nodes)
 implicit none
 character(len=*), intent(in) :: filename
 integer, intent(out) :: nnode, nelem
 real(8), allocatable, intent(out) :: coords(:,:)
 integer, allocatable, intent(out) :: conn(:,:)
 integer :: max_nodes,max_nset_size,max_nsets
 real(8) :: x,y,z
 character(len=256) :: line
 character(len=64), allocatable :: nset_names(:)
 integer, allocatable :: nset_sizes(:),nset_nodes(:,:)
 integer :: i, j, k, id, n(10)
 integer :: current_nset, node_id, ierr
 logical :: reading_nodes, reading_elems, reading_nset
 character(len=64) :: temp_name

 nnode=0
 nelem=0
 max_nsets=0
 max_nset_size=0
 reading_nodes=.false.
 reading_elems=.false.
 reading_nset=.false.
 open(unit=10, file=filename, status='old')
 do
  read(10,'(A)', end=100) line
! Détection des sections
  if (index(line, "*NODE") > 0) then
   reading_nodes = .true.
   reading_elems = .false.
   reading_nset  = .false.
   cycle
  else if (index(line, "*ELEMENT") > 0 .and. index(line, "C3D8") > 0) then
   reading_nodes = .false.
   reading_elems = .true.
   reading_nset  = .false.
   nnode_per_element=8
   cycle
  else if (index(line, "*ELEMENT") > 0 .and. index(line, "C3D4") > 0) then
   reading_nodes = .false.
   reading_elems = .true.
   reading_nset  = .false.
   nnode_per_element=4
   cycle
  else if (index(line, "*ELEMENT") > 0 .and. index(line, "C3D10") > 0) then
   reading_nodes = .false.
   reading_elems = .true.
   reading_nset  = .false.
   nnode_per_element=10
   cycle
  else if (index(line, "*NSET") > 0) then
   reading_nodes = .false.
   reading_elems = .false.
   reading_nset  = .true.
! Extraire le nom du NSET
   call parse_nset_name(line, temp_name)
   max_nsets = max_nsets + 1
   max_nodes = 0
   cycle
  else if (index(line, "*") == 1) then
   reading_nodes = .false.
   reading_elems = .false.
   reading_nset  = .false.
   cycle
  end if
! Lecture des données
  if (reading_nodes) then
   write(*,'("READING NODE")')
   read(line,*,iostat=ierr) id,x,y,z
   if (ierr.eq.0) then
    nnode=nnode+1
   endif
   write(*,'("READING NODE",i4)') id
  else if (reading_elems) then
   write(*,'("READING ELEMS")')
   read(line,*,iostat=ierr) id, (n(inode),inode=1,nnode_per_element)
   if (ierr.eq.0) then
    nelem=nelem+1
   endif
  else if (reading_nset) then
   write(*,'("READING NSET")')
   write(*,'(i4)') count_commas(line)
   write(*,'("READING LINE in NSET")')
   do i = 1, count_commas(line)
    read(line,*) node_id
    max_nodes = max_nodes + 1
   end do
   max_nset_size=max(max_nset_size,max_nodes)
   write(*,*) max_nsets,max_nodes
  end if
 end do
100 continue
 close(10)
 print *, "Number of nodes:", nnode
 print *, "Number of elements:", nelem
 print *, "Number of nodes per element:", nnode_per_element
 print *, "Number of NSETs:", max_nsets
 print *, "Max number of nodes / NSET:", max_nset_size
 allocate(coords(3, nnode))
 allocate(conn(nnode_per_element, nelem))
 allocate(nset_sizes(max_nsets),nset_nodes(max_nsets, max_nset_size),nset_names(max_nsets))
 nnode=0
 nelem=0
 current_nset=0
 nset_sizes=0
 reading_nodes=.false.
 reading_elems=.false.
 reading_nset=.false.
 open(unit=10, file=filename, status='old')
 do
  read(10,'(A)', end=200) line
! Détection des sections
  if (index(line, "*NODE") > 0) then
   write(*,'("NODES")')
   reading_nodes = .true.
   reading_elems = .false.
   reading_nset  = .false.
   cycle
  else if (index(line, "*ELEMENT") > 0 .and. index(line, "C3D8") > 0) then
   reading_nodes = .false.
   reading_elems = .true.
   reading_nset  = .false.
   nnode_per_element=8
   nntg_per_element=8
   cycle
  else if (index(line, "*ELEMENT") > 0 .and. index(line, "C3D4") > 0) then
   reading_nodes = .false.
   reading_elems = .true.
   reading_nset  = .false.
   nnode_per_element=4
   nntg_per_element=1
   cycle
  else if (index(line, "*ELEMENT") > 0 .and. index(line, "C3D10") > 0) then
   reading_nodes = .false.
   reading_elems = .true.
   reading_nset  = .false.
   nnode_per_element=10
   nntg_per_element=4
   cycle
  else if (index(line, "*NSET") > 0) then
   write(*,'("NSET")')
   reading_nodes = .false.
   reading_elems = .false.
   reading_nset  = .true.
! Extraire le nom du NSET
   call parse_nset_name(line, temp_name)
   current_nset = current_nset + 1
   nset_names(current_nset) = temp_name
   cycle
  else if (index(line, "*") == 1) then
   reading_nodes = .false.
   reading_elems = .false.
   reading_nset  = .false.
   cycle
  end if
! Lecture des données
  if (reading_nodes) then
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
  else if (reading_elems) then
   write(*,'("READING ELEMS")')
   read(line,*,iostat=ierr) id, (n(inode),inode=1,nnode_per_element)
   if (ierr.eq.0) then
    nelem=nelem+1
    do inode=1,nnode_per_element
     conn(inode,nelem)=n(inode)
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
    j=index(line, ",")
    read(line(k:j-1),*) node_id
    line(j:j)=" "
    k=j+1
    nset_sizes(current_nset)=nset_sizes(current_nset)+1
    nset_nodes(current_nset,nset_sizes(current_nset))=node_id
   end do
  end if
 end do

200 continue
 close(10)

 do id=1,nnode
  write(*,'("Node",i4," Coords",3es12.4)') id,coords(:,id)
 enddo

 do id=1,nelem
  write(*,'("Element",i4," Nodes",6i4)') id,conn(:,id)
 enddo

 do id=1,current_nset
  write(*,'("Nset",i4," Size",i4)') id,nset_sizes(id)
  do i=1,nset_sizes(id)
   write(*,'("    Node:",i4)') nset_nodes(id,i)
  enddo
 enddo

end subroutine
! *******************************************************************
! *** subroutine parse_nset_name
! *******************************************************************
subroutine parse_nset_name(line, name)
 character(len=*), intent(in) :: line
 character(len=*), intent(out) :: name
 integer :: p1
 name = ""
 p1 = index(line, "NSET=")
 if (p1 > 0) then
  read(line(p1+5:), '(A)') name
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
                  ,internalenergy_old(intg),internalenergy_new(intg),dissipatedenergy_old(intg),dissipatedenergy_new(intg),modulus_old(intg))
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
 type(Tensor2s) :: dstress

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
 dissipatedenergy_new=dissipatedenergy_old+0.5d0*((stress_old+stress_new).ddot.(plastic_strain_new-plastic_strain_old))-(0.5d0*(yield_old+yield_new)-threshold)*(cumulated_new-cumulated_old)
! *** calculate final internal energy density
 internalenergy_new=0.5d0*(stress_new.ddot.elastic_strain_new)+hardening/(coefficient+1.d0)*cumulated_new**(coefficient+1.d0)
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
