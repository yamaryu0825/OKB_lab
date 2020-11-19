
!**************************�p�����[�^�̒�`���W���[��*****************************
module para
		implicit none
	!��͖��ɕύX����GA�i�[�p�p�����[�^
		integer,parameter :: ng = 5															!�����^�萔ng=5���`
		integer,parameter :: nn = 100 													!1���㓖����̌̐�nn�������ŕύX
	!�ٕ����@�ۂ̒e���萔 �@�ێ厲��1 �����ʂ�2-3��
		double precision :: E1_f(nn) 														!= 230.d0 [GPa]
		double precision :: E2_f(nn)  													!= 17.5d0 [GPa]
		double precision :: G12_f(nn)														!= 37.d0  [GPa]
  	double precision :: G23_f(nn)														!�����ʓ�
		double precision :: v12_f(nn)														!= 0.46d0
		double precision :: v23_f(nn) 													!= 0.17d0
	!�����������̒e���萔
		double precision, parameter :: E_m =  1.8d0 						![GPa]�����̃����O��
		double precision, parameter :: v_m = 	0.33d0						!	[-] �����̃|�A�\����
  	double precision, parameter :: G_m  = E_m/(2*(1+v_m))   ![GPa]���e���W��
	!�Y�f�@�ۂ̑̐ϊܗL��
	!double precision, parameter :: vf = 0.17d0 !(TANK��)
		double precision :: Vf(nn)															!Vf��ϐ��ɕύX
	!TANC�ނ��v�Z����ۂ̃p�����[�^
		integer :: N = 720																			!	�w��
		double precision, parameter :: pi = 3.14159265359				!	��
end module para
!**************************�p�����[�^�̒�`���W���[���̏I��***********************


!**************************�v�Z�p�ϐ��̒�`���W���[��*****************************
module variables
		use para																		!�p�����[�^��`
		double precision :: E1,E2,G12,G23,v12,v23		!�o�͗p�z��
		double precision :: EGv_yomu(nn,ng)   			!�@�ےe���萔�Ǎ��p�z��
		double precision :: ec(nn,9)				 				!��������(�o�͌���)�i�[�p�z��
		double precision, dimension(6,6) :: D, C		!use eshelby
		double precision, dimension(6,6) :: u, pred	!use inverse_D
		integer :: i, j, k, l, m, ii								!�J�E���^�[�ϐ��̐錾
end module variables
!**************************�v�Z�p�ϐ��̒�`���W���[���̏I��***********************


!**************************EMT_layerwise�̃��C���v���O����***********************
PROGRAM layerwise
		use para							!�p�����[�^��`
		use variables					!�ϐ���`
		call input						!�f�[�^�Ǎ�
	!----------------------1~nn�܂ł̊e�̂̓����������v�Z------------------------
		do ii = 1,nn
			call eshelby				!UD��CFRP�ϑw�̒e�������}�g���N�X�̍쐬
			call layer					!TANC�ނ̒e�������}�g���N�X�̍쐬
		end do
	!-----------------------���������v�Z�̏I��-------------------------------------
		call output						!�f�[�^�o��
		stop
END PROGRAM
!**************************EMT_layerwise�̃��C���v���O�����̏I��*****************


!**************************�f�[�^�Ǎ��p�T�u���[�`��******************************
SUBROUTINE input
		use para
		use variables
	!make_ini.f90�Ń����_���������ꂽ�@�ےe���萔
		open(10,file ='elastic_modulus.txt',status = 'old')
	!i�s�ڂ̐��l��ǂݍ��݁AEGv_yomu�AVf(i)�Ɋi�[
		do i = 1,nn
	    read(10,'(F15.5,F15.5,F15.6,F15.6,F15.5,F15.6)') EGv_yomu(i,:), Vf(i)
		end do
 	!�ǂݍ��񂾒e���萔���v�Z�p�ϐ�E1�`G23�Ɋi�[
		do i = 1,nn
  		E1_f(i) = EGv_yomu(i,1)
  		E2_f(i) = EGv_yomu(i,2)
  		v12_f(i) = EGv_yomu(i,3)
  		v23_f(i) = EGv_yomu(i,4)
  		G12_f(i) = EGv_yomu(i,5)
  		G23_f(i) = E2_f(i)/(2.d0*(v23_f(i)+1.d0))
		end do
		close(10)
END SUBROUTINE
!**************************�f�[�^�Ǎ��p�T�u���[�`���̏I��*************************


!**************************�G�V�F���r�[�̊֐�(������ޓ������������߂�)***********
SUBROUTINE eshelby
		use para
		use variables
	!�ϐ���`
		double precision, dimension(6,6) :: Dm, Df, S, Dk	!�e�e�������}�g���N�X����уe���\��
		double precision, dimension(6,6) :: D1, D2, D3		!�e�v�Z��
		double precision :: delta
	!�Y�f�@�ۂ̃}�g���N�X�t�@�C�����J���H
		open(10,file='matrix.csv',status='replace')
	!�Y�f�@�ۂ̍����}�g���N�X
		vf21 = v12_f(ii) * E2_f(ii) / E1_f(ii)
		Gf23 = E2_f(ii)/(2*(1+v23_f(ii)))
		delta = (1.0-(2.0*v12_f(ii)*vf21)-(v23_f(ii)*v23_f(ii))-(2.0*v12_f(ii)*v23_f(ii)*vf21)) / (E1_f(ii)*E2_f(ii)*E2_f(ii))
		Df = 0.0
		Df(1,1) = (1.0-(v23_f(ii)*v23_f(ii))) / (E2_f(ii)*E2_f(ii)*delta)
		Df(1,2) = (v12_f(ii)+(v12_f(ii)*v23_f(ii))) / (E1_f(ii)*E2_f(ii)*delta)
		Df(1,3) = Df(1,2)
		Df(2,1) = Df(1,2)
		Df(2,2) = (1.0-(v12_f(ii)*vf21)) / (E1_f(ii)*E2_f(ii)*delta)
		Df(2,3) = (v23_f(ii)+(v12_f(ii)*vf21)) / (E1_f(ii)*E2_f(ii)*delta)
		Df(3,1) = Df(1,2)
		Df(3,2) = Df(2,3)
		Df(3,3) = Df(2,2)
		Df(4,4) = E2_f(ii) / (2.0*(1.0+v23_f(ii)))
		Df(5,5) = G12_f(ii)
		Df(6,6) = G12_f(ii)
		Dk = Df
	!���W�ϊ�
		Df(1,1) = Dk(3,3); Df(2,1) = Dk(2,3)
		Df(1,2) = Dk(3,2); Df(2,2) = Dk(2,2)
		Df(1,3) = Dk(3,1); Df(2,3) = Dk(2,1)
		Df(3,1) = Dk(1,3); Df(4,4) = Dk(6,6)
		Df(3,2) = Dk(1,2); Df(5,5) = Dk(5,5)
		Df(3,3) = Dk(1,1); Df(6,6) = Dk(4,4)
	!�X�V
		Dk = 0.0
	!�����̍����}�g���N�X
		Dm = 0.0
		Dm(1,1) = (E_m*(1.0-v_m)) / ((1+v_m)*(1.0-(2.0*v_m)))
		Dm(1,2) = (E_m*v_m) / ((1+v_m)*(1.0-(2.0*v_m)))
		Dm(1,3) = Dm(1,2)
		Dm(2,1) = Dm(1,2)
		Dm(2,2) = Dm(1,1)
		Dm(2,3) = Dm(1,2)
		Dm(3,1) = Dm(1,2)
		Dm(3,2) = Dm(1,2)
		Dm(3,3) = Dm(1,1)
		Dm(4,4) = E_m / (2.0*(1.0+v_m))
		Dm(5,5) = Dm(4,4)
		Dm(6,6) = Dm(4,4)
	!�G�V�F���r�[�e���\��
		S = 0.0
		S(1,1) = (5.0-(4.0*v_m)) / (8.0*(1.0-v_m))
		S(1,2) = ((4.0*v_m)-1.0) / (8.0*(1.0-v_m))
		S(1,3) = v_m / (2.0*(1.0-v_m))
		S(2,1) = S(1,2)
		S(2,2) = S(1,1)
		S(2,3) = S(1,3)
		S(3,1) = 0.0
		S(3,2) = 0.0
		S(3,3) = 0.0
		S(4,4) = 0.5
		S(5,5) = 0.5
		S(6,6) = (3.0-(4.0*v_m)) / (4.0*(1.0-v_m))
	!------------------������ޓ��������}�g���N�X�̍쐬-----------------------------
	!��1��
		D1 = Dm
	!��2��
		pred = (1-vf(ii))*(Df-Dm)
		pred = MATMUL(pred,S)
		pred = pred + Dm
		call inverse_D
		D2 = u
	!��3��
		D3 = Df - Dm
		D3 = MATMUL(D3,S)
		D3 = D3 + Dm
		D3 = (1-vf(ii))*D3
		D3 = D3 + (vf(ii)*Df)
		D = MATMUL(D1,D2)
		D = MATMUL(D,D3)
		Dk = D
	!���W�ϊ�
		D(1,1) = Dk(3,3); D(2,1) = Dk(2,3)
		D(1,2) = Dk(3,2); D(2,2) = Dk(2,2)
		D(1,3) = Dk(3,1); D(2,3) = Dk(2,1)
		D(3,1) = Dk(1,3); D(4,4) = Dk(6,6)
		D(3,2) = Dk(1,2); D(5,5) = Dk(5,5)
		D(3,3) = Dk(1,1); D(6,6) = Dk(4,4)
	!�R���v���C�A���X�}�g���N�X�̍쐬-
		pred = D
		call inverse_D
		C = u
	!�e���萔�ւ̕ϊ�
		ec(ii,1) =   1.0 / C(1,1)		!	E11	[GPa]
		ec(ii,2) =   1.0 / C(2,2)		!	E22 [GPa]
		ec(ii,3) =   1.0 / C(3,3)		!	E33 [GPa]
		ec(ii,4) = - ec(ii,1) * C(1,2)	!	V12
		ec(ii,5) = - ec(ii,1) * C(1,3)	!	V13
		ec(ii,6) = - ec(ii,2) * C(2,3)	!	V23
		ec(ii,7) =   1.0 / C(6,6)		!	G12 [GPa]
		ec(ii,8) =   1.0 / C(5,5)		!	G31 [GPa]
		ec(ii,9) =   1.0 / C(4,4)		!	G23 [GPa]
	!------------------������ޓ��������}�g���N�X�̍쐬�̏I��-----------------------
	!------------------�G�V�F���r�[�̌v�Z���ʂ̏o��---------------------------------
		write(10,*) " The elastic stiffness matrix of unidirectional CFRP "
		do i = 1, 6
			write(10,10) D(i,:)
		enddo
		write(10,*) " The compliance matrix of unidirectional CFRP "
		do i = 1, 6
			write(10,10) C(i,:)
		enddo
		10 format(6f14.6)
		close(10)
	!------------------�G�V�F���r�[�̌v�Z���ʂ̏o�͂̏I��---------------------------
		return
END SUBROUTINE
!**************************�G�V�F���r�[�̊֐�(������ޓ������������߂�)�̏I��******


!-----���C���[���C�Y�̊֐�(TANK�ނ̒e���萔�����߂�)------------------------------------------
SUBROUTINE layer
use para
use variables
!	��]��p�f
	double precision, dimension(6,6) :: T_s, T_e
!	�p�x�̂��������CFRP�̒e�������}�g���N�X
	double precision, dimension(6,6) :: Dcr
!	TANC�ނ̒e�������}�g���N�X
	double precision, dimension(6,6) :: Dc
!	TANC�ނ̃R���v���C�A���X�}�g���N�X
	double precision, dimension(6,6) :: Cc

	open(100,file='matrix_tanc.csv',status='replace')

!	�e�}�g���N�X�̏�����
	T_s = 0.0
	Dcr = 0.0
	Dc  = 0.0
!	�p�xtheta�̏����ݒ�
	theta = 0.0
!	�w���Ƃ̉�]�s��̍쐬����у}�g���N�X�̌v�Z
	do i = 1, N+1
!		��]��p�f�̍쐬
		T_s(1,1) =   cos(theta)**2.0
		T_s(1,2) =   sin(theta)**2.0
		T_s(1,6) =   sin(2.0*theta)
		T_s(2,1) =   sin(theta)**2.0
		T_s(2,2) =   cos(theta)**2.0
		T_s(2,6) = - sin(2.0*theta)
		T_s(3,3) =   1.0
		T_s(4,4) =   cos(theta)
		T_s(4,5) = - sin(theta)
		T_s(5,4) =   sin(theta)
		T_s(5,5) =   cos(theta)
		T_s(6,1) = - sin(2.0*theta) / 2.0
		T_s(6,2) =   sin(2.0*theta) / 2.0
		T_s(6,6) =   cos(2.0*theta)
		T_e = TRANSPOSE(T_s)
!		��]��p�f�ɂ��v�Z
		Dcr = MATMUL(T_s,D)
		Dcr = MATMUL(Dcr,T_e)
!		�e�������}�g���N�X�̍���
		Dc = Dc + Dcr
!		�p�xtheta�̍X�V
		theta = theta + (pi/N)
!		print *, " N FINISH "
	enddo
!	�e�������}�g���N�X�̕��ω�
	Dc = Dc / (N+1)
!	�R���v���C�A���X�}�g���N�X�̍쐬
	pred = Dc
	call inverse_D
	Cc = u
!�e���萔�ւ̕ϊ�
	ec(ii,1) =   1.0 / Cc(1,1)		!	E11	[GPa]
	ec(ii,2) =   1.0 / Cc(2,2)		!	E22 [GPa]
	ec(ii,3) =   1.0 / Cc(3,3)		!	E33 [GPa]
	ec(ii,4) = - ec(ii,1) * Cc(1,2)	!	V12
	ec(ii,5) = - ec(ii,1) * Cc(1,3)	!	V13
	ec(ii,6) = - ec(ii,2) * Cc(2,3)	!	V23
	ec(ii,7) =   1.0 / Cc(6,6)		!	G12 [GPa]
	ec(ii,8) =   1.0 / Cc(5,5)		!	G31 [GPa]
	ec(ii,9) =   1.0 / Cc(4,4)		!	G23 [GPa]
!-----------------------------------------------------------


!------�����o��----------------------------------------------
	write(100,*) " The elastic stiffness matrix of TANC "
	do i = 1, 6
		write(100,10) Dc(i,:)
	enddo
	write(100,*) " The compliance matrix of TANC "
	do i = 1, 6
		write(100,10) Cc(i,:)
	enddo
	10 format(6f14.6)

	close(100)
	close(1000)
!-----------------------------------------------------------


return
END SUBROUTINE
!-----���C���[���C�Y�̊֐��̏I��----------------------------------------------





!----TANK�ނ̒e���萔���o��----------------------------------------
SUBROUTINE output
use para
use variables

!�y����zec(ii,:)��ec(ii,:)�ɕύX

!�o�͗p�t�@�C�����Ăяo��
open(20,file ='composite_GPa.txt',status = 'replace')   !�v�Z��̒e���萔���i�[(����p�A�P��GPa)
open(30,file ='composite_input.txt',status = 'replace') !�v�Z��̒e���萔���i�[(Abaqus�p�A�P�� MPa) ���݂��Ȃ��ꍇ�͐V�K�쐬


!�o�͎��̈�s�ڂ̍��ږ�
	write(20,'(7X,A,7X,A,8X,A,6X,A,11X,A,11X,A,11X,A,8X,A,7X,A,7X,A)')&
	& 'E11(GPa)','E22(GPa)','E33(GPa)','��12','��13','��23','G12(GPa)','G13(GPa)','G23(GPa)','Vf'

do i =1,nn
	write(20,'(F15.5,F15.5,F15.5,F15.6,F15.6,F15.6,F15.5,F15.5,F15.5,F15.5)') ec(i,:),Vf(i)
end do

	!GPa��MPa�\�L�ɕύX
do i =1,nn
	  ec(i,1) = ec(i,1)*1.d+03
	  ec(i,2) = ec(i,2)*1.d+03
		ec(i,3) = ec(i,3)*1.d+03
	  ec(i,7) = ec(i,7)*1.d+03
	  ec(i,8) = ec(i,8)*1.d+03
		ec(i,9) = ec(i,9)*1.d+03
end do
	100 format(F8.1,A,X,F8.1,A,X,F8.1,A,X,F8.6,A,X,F8.6,A,X,F8.6,A,X,F8.1,A,X,F8.1,A,X,F8.1)

	!Abaqus�p��MPa�\�L�ɂ����e���萔�𑕒u�ԍ�30�ɓ���
do i =1,nn
	write(30,100)  ec(i,1),',',ec(i,2),',',ec(i,3),',',ec(i,4),',',ec(i,5),',',ec(i,6),',',ec(i,7),',',ec(i,8),',',&
	& ec(i,9)
end do

close(20)
close(30)

!   print *, "---Elastic modulus of TANK material---"
!		write(*,*) 'E11=',ec(ii,1)
!		write(*,*) 'E22=',ec(ii,2)
!		write(*,*) 'V12=',ec(ii,4)
!		write(*,*) 'V13=',ec(ii,5)
!		write(*,*) 'V23=',ec(ii,6)
!		write(*,*) 'G12=',ec(ii,7)
!		write(*,*) 'G13=',ec(ii,8)
!		write(*,*) 'G23=',ec(ii,9)
!	print *, "--------------------------------------"
!-----------------------------------------------------------

END SUBROUTINE












!----�t�s������߂�֐�----------------------------------------------
SUBROUTINE inverse_D
use para
use variables
	double precision, dimension(6) :: h1_D, h2_D
	double precision :: max_in
	integer :: imax

	u = 0.0
	do i = 1, 6
		u(i,i) = 1.0
	enddo

	do i = 1, 6
		imax = i
		do j = i+1, 6
			if (abs(pred(j,1)).gt.abs(pred(imax,1))) then
				imax = j
				h1_D(:) = pred(imax,:)
				pred(j,:) = pred(i,:)
				pred(i,:) = h1_D(:)
				h2_D(:) = u(imax,:)
				u(j,:) = u(i,:)
				u(i,:) = h2_D(:)
			endif
		enddo
		max_in = pred(i,i)
		pred(i,:) = pred(i,:) / max_in
		u(i,:) = u(i,:) / max_in
		do j = 1, 6
			if (j.ne.i) then
				max_in = pred(j,i) / pred(i,i)
				pred(j,:) = pred(j,:) - (max_in*pred(i,:))
				u(j,:) = u(j,:) - (max_in*u(i,:))
			endif
		enddo
	enddo

return
END SUBROUTINE
!-----------------------------------------------------------
