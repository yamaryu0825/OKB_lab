
!**************************パラメータの定義モジュール*****************************
module para
		implicit none
	!解析毎に変更するGA格納用パラメータ
		integer,parameter :: ng = 5															!整数型定数ng=5を定義
		integer,parameter :: nn = 100 													!1世代当たりの個体数nnをここで変更
	!異方性繊維の弾性定数 繊維主軸は1 等方面は2-3面
		double precision :: E1_f(nn) 														!= 230.d0 [GPa]
		double precision :: E2_f(nn)  													!= 17.5d0 [GPa]
		double precision :: G12_f(nn)														!= 37.d0  [GPa]
  	double precision :: G23_f(nn)														!等方面内
		double precision :: v12_f(nn)														!= 0.46d0
		double precision :: v23_f(nn) 													!= 0.17d0
	!等方性樹脂の弾性定数
		double precision, parameter :: E_m =  1.8d0 						![GPa]樹脂のヤング率
		double precision, parameter :: v_m = 	0.33d0						!	[-] 樹脂のポアソン比
  	double precision, parameter :: G_m  = E_m/(2*(1+v_m))   ![GPa]横弾性係数
	!炭素繊維の体積含有率
	!double precision, parameter :: vf = 0.17d0 !(TANK材)
		double precision :: Vf(nn)															!Vfを変数に変更
	!TANC材を計算する際のパラメータ
		integer :: N = 720																			!	層数
		double precision, parameter :: pi = 3.14159265359				!	π
end module para
!**************************パラメータの定義モジュールの終了***********************


!**************************計算用変数の定義モジュール*****************************
module variables
		use para																		!パラメータ定義
		double precision :: E1,E2,G12,G23,v12,v23		!出力用配列
		double precision :: EGv_yomu(nn,ng)   			!繊維弾性定数読込用配列
		double precision :: ec(nn,9)				 				!等価剛性(出力結果)格納用配列
		double precision, dimension(6,6) :: D, C		!use eshelby
		double precision, dimension(6,6) :: u, pred	!use inverse_D
		integer :: i, j, k, l, m, ii								!カウンター変数の宣言
end module variables
!**************************計算用変数の定義モジュールの終了***********************


!**************************EMT_layerwiseのメインプログラム***********************
PROGRAM layerwise
		use para							!パラメータ定義
		use variables					!変数定義
		call input						!データ読込
	!----------------------1~nnまでの各個体の等価剛性を計算------------------------
		do ii = 1,nn
			call eshelby				!UD材CFRP積層板の弾性剛性マトリクスの作成
			call layer					!TANC材の弾性剛性マトリクスの作成
		end do
	!-----------------------等価剛性計算の終了-------------------------------------
		call output						!データ出力
		stop
END PROGRAM
!**************************EMT_layerwiseのメインプログラムの終了*****************


!**************************データ読込用サブルーチン******************************
SUBROUTINE input
		use para
		use variables
	!make_ini.f90でランダム生成された繊維弾性定数
		open(10,file ='elastic_modulus.txt',status = 'old')
	!i行目の数値を読み込み、EGv_yomu、Vf(i)に格納
		do i = 1,nn
	    read(10,'(F15.5,F15.5,F15.6,F15.6,F15.5,F15.6)') EGv_yomu(i,:), Vf(i)
		end do
 	!読み込んだ弾性定数を計算用変数E1～G23に格納
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
!**************************データ読込用サブルーチンの終了*************************


!**************************エシェルビーの関数(一方向材等価剛性を求める)***********
SUBROUTINE eshelby
		use para
		use variables
	!変数定義
		double precision, dimension(6,6) :: Dm, Df, S, Dk	!各弾性剛性マトリクスおよびテンソル
		double precision, dimension(6,6) :: D1, D2, D3		!各計算項
		double precision :: delta
	!炭素繊維のマトリクスファイルを開く？
		open(10,file='matrix.csv',status='replace')
	!炭素繊維の剛性マトリクス
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
	!座標変換
		Df(1,1) = Dk(3,3); Df(2,1) = Dk(2,3)
		Df(1,2) = Dk(3,2); Df(2,2) = Dk(2,2)
		Df(1,3) = Dk(3,1); Df(2,3) = Dk(2,1)
		Df(3,1) = Dk(1,3); Df(4,4) = Dk(6,6)
		Df(3,2) = Dk(1,2); Df(5,5) = Dk(5,5)
		Df(3,3) = Dk(1,1); Df(6,6) = Dk(4,4)
	!更新
		Dk = 0.0
	!樹脂の剛性マトリクス
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
	!エシェルビーテンソル
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
	!------------------一方向材等価剛性マトリクスの作成-----------------------------
	!第1項
		D1 = Dm
	!第2項
		pred = (1-vf(ii))*(Df-Dm)
		pred = MATMUL(pred,S)
		pred = pred + Dm
		call inverse_D
		D2 = u
	!第3項
		D3 = Df - Dm
		D3 = MATMUL(D3,S)
		D3 = D3 + Dm
		D3 = (1-vf(ii))*D3
		D3 = D3 + (vf(ii)*Df)
		D = MATMUL(D1,D2)
		D = MATMUL(D,D3)
		Dk = D
	!座標変換
		D(1,1) = Dk(3,3); D(2,1) = Dk(2,3)
		D(1,2) = Dk(3,2); D(2,2) = Dk(2,2)
		D(1,3) = Dk(3,1); D(2,3) = Dk(2,1)
		D(3,1) = Dk(1,3); D(4,4) = Dk(6,6)
		D(3,2) = Dk(1,2); D(5,5) = Dk(5,5)
		D(3,3) = Dk(1,1); D(6,6) = Dk(4,4)
	!コンプライアンスマトリクスの作成-
		pred = D
		call inverse_D
		C = u
	!弾性定数への変換
		ec(ii,1) =   1.0 / C(1,1)		!	E11	[GPa]
		ec(ii,2) =   1.0 / C(2,2)		!	E22 [GPa]
		ec(ii,3) =   1.0 / C(3,3)		!	E33 [GPa]
		ec(ii,4) = - ec(ii,1) * C(1,2)	!	V12
		ec(ii,5) = - ec(ii,1) * C(1,3)	!	V13
		ec(ii,6) = - ec(ii,2) * C(2,3)	!	V23
		ec(ii,7) =   1.0 / C(6,6)		!	G12 [GPa]
		ec(ii,8) =   1.0 / C(5,5)		!	G31 [GPa]
		ec(ii,9) =   1.0 / C(4,4)		!	G23 [GPa]
	!------------------一方向材等価剛性マトリクスの作成の終了-----------------------
	!------------------エシェルビーの計算結果の出力---------------------------------
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
	!------------------エシェルビーの計算結果の出力の終了---------------------------
		return
END SUBROUTINE
!**************************エシェルビーの関数(一方向材等価剛性を求める)の終了******


!-----レイヤーワイズの関数(TANK材の弾性定数を求める)------------------------------------------
SUBROUTINE layer
use para
use variables
!	回転作用素
	double precision, dimension(6,6) :: T_s, T_e
!	角度のついた一方向CFRPの弾性剛性マトリクス
	double precision, dimension(6,6) :: Dcr
!	TANC材の弾性剛性マトリクス
	double precision, dimension(6,6) :: Dc
!	TANC材のコンプライアンスマトリクス
	double precision, dimension(6,6) :: Cc

	open(100,file='matrix_tanc.csv',status='replace')

!	各マトリクスの初期化
	T_s = 0.0
	Dcr = 0.0
	Dc  = 0.0
!	角度thetaの初期設定
	theta = 0.0
!	層ごとの回転行列の作成およびマトリクスの計算
	do i = 1, N+1
!		回転作用素の作成
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
!		回転作用素による計算
		Dcr = MATMUL(T_s,D)
		Dcr = MATMUL(Dcr,T_e)
!		弾性剛性マトリクスの合成
		Dc = Dc + Dcr
!		角度thetaの更新
		theta = theta + (pi/N)
!		print *, " N FINISH "
	enddo
!	弾性剛性マトリクスの平均化
	Dc = Dc / (N+1)
!	コンプライアンスマトリクスの作成
	pred = Dc
	call inverse_D
	Cc = u
!弾性定数への変換
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


!------書き出し----------------------------------------------
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
!-----レイヤーワイズの関数の終了----------------------------------------------





!----TANK材の弾性定数を出力----------------------------------------
SUBROUTINE output
use para
use variables

!【次回】ec(ii,:)をec(ii,:)に変更

!出力用ファイルを呼び出す
open(20,file ='composite_GPa.txt',status = 'replace')   !計算後の弾性定数を格納(見る用、単位GPa)
open(30,file ='composite_input.txt',status = 'replace') !計算後の弾性定数を格納(Abaqus用、単位 MPa) 存在しない場合は新規作成


!出力時の一行目の項目名
	write(20,'(7X,A,7X,A,8X,A,6X,A,11X,A,11X,A,11X,A,8X,A,7X,A,7X,A)')&
	& 'E11(GPa)','E22(GPa)','E33(GPa)','ν12','ν13','ν23','G12(GPa)','G13(GPa)','G23(GPa)','Vf'

do i =1,nn
	write(20,'(F15.5,F15.5,F15.5,F15.6,F15.6,F15.6,F15.5,F15.5,F15.5,F15.5)') ec(i,:),Vf(i)
end do

	!GPaをMPa表記に変更
do i =1,nn
	  ec(i,1) = ec(i,1)*1.d+03
	  ec(i,2) = ec(i,2)*1.d+03
		ec(i,3) = ec(i,3)*1.d+03
	  ec(i,7) = ec(i,7)*1.d+03
	  ec(i,8) = ec(i,8)*1.d+03
		ec(i,9) = ec(i,9)*1.d+03
end do
	100 format(F8.1,A,X,F8.1,A,X,F8.1,A,X,F8.6,A,X,F8.6,A,X,F8.6,A,X,F8.1,A,X,F8.1,A,X,F8.1)

	!Abaqus用にMPa表記にした弾性定数を装置番号30に入力
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












!----逆行列を求める関数----------------------------------------------
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
