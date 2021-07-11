%Master Thesis R06323052 ����� 
%calibration 

%home production
t_ness = [513 522 527];       % ���n�ɶ�-���\�ɶ�:[�ίv(�t�Ⱥ�)+�d�~�B�N�D1�B�۸˩ΤƧ�+���\(�t�Y�d�])]-[���\(�t�Y�d�])]  from Data_Time(new)
t_free_work = 24-t_ness/60 ;  % �ۥѮɶ�+�u�@�ɶ�
k0 = zeros(1,3) ; 
k1 = zeros(1,3) ; 
k2 = zeros(1,3) ; 
x0 = [194 180 167];
x1 = [323 313 295];
x2 = [382 366 339];
x0 = x0/60./t_free_work;
x1 = x1/60./t_free_work;
x2 = x2/60./t_free_work;
for i = 1:3
  k0(i) = x0(i);
  k1(i) = x1(i)-k0(i);
  f=@(kk)(k0(i)+k1(i)*2^kk-x2(i));
  [kk, fval_k2]=fsolve(f,0.4844);
  k2(i) = kk;  
end


%price parameter
pb = [0.0982 0.0995 0.1256]; %from calibration ��������ܼ�
pe = [0.0760 0.0611 0.0509]; %from calibration ��������ܼ�
am = ones(1,4)             ; %from calibration ��������ܼ� 
c_ = [0.5723 0.5967 0.5740]; %from calibration ��������ܼ�
%production side parameter
n_L = [5008525 5017600 4462241]; % from �U�@�N�����ͨ|�`��
n_H = [ 534096  767832 1093217]; % from �U�@�N�����ͨ|�`��
N_L = n_L./(n_L+n_H) ; %�зǤ�N_L
N_H = n_H./(n_L+n_H) ; %�зǤ�N_H
%�����ܼ�
b_L      = [2.333 2.189 1.874];  %from �U�@�N�����ͨ|�`��(64)
b_H      = [1.691 1.599 1.636];  %from �U�@�N�����ͨ|�`��(64)
b_total  = [2.309 2.154 1.923];  %from �U�@�N�����ͨ|�`��(64)
wf_L     = [0.4499 0.4749 0.4870];  %from calibration ��������ܼ�
wf_H     = [0.9290 0.9165 0.8721];  %from calibration ��������ܼ�
R        = [0.435 0.432 0.431];  %from PSFD (L->L percentage)
tf_L     = [135 156 172]; %from Human_newnew (min./day)
tf_H     = [226 241 249]; %from Human_newnew (min./day)
tf_L     = tf_L/60./t_free_work ; %�зǤ�(���⦨�p�ɫ�A�b���⦨�ۥѮɶ��e��)
tf_H     = tf_H/60./t_free_work ; %�зǤ�(���⦨�p�ɫ�A�b���⦨�ۥѮɶ��e��)
Z_L      = N_L.*tf_L ; %L_type labor force
Z_H      = N_H.*tf_H ; %H_type labor force
% Assume
rho = 0.401 ; % from Krusell(2000)

%calibration for A pi
A  = zeros(1,3);
pi = zeros(1,3);
for i = 1:3
x0_t = [1.1135 0.5762];
[x,fval_A_pi]=fsolve(@(x)myfun_A_pi(x,Z_L(i),Z_H(i),wf_L(i),wf_H(i),rho),x0_t);
A(i)  = x(1) ; 
pi(i) = x(2) ;
end

%calibration for gamma phi eta betta
phi = zeros(1,3);
eta = zeros(1,3);
for i = 1:3
x0_t = [b_L(i) b_L(i) b_H(i) 0.3448 0.1 ];
N = 1000000 ;
options = optimset('MaxFunEvals',N,'MaxIter',N);
[x,fval]=fsolve(@(x)myfun_5(x,pb(i),pe(i),am(i),c_(i),k0(i),k1(i),k2(i),wf_L(i),wf_H(i),R(i),b_L(i),b_H(i)),x0_t,options);
    phi(i) = x(4);
    eta(i) = x(5);
end
%%

%caliabration_check
bd_L_m      = zeros(1,3);
be_L_m      = zeros(1,3);
b_L_m       = zeros(1,3);
b_H_m       = zeros(1,3);
b_total_m   = zeros(1,3);

j = 1 ; %preference parameter ����ܼ�
for i = 1:3
x0 = [b_L(i), b_L(i), b_H(i)] ;
N = 1000000 ;
options = optimset('MaxFunEvals',N,'MaxIter',N);
[x,fval_m]=fsolve(@(x)myfun_5_c(x,pb(i),pe(i),am(i),c_(i),k0(i),k1(i),k2(i),wf_L(i),wf_H(i),phi(j),eta(j)),x0,options);
    bd_L_m(i)    = x(1) ;
    be_L_m(i)    = x(2) ;
    b_L_m(i)     = R(i)*bd_L_m(i)+(1-R(i))*be_L_m(i) ;
    b_H_m(i)     = x(3) ;
    b_total_m(i) = N_L(i)*(b_L_m(i))+N_H(i)*b_H_m(i) ;
 end
disp('     phi(1)      phi(2)     phi(3)')
disp(phi)
disp('     eta(1)      eta(2)     eta(3)')
disp(eta)