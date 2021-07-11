%Master Thesis R06323052 朝霍 
%calibration 

%home production
t_sleep         = [513 522 527]; %何痸(と何) from Data_Time(new)
t_sleep_shower  = [564 574 578]; %何痸(と何)+縟瑍∟疍帝杆 from Data_Time(new)
t_sleep_eat     = [591 600 605]; %何痸(と何)+ノ繺(甦) from Data_Time(new)
t_sleep_shower_eat = [642 652 656]; %ゲ璶丁: 何痸(と何)+縟瑍∟疍帝杆┪てЗ+ノ繺(甦) from Data_Time(new)
t_ness = t_sleep_eat; 

t_free_work = 24-t_ness/60 ;  % パ丁+丁
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
%part 1 original(including 1822)
pb_a_o = [0.0982 0.0995 0.1256]; % a: average(だ糷舦)
pe_a_o = [0.0760 0.0611 0.0509];
c_a_o  = [0.5723 0.5967 0.5740];
pb_t_o = [0.0982 0.0995 0.1256];
pe_t_o = [0.0760 0.0611 0.0509];
c_t_o  = [0.7350 0.6953 0.6390];
%part 2 new(excluding 1822)
pb_a_n = [0.0965 0.0954 0.1251]; % a: average(だ糷舦)
pe_a_n = [0.0933 0.0879 0.0790];
c_a_n  = [0.5820 0.6061 0.5796];
pb_t_n = [0.0965 0.0954 0.1251];
pe_t_n = [0.1038 0.0941 0.0833];
c_t_n  = [0.6978 0.6745 0.6264];

pb = pb_a_o ;    %from calibration 基闽跑计
pe = pe_a_o ;    %from calibration 基闽跑计
c_ = c_a_o ;     %from calibration 基闽跑计
am = ones(1,4) ; %from calibration 基闽跑计 

%production side parameter
n_L = [5008525 5017600 4462241]; % from キАネ▅羆计
n_H = [ 534096  767832 1093217]; % from キАネ▅羆计
N_L = n_L./(n_L+n_H) ; %夹非てN_L
N_H = n_H./(n_L+n_H) ; %夹非てN_H
%ずネ跑计
b_L      = [2.333 2.189 1.874];  %from キАネ▅羆计(64)
b_H      = [1.691 1.599 1.636];  %from キАネ▅羆计(64)
b_total  = [2.309 2.154 1.923];  %from キАネ▅羆计(64)
wf_L     = [0.4499 0.4749 0.4870];  %from calibration 基闽跑计
wf_H     = [0.9290 0.9165 0.8721];  %from calibration 基闽跑计
R        = [0.435 0.432 0.431];  %from PSFD (L->L percentage)
tf_L     = [135 156 172]; %from Human_newnew (min./day)
tf_H     = [226 241 249]; %from Human_newnew (min./day)
tf_L     = tf_L/60./t_free_work ; %夹非て(传衡Θ传衡Θパ丁ゑ)
tf_H     = tf_H/60./t_free_work ; %夹非て(传衡Θ传衡Θパ丁ゑ)
Z_L      = N_L.*tf_L ; %L_type labor force
Z_H      = N_H.*tf_H ; %H_type labor force
% Assume
rho = 0.401 ; % from Krusell(2000)


%calibration for gamma phi eta betta
bd_L_m      = zeros(1,3);
be_L_m      = zeros(1,3);
b_L_m       = zeros(1,3);
b_H_m       = zeros(1,3);
b_total_m   = zeros(1,3);
tfd_L_m     = zeros(1,3);
tfe_L_m     = zeros(1,3);
tf_L_m      = zeros(1,3);
tf_H_m      = zeros(1,3);

alpha  = zeros(1,3);
gamma  = zeros(1,3);
for i = 1:3
x0_t = [b_L(i) b_L(i) b_H(i) R(i) 1.1418 1.3558 0.7686];
N = 10000000 ;
options = optimset('MaxFunEvals',N,'MaxIter',N,'Tolfun',1e-5);
[x,fva_cal]=fsolve(@(x)myfun_4(x,pb(i),pe(i),am(i),c_(i),k0(i),k1(i),k2(i),wf_L(i),wf_H(i),b_L(i),b_H(i),R(i)),x0_t,options);
    bd_L_m(i)    = x(1) ;
    be_L_m(i)    = x(2) ;
    b_L_m(i)     = R(i)*bd_L_m(i)+(1-R(i))*be_L_m(i) ;
    b_H_m(i)     = x(3) ;
    b_total_m(i) = N_L(i)*(b_L_m(i))+N_H(i)*b_H_m(i) ;
    alpha(i)     = x(4) ;
    gamma(i)     = x(5) ;
    tfd_L_m(i)   = (1-k0(i)-((gamma(i)*pb(i)*bd_L_m(i)+(gamma(i)*k2(i)+alpha(i))*wf_L(i)*k1(i)*(bd_L_m(i))^(k2(i)))/(alpha(i)*wf_L(i))));
    tfe_L_m(i)   = (1-k0(i)-((gamma(i)*(pb(i)+pe(i))*be_L_m(i)+(gamma(i)*k2(i)+alpha(i))*wf_L(i)*k1(i)*(be_L_m(i))^(k2(i)))/(alpha(i)*wf_L(i))));
    tf_L_m(i)    = R(i)*tfd_L_m(i)+(1-R(i))*tfe_L_m(i) ;
    tf_H_m(i)    = (1-k0(i)-((gamma(i)*(pb(i)+pe(i))*b_H_m(i)+(gamma(i)*k2(i)+alpha(i))*wf_H(i)*k1(i)*(b_H_m(i))^(k2(i)))/(alpha(i)*wf_H(i))));
end


%%
%caliabration_check
bd_L_m      = zeros(1,3);
be_L_m      = zeros(1,3);
b_L_m       = zeros(1,3);
b_H_m       = zeros(1,3);
b_total_m   = zeros(1,3);
tfd_L_m     = zeros(1,3);
tfe_L_m     = zeros(1,3);
tf_L_m      = zeros(1,3);
tf_H_m      = zeros(1,3);

j = 1 ; %preference parameter 匡跑计
for i = 1:3
x0 = [b_L(i), b_L(i), b_H(i), R(i), wf_L(i), wf_H(i)] ;
N = 1000000 ;
options = optimset('MaxFunEvals',N,'MaxIter',N);
[x,fval_m]=fsolve(@(x)myfun_4_c(x,pb(i),pe(i),am(i),c_(i),k0(i),k1(i),k2(i),wf_L(i),wf_H(i),alpha(j),gamma(j)),x0,options);
    bd_L_m(i)    = x(1) ;
    be_L_m(i)    = x(2) ;
    b_L_m(i)     = R(i)*bd_L_m(i)+(1-R(i))*be_L_m(i) ;
    b_H_m(i)     = x(3) ;
    b_total_m(i) = N_L(i)*(b_L_m(i))+N_H(i)*b_H_m(i) ;
    tfd_L_m(i)   = (1-k0(i)-((gamma(j)*pb(i)*bd_L_m(i)+(gamma(j)*k2(i)+alpha(j))*wf_L(i)*k1(i)*(bd_L_m(i))^(k2(i)))/(alpha(j)*wf_L(i))));
    tfe_L_m(i)   = (1-k0(i)-((gamma(j)*(pb(i)+pe(i))*be_L_m(i)+(gamma(j)*k2(i)+alpha(j))*wf_L(i)*k1(i)*(be_L_m(i))^(k2(i)))/(alpha(j)*wf_L(i))));
    tf_L_m(i)    = R(i)*tfd_L_m(i)+(1-R(i))*tfe_L_m(i) ;
    tf_H_m(i)    = (1-k0(i)-((gamma(j)*(pb(i)+pe(i))*b_H_m(i)+(gamma(j)*k2(i)+alpha(j))*wf_H(i)*k1(i)*(b_H_m(i))^(k2(i)))/(alpha(j)*wf_H(i))));
end

disp('     eta(1)      eta(2)     eta(3)')
disp(eta)
disp('     phi(1)      phi(2)     phi(3)')
disp(phi)
disp('    beta(1)   beta(2)   beta(3)')
disp(beta)
disp('    R_m(1)    R_m(2)   R_m(3)')
disp(R_m)