%Master Thesis R06323052 朝霍 
%calibration
clear
clc

%home production
t_sleep         = [513 522 527];    %何痸(と何) from Data_Time(new)
t_sleep_shower  = [564 574 578];    %何痸(と何)+縟瑍∟疍帝杆 from Data_Time(new)
t_sleep_eat     = [591 600 605];    %何痸(と何)+ノ繺(甦) from Data_Time(new)
t_sleep_shower_eat = [642 652 656]; %ゲ璶丁: 何痸(と何)+縟瑍∟疍帝杆┪てЗ+ノ繺(甦) from Data_Time(new)
t_ness = t_sleep_shower_eat ; 

t_free_work = 24-t_ness/60 ;  % パ丁+丁
k0 = zeros(1,3) ; 
k1 = zeros(1,3) ; 
k2 = zeros(1,3) ; 
x0 = [194 180 167]; % from Data_Woman_new
x1 = [323 313 295]; % from Data_Woman_new
x2 = [382 366 339]; % from Data_Woman_new
x0 = x0/60./t_free_work;
x1 = x1/60./t_free_work;
x2 = x2/60./t_free_work;
for i = 1:3
  k0(i) = x0(i);
  k1(i) = x1(i)-k0(i);
  f=@(kk)(k0(i)+k1(i)*2^kk-x2(i));
  [kk, fval_k2]=fsolve(f,0.4845);
  k2(i) = kk;  
end


%price parameter
%part 1 original(including 1822) 
pb_o  = [0.0859 0.1085 0.1581]; 
pe_o  = [0.0651 0.0652 0.0673];
c_ao  = [0.5714 0.6774 0.7497]; % a: average(だ糷舦) ((ネ秨綪-场だ)/Θ*2)
c_to  = [0.6314 0.7393 0.8049]; % t: total(场キА)   ((ネ秨綪-场だ)/Θ*2)
c_aom = [0.5788 0.7046 0.7787]; % a: average(だ糷舦) (度σ納(<23烦)+ひヾ產畑)
c_tom = [0.6535 0.7796 0.8465]; % t: total(场キА)   (度σ納(<23烦)+ひヾ產畑)
%part 2 new(excluding 1822)
pb_n =  [0.0711 0.0957 0.1368];
pe_an = [0.0868 0.0965 0.1047]; % a: average(だ糷舦)
c_an  = [0.5773 0.6711 0.7554]; % a: average(だ糷舦)
pe_tn = [0.0920 0.1020 0.1098]; % t: total(场キА)
c_tn  = [0.6341 0.7312 0.8108]; % t: total(场キА)

pb = pb_n ;      %from 基闽跑计
pe = pe_tn ;      %from 基闽跑计
c_ = c_tn ;     %from 基闽跑计

%production side parameter
n_L = [5749031 5731445 5053088]; % from キАネ▅羆计
n_H = [ 693625 1006389 1428377]; % from キАネ▅羆计
N_L = n_L./(n_L+n_H) ; %夹非てN_L
N_H = n_H./(n_L+n_H) ; %夹非てN_H
%ずネ跑计
b_L      = [2.333 2.189 1.874];  %from キАネ▅羆计(64)
b_H      = [1.691 1.599 1.636];  %from キАネ▅羆计(64)
b_total  = [2.309 2.154 1.923];  %from キАネ▅羆计(64)
wf_L     = [0.4758 0.5625 0.6490];  %from calibration 基闽跑计
wf_H     = [0.9429 1.0197 1.1240];  %from calibration 基闽跑计
am       = [1.0000 1.1461 1.3091];  %from calibration 基闽跑计

R        = [0.435 0.432 0.431];  %from PSFD (L->L percentage)
tf_L     = [134 161 185]; %from Data_Woman_work (min./day)
tf_H     = [228 236 245]; %from Data_Woman_work (min./day)
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
R_m         = zeros(1,3);
b_total_m   = zeros(1,3);
tfd_L_m     = zeros(1,3);
tfe_L_m     = zeros(1,3);
tf_L_m      = zeros(1,3);
tf_H_m      = zeros(1,3);
Z_L_m       = zeros(1,3);
Z_H_m       = zeros(1,3);

alpha  = zeros(1,3);
beta   = zeros(1,3);
gamma  = zeros(1,3);
eta    = zeros(1,3);
A      = zeros(1,3);
pi     = zeros(1,3);
for i = 1:3
x0_t = [b_L(i) b_L(i) b_H(i) R(i) 0.5438 0.8351 0.3459 1.4195 0.5349];
N = 100000 ;
options = optimset('MaxFunEvals', N, 'MaxIter', N, 'Algorithm', 'levenberg-marquardt');
[x,fval_cali]=fsolve(@(x)myfun_n(x,pb(i),pe(i),am(i),c_(i),k0(i),k1(i),k2(i),wf_L(i),wf_H(i),tf_L(i),tf_H(i),b_L(i),b_H(i),N_L(i),N_H(i),rho),x0_t,options);
    bd_L_m(i)    = x(1) ;
    be_L_m(i)    = x(2) ;
    R_m(i)       = x(3) ;
    b_L_m(i)     = R_m(i)*bd_L_m(i)+(1-R_m(i))*be_L_m(i) ;
    b_total_m(i) = N_L(i)*(b_L_m(i))+N_H(i)*b_H(i) ;
    gamma(i)     = x(4) ;
    alpha(i)     = x(5) ;
    eta(i)       = x(6) ;
    beta(i)      = x(7) ;
    A(i)         = x(8) ; 
    pi(i)        = x(9) ; 
    tfd_L_m(i)   = (1-k0(i)-((gamma(i)*pb(i)*bd_L_m(i)+(gamma(i)*k2(i)+alpha(i))*wf_L(i)*k1(i)*bd_L_m(i)^k2(i))/(alpha(i)*wf_L(i))));
    tfe_L_m(i)   = (1-k0(i)-((gamma(i)*(pb(i)+pe(i))*be_L_m(i)+(gamma(i)*k2(i)+alpha(i))*wf_L(i)*k1(i)*be_L_m(i)^k2(i))/(alpha(i)*wf_L(i))));
    tf_L_m(i)    = R_m(i)*tfd_L_m(i)+(1-R_m(i))*tfe_L_m(i) ;
    tf_H_m(i)    = (1-k0(i)-((gamma(i)*(pb(i)+pe(i))*b_H(i)+(gamma(i)*k2(i)+alpha(i))*wf_H(i)*k1(i)*b_H(i)^k2(i))/(alpha(i)*wf_H(i))));
    Z_L_m(i)     = N_L(i)*tf_L_m(i);
    Z_H_m(i)     = N_H(i)*tf_H_m(i);
calibration_P    = [alpha(i) eta(i) beta(i) gamma(i) A(i) pi(i) R_m(i)];
disp('     alpha        eta         beta     gamma        A            pi          R_m')
disp(calibration_P )
end
%calibration_P    = [alpha(1) eta(1) beta(1) gamma(1) A(1) pi(1) R_m(1)];
%%disp('       alpha          eta        beta     gamma        A            pi          R_m')
%%disp(calibration_P )

%%
%caliabration_check
bd_L_m      = zeros(1,3);
be_L_m      = zeros(1,3);
b_L_m       = zeros(1,3);
b_H_m       = zeros(1,3);
b_total_m   = zeros(1,3);
wf_L_m      = zeros(1,3);
wf_H_m      = zeros(1,3);
R_m         = zeros(1,3);
tfd_L_m     = zeros(1,3);
tfe_L_m     = zeros(1,3);
tf_L_m      = zeros(1,3);
tf_H_m      = zeros(1,3);
Z_L_m       = zeros(1,3);
Z_H_m       = zeros(1,3);

j = 1 ; %preference parameter 匡跑计
for i = 1:3
x0 = [b_L(i), b_L(i), b_H(i), R(i), wf_L(i), wf_H(i)] ;
N = 1000000 ;
options = optimset('MaxFunEvals',N,'MaxIter',N, 'Algorithm', 'levenberg-marquardt');
[x,fval_m]=fsolve(@(x)myfun_n_c(x,pb(i),pe(i),am(i),c_(i),k0(i),k1(i),k2(i),N_L(i),N_H(i),rho,alpha(j),eta(j),gamma(j),beta(j),A(j),pi(j)),x0,options);

    bd_L_m(i)    = x(1) ;
    be_L_m(i)    = x(2) ;
    b_L_m(i)     = x(4)*bd_L_m(i)+(1-x(4))*be_L_m(i) ;
    b_H_m(i)     = x(3) ;
    R_m(i)       = x(4) ;
    b_total_m(i) = N_L(i)*(b_L_m(i))+N_H(i)*b_H_m(i) ;
    wf_L_m(i)    = x(5) ;
    wf_H_m(i)    = x(6) ;
    tfd_L_m(i)   = (1-k0(i)-((gamma(i)*pb(i)*bd_L_m(i)+(gamma(i)*k2(i)+alpha(i))*wf_L(i)*k1(i)*bd_L_m(i)^k2(i))/(alpha(i)*wf_L(i))));
    tfe_L_m(i)   = (1-k0(i)-((gamma(i)*(pb(i)+pe(i))*be_L_m(i)+(gamma(i)*k2(i)+alpha(i))*wf_L(i)*k1(i)*be_L_m(i)^k2(i))/(alpha(i)*wf_L(i))));
    tf_L_m(i)    = R_m(i)*tfd_L_m(i)+(1-R_m(i))*tfe_L_m(i) ;
    tf_H_m(i)    = (1-k0(i)-((gamma(i)*(pb(i)+pe(i))*b_H_m(i)+(gamma(i)*k2(i)+alpha(i))*wf_H(i)*k1(i)*b_H_m(i)^k2(i))/(alpha(i)*wf_H(i))));
    Z_L_m(i)     = N_L(i)*tf_L_m(i);
    Z_H_m(i)     = N_H(i)*tf_H_m(i);
end

disp('    R_m(1)    R_m(2)   R_m(3)')
disp(R_m)