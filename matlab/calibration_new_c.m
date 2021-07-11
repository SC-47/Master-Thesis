%Master Thesis R06323052 朝霍 
%calibration 
clear
clc

%home production
t_sleep         = [513 522 527];    %何痸(と何) from Data_Time(new)
t_sleep_shower  = [564 574 578];    %何痸(と何)+縟瑍∟疍帝杆 from Data_Time(new)
t_sleep_eat     = [591 600 605];    %何痸(と何)+ノ繺(甦) from Data_Time(new)
t_sleep_shower_eat = [642 652 656]; %ゲ璶丁: 何痸(と何)+縟瑍∟疍帝杆┪てЗ+ノ繺(甦) from Data_Time(new)
t_ness = t_sleep_shower ; 

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
pb_o  = [0.0921 0.1225 0.1737]; 
pe_o  = [0.0676 0.0698 0.0667];
c_ao  = [0.6137 0.7217 0.8182]; % a: average(だ糷舦) ((ネ秨綪-场だ)/Θ*2)
c_to  = [0.6876 0.7986 0.8889]; % t: total(场キА)   ((ネ秨綪-场だ)/Θ*2)
c_aom = [0.6256 0.7528 0.8537]; % a: average(だ糷舦) (度σ納(<23烦)+ひヾ產畑)
c_tom = [0.7162 0.8468 0.9406]; % t: total(场キА)   (度σ納(<23烦)+ひヾ產畑)
%part 2 new(excluding 1822)
pb_n =  [0.0763 0.1062 0.1504];
pe_an = [0.0903 0.1016 0.1120]; % a: average(だ糷舦)
c_an  = [0.6206 0.7196 0.8191]; % a: average(だ糷舦)
pe_tn = [0.0967 0.1084 0.1184]; % t: total(场キА)
c_tn  = [0.6906 0.7945 0.8902]; % t: total(场キА)

pb = pb_o ;      %from 基闽跑计
pe = pe_o ;      %from 基闽跑计
c_ = c_ao ;     %from 基闽跑计

%production side parameter
n_L = [5749031 5731445 5053088]; % from キАネ▅羆计
n_H = [ 693625 1006389 1428377]; % from キАネ▅羆计
N_L = n_L./(n_L+n_H) ; %夹非てN_L
N_H = n_H./(n_L+n_H) ; %夹非てN_H
%ずネ跑计
b_L      = [2.333 2.189 1.874];  %from キАネ▅羆计(64)
b_H      = [1.691 1.599 1.636];  %from キАネ▅羆计(64)
b_total  = [2.309 2.154 1.923];  %from キАネ▅羆计(64)
wf_L     = [0.5182 0.5990 0.6938];  %from calibration 基闽跑计
wf_H     = [1.0034 1.1081 1.2063];  %from calibration 基闽跑计
am       = [1.0000 1.1527 1.3158];  %from calibration 基闽跑计

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

phi    = zeros(1,3);
eta    = zeros(1,3);
beta   = zeros(1,3);
gamma  = zeros(1,3);
A      = zeros(1,3);
pi     = zeros(1,3);
for i = 1:3
x0_t = [b_L(i) b_L(i) b_H(i) R(i) 1.1418 1.3558 0.7686 1.1885 0.5454];
N = 100000 ;
options = optimset('MaxFunEvals', N, 'MaxIter', N, 'Algorithm', 'levenberg-marquardt');
[x,fval_cali]=fsolve(@(x)myfun(x,pb(i),pe(i),am(i),c_(i),k0(i),k1(i),k2(i),wf_L(i),wf_H(i),tf_L(i),tf_H(i),b_L(i),b_H(i),N_L(i),N_H(i),rho),x0_t,options);
    bd_L_m(i)    = x(1) ;
    be_L_m(i)    = x(2) ;
    R_m(i)       = x(3) ;
    b_L_m(i)     = R_m(i)*bd_L_m(i)+(1-R_m(i))*be_L_m(i) ;
    b_total_m(i) = N_L(i)*(b_L_m(i))+N_H(i)*b_H(i) ;
    gamma(i)     = x(4) ;
    phi(i)       = x(5) ;
    eta(i)       = x(6) ;
    A(i)         = x(7) ; 
    pi(i)        = x(8) ; 
    tfd_L_m(i)   = (1-k0(i)-((gamma(i)*pb(i)*bd_L_m(i)+(gamma(i)*k2(i)+(phi(i)*bd_L_m(i)^(1-eta(i))))*wf_L(i)*k1(i)*(bd_L_m(i))^(k2(i)))/(wf_L(i)*(phi(i)*bd_L_m(i)^(1-eta(i)))))) ;
    tfe_L_m(i)   = (1-k0(i)-((gamma(i)*(pb(i)+pe(i))*be_L_m(i)+(gamma(i)*k2(i)+(phi(i)*be_L_m(i)^(1-eta(i))))*wf_L(i)*k1(i)*(be_L_m(i))^(k2(i)))/(wf_L(i)*(phi(i)*be_L_m(i)^(1-eta(i)))))) ;
    tf_L_m(i)    = R_m(i)*tfd_L_m(i)+(1-R_m(i))*tfe_L_m(i) ;
    tf_H_m(i)    = (1-k0(i)-((gamma(i)*(pb(i)+pe(i))*b_H(i)+(gamma(i)*k2(i)+(phi(i)*b_H(i)^(1-eta(i))))*wf_H(i)*k1(i)*(b_H(i))^(k2(i)))/(wf_H(i)*(phi(i)*b_H(i)^(1-eta(i)))))) ;
    Z_L_m(i)     = N_L(i)*tf_L_m(i);
    Z_H_m(i)     = N_H(i)*tf_H_m(i);
    beta(i)      = -((1+gamma(i))*log(((be_L_m(i)/bd_L_m(i))^(1-eta(i)))*(pb(i)*bd_L_m(i)+wf_L(i)*k1(i)*k2(i)*(bd_L_m(i))^(k2(i)))/((pb(i)+pe(i))*be_L_m(i)+wf_L(i)*k1(i)*k2(i)*(be_L_m(i))^(k2(i))))+(phi(i)/(1-eta(i)))*(bd_L_m(i)^(1-eta(i))-be_L_m(i)^(1-eta(i))))/log(wf_L(i)/wf_H(i)) ;
%%calibration_P    = [phi(i) eta(i) beta(i) gamma(i) A(i) pi(i) R_m(i)];
%%disp('       phi          eta        beta     gamma        A            pi          R_m')
%%disp(calibration_P )
end
k = 1;
calibration_P    = [phi(k) eta(k) beta(k) gamma(k) A(k) pi(k) R_m(k)];
disp('       phi          eta        beta     gamma        A            pi          R_m')
disp(calibration_P )
original = [bd_L_m(k) be_L_m(k) b_L_m(k) b_H(k) b_total_m(k) R_m(k)];
disp('      bd_L       be_L        b_L         b_H      b_total        R')
disp(original)


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
for i = 2:2
x0 = [b_L(i), b_L(i), b_H(i), R(i), wf_L(i), wf_H(i)] ;
N = 1000000 ;
options = optimset('MaxFunEvals',N,'MaxIter',N, 'Algorithm', 'levenberg-marquardt');
[x,fval_m]=fsolve(@(x)myfun_3_c(x,pb(i),pe(i),am(i),c_(i),k0(1),k1(1),k2(1),N_L(i),N_H(i),rho,phi(j),eta(j),gamma(j),beta(j),A(j),pi(j)),x0,options);
    bd_L_m(i)    = x(1) ;
    be_L_m(i)    = x(2) ;
    b_L_m(i)     = x(4)*bd_L_m(i)+(1-x(4))*be_L_m(i) ;
    b_H_m(i)     = x(3) ;
    R_m(i)       = x(4) ;
    b_total_m(i) = N_L(i)*(b_L_m(i))+N_H(i)*b_H_m(i) ;
    wf_L_m(i)    = x(5) ;
    wf_H_m(i)    = x(6) ;
    tfd_L_m(i)   = (1-k0(i)-((gamma(j)*pb(i)*bd_L_m(i)+(gamma(j)*k2(i)+(phi(j)*bd_L_m(i)^(1-eta(j))))*wf_L(i)*k1(i)*(bd_L_m(i))^(k2(i)))/(wf_L_m(i)*(phi(j)*bd_L_m(i)^(1-eta(j))))));
    tfe_L_m(i)   = (1-k0(i)-((gamma(j)*(pb(i)+pe(i))*be_L_m(i)+(gamma(j)*k2(i)+(phi(j)*be_L_m(i)^(1-eta(j))))*wf_L_m(i)*k1(i)*(be_L_m(i))^(k2(i)))/(wf_L_m(i)*(phi(j)*be_L_m(i)^(1-eta(j))))));
    tf_L_m(i)    = R_m(i)*tfd_L_m(i)+(1-R_m(i))*tfe_L_m(i);
    tf_H_m(i)    = (1-k0(i)-((gamma(j)*(pb(i)+pe(i))*b_H_m(i)+(gamma(j)*k2(i)+(phi(j)*b_H_m(i)^(1-eta(j))))*wf_H_m(i)*k1(i)*(b_H_m(i))^(k2(i)))/(wf_H(i)*(phi(j)*b_H_m(i)^(1-eta(j))))));
    Z_L_m(i)     = N_L(i)*tf_L_m(i);
    Z_H_m(i)     = N_H(i)*tf_H_m(i);
    decomposition = [bd_L_m(i) be_L_m(i) b_L_m(i) b_H(i) b_total_m(i) R_m(i)];
    disp('      bd_L       be_L        b_L          b_H      b_total        R')
    disp(decomposition)
end