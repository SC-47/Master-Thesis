%Master Thesis R06323052 陳思齊 
%calibration (wage deifinition is only wage(excluding owner))
clear
clc

%price parameter

%1 definition(受雇人員報酬)
am_1    = [1.0000 1.1527 1.3158];  %from calibration 價格相關變數(受雇人員報酬)
wf_L_1  = [0.5182 0.5990 0.6938];  %from calibration 價格相關變數(受雇人員報酬)
wf_H_1  = [1.0034 1.1081 1.2063];  %from calibration 價格相關變數(受雇人員報酬)
%original(including 1822) 
pb_1_o  = [0.0921 0.1225 0.1737]; 
pe_1_o  = [0.0676 0.0698 0.0667];
c_1_ao  = [0.6137 0.7217 0.8182]; % a: average(分層加權) ((生活開銷-子女部分)/成年人*2)
c_1_to  = [0.6876 0.7986 0.8889]; % t: total(全部平均)   ((生活開銷-子女部分)/成年人*2)
c_1_aom = [0.6256 0.7528 0.8537]; % a: average(分層加權) (僅考慮子女(<23歲)+夫妻家庭)
c_1_tom = [0.7162 0.8468 0.9406]; % t: total(全部平均)   (僅考慮子女(<23歲)+夫妻家庭)
%new(excluding 1822)
pb_1_n  = [0.0763 0.1062 0.1504];
pe_1_an = [0.0903 0.1016 0.1120]; % a: average(分層加權)
c_1_an  = [0.6206 0.7196 0.8191]; % a: average(分層加權)
pe_1_tn = [0.0967 0.1084 0.1184]; % t: total(全部平均)
c_1_tn  = [0.6906 0.7945 0.8902]; % t: total(全部平均)

%2 definition(受雇人員報酬+產業主所得)
am_2   = [1.0000 1.1461 1.3091];  %from calibration 價格相關變數(受雇人員報酬+產業主所得)
wf_L_2 = [0.4758 0.5625 0.6490];  %from calibration 價格相關變數(受雇人員報酬+產業主所得)
wf_H_2 = [0.9429 1.0197 1.1240];  %from calibration 價格相關變數(受雇人員報酬+產業主所得)
%original(including 1822) 
pb_2_o  = [0.0859 0.1085 0.1581]; 
pe_2_o  = [0.0622 0.0633 0.0660];
c_2_ao  = [0.5714 0.6774 0.7497]; % a: average(分層加權) ((生活開銷-子女部分)/成年人*2)
c_2_to  = [0.6314 0.7393 0.8049]; % t: total(全部平均)   ((生活開銷-子女部分)/成年人*2)
c_2_aom = [0.5788 0.7046 0.7787]; % a: average(分層加權) (僅考慮子女(<23歲)+夫妻家庭)
c_2_tom = [0.6535 0.7796 0.8465]; % t: total(全部平均)   (僅考慮子女(<23歲)+夫妻家庭)
%new(excluding 1822)
pb_2_n  = [0.0711 0.0957 0.1368];
pe_2_an = [0.0839 0.0946 0.1035]; % a: average(分層加權)
c_2_an  = [0.5773 0.6711 0.7554]; % a: average(分層加權)
pe_2_tn = [0.0891 0.1000 0.1085]; % t: total(全部平均)
c_2_tn  = [0.6341 0.7312 0.8108]; % t: total(全部平均)

%control part
am   =    am_2 ; 
wf_L =  wf_L_2 ; 
wf_H =  wf_H_2 ;
pb   = pb_2_n ;     
pe   = pe_2_tn ;    
c_   = c_2_tn ;    


%home production
t_sleep            = [513 522 527]; %睡眠(含午睡) from Data_Time(new)
t_sleep_shower     = [564 574 578]; %睡眠(含午睡)+盥洗、沐浴、著裝 from Data_Time(new)
t_sleep_eat        = [591 600 605]; %睡眠(含午睡)+用餐(含吃宵夜) from Data_Time(new)
t_sleep_shower_eat = [642 652 656]; %睡眠(含午睡)+盥洗、沐浴、著裝或化妝+用餐(含吃宵夜)(必要時間) from Data_Time(new)
t_ness = t_sleep_eat ; 

t_free_work = 24-t_ness/60 ;  % 自由時間+工作時間
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


%preference parameter

%fertility related
b_L      = [2.333 2.189 1.874];  %from 各世代平均生育總數(64)
b_H      = [1.691 1.599 1.636];  %from 各世代平均生育總數(64)
b_total  = [2.309 2.154 1.923];  %from 各世代平均生育總數(64)
R        = [0.435 0.432 0.431];  %from PSFD (L->L percentage)

%production side parameter
tf_L = [134 161 185]; %from Data_Woman_work (min./day)
tf_H = [228 236 245]; %from Data_Woman_work (min./day)
tf_L = tf_L/60./t_free_work ; %標準化(換算成小時後，在換算成自由時間占比)
tf_H = tf_H/60./t_free_work ; %標準化(換算成小時後，在換算成自由時間占比)
n_L  = [5749031 5731445 5053088]; % from 各世代平均生育總數
n_H  = [ 693625 1006389 1428377]; % from 各世代平均生育總數
N_L  = n_L./(n_L+n_H) ; %標準化N_L
N_H  = n_H./(n_L+n_H) ; %標準化N_H
Z_L  = N_L.*tf_L ; %L_type labor force
Z_H  = N_H.*tf_H ; %H_type labor force

% Assume
rho = 0.401 ; % from Krusell(2000)



%calibration for gamma phi eta  & part of beta A pi
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
beta        = zeros(1,3);
A           = zeros(1,3);
pi          = zeros(1,3);

for i = 1:1
x0_t = [b_L(i) b_L(i) b_H(i) R(i) 1.1418 1.3558 0.7686 1.1885 0.5454]; % bd_L be_L b_H R gamma phi eta A pi
N = 100000 ;
options = optimset('MaxFunEvals', N, 'MaxIter', N, 'Algorithm', 'levenberg-marquardt');
[x,fval_cali]=fsolve(@(x)myfun_2a(x,pb(i),pe(i),am(i),c_(i),k0(i),k1(i),k2(i),wf_L(i),wf_H(i),tf_L(i),tf_H(i),b_L(i),b_H(i),N_L(i),N_H(i),rho),x0_t,options);
    bd_L_m(i)    = x(1) ;
    be_L_m(i)    = x(2) ;
    R_m(i)       = x(3) ;
    b_L_m(i)     = R_m(i)*bd_L_m(i)+(1-R_m(i))*be_L_m(i) ;
    b_total_m(i) = N_L(i)*b_L_m(i)+N_H(i)*b_H(i) ;
    gamma        = x(4) ;
    alpha_L      = x(5) ;
    alpha_H      = x(6) ;
    A(i)         = x(7) ; 
    pi(i)        = x(8) ; 
    tfd_L_m(i)   = (1-k0(i)-((gamma*pb(i)*bd_L_m(i)+(gamma*k2(i)+alpha_L)*wf_L(i)*k1(i)*bd_L_m(i)^k2(i))/(alpha_L*wf_L(i)))) ;
    tfe_L_m(i)   = (1-k0(i)-((gamma*(pb(i)+pe(i))*be_L_m(i)+(gamma*k2(i)+alpha_L)*wf_L(i)*k1(i)*be_L_m(i)^k2(i))/(alpha_L*wf_L(i)))) ;
    tf_L_m(i)    = R_m(i)*tfd_L_m(i)+(1-R_m(i))*tfe_L_m(i) ;
    tf_H_m(i)    = (1-k0(i)-((gamma*(pb(i)+pe(i))*b_H(i)+(gamma*k2(i)+alpha_H)*wf_H(i)*k1(i)*b_H(i)^k2(i))/(alpha_H*wf_H(i)))) ;
    Z_L_m(i)     = N_L(i)*tf_L_m(i) ;
    Z_H_m(i)     = N_H(i)*tf_H_m(i) ;
    beta(i)      = (-((1+gamma)*log((pb(i)*bd_L_m(i)+wf_L(i)*k1(i)*k2(i)*(bd_L_m(i))^(k2(i)))/((pb(i)+pe(i))*be_L_m(i)+wf_L(i)*k1(i)*k2(i)*(be_L_m(i))^(k2(i))))+alpha_L*log(bd_L_m(i)/be_L_m(i)))/(log(wf_L(i)/wf_H(i)))) ;
calibration_P1   = [alpha_L alpha_H gamma];
disp('     alp_L      alp_H   gamma')
disp(calibration_P1 )
calibration_P2   = [beta(i) A(i) pi(i)];
disp('      beta          A            pi')
disp(calibration_P2 )
disp('      R_m')
disp(R_m(i))
end



%caliabration for part of beta A pi & part of pb pe c_
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

beta   = zeros(1,3);
A      = zeros(1,3);
pi     = zeros(1,3);
pb_m   = zeros(1,3);
pe_m   = zeros(1,3);
c_m    = zeros(1,3);

for i = 1:2
    if i == 2
        k0(2)=k0(2)*1.03;
        k1(2)=k1(2)*1.03;
    end
x0_t = [b_L(i) b_L(i) R(i) pb(i) pe(i) c_(i) 1.1885 0.5454 ];  
N = 100000 ;
options = optimset('MaxFunEvals', N, 'MaxIter', N, 'Algorithm', 'levenberg-marquardt');
[x,fval_cali]=fsolve(@(x)myfun_2a2c(x,alpha_L,alpha_H,am(i),gamma,k0(i),k1(i),k2(i),wf_L(i),wf_H(i),b_L(i),b_H(i),N_L(i),N_H(i),rho,tf_L(i),tf_H(i)),x0_t,options);
    bd_L_m(i)    = x(1) ;
    be_L_m(i)    = x(2) ;
    R_m(i)       = x(3) ;
    b_L_m(i)     = R_m(i)*bd_L_m(i)+(1-R_m(i))*be_L_m(i) ;
    b_total_m(i) = N_L(i)*(b_L_m(i))+N_H(i)*b_H(i) ;
    pb_m(i)      = x(4) ;
    pe_m(i)      = x(5) ;
    c_m(i)       = x(6) ;
    A(i)         = x(7) ; 
    pi(i)        = x(8) ; 
    tfd_L_m(i)   = (1-k0(i)-((gamma*pb_m(i)*bd_L_m(i)+(gamma*k2(i)+alpha_L)*wf_L(i)*k1(i)*bd_L_m(i)^k2(i))/(alpha_L*wf_L(i)))) ;
    tfe_L_m(i)   = (1-k0(i)-((gamma*(pb_m(i)+pe_m(i))*be_L_m(i)+(gamma*k2(i)+alpha_L)*wf_L(i)*k1(i)*be_L_m(i)^k2(i))/(alpha_L*wf_L(i)))) ;
    tf_L_m(i)    = R_m(i)*tfd_L_m(i)+(1-R_m(i))*tfe_L_m(i) ;
    tf_H_m(i)    = (1-k0(i)-((gamma*(pb_m(i)+pe_m(i))*b_H(i)+(gamma*k2(i)+alpha_H)*wf_H(i)*k1(i)*b_H(i)^k2(i))/(alpha_H*wf_H(i)))) ;
    Z_L_m(i)     = N_L(i)*tf_L_m(i) ;
    Z_H_m(i)     = N_H(i)*tf_H_m(i) ;
    beta(i)      = (-((1+gamma)*log((pb_m(i)*bd_L_m(i)+wf_L(i)*k1(i)*k2(i)*(bd_L_m(i))^(k2(i)))/((pb_m(i)+pe_m(i))*be_L_m(i)+wf_L(i)*k1(i)*k2(i)*(be_L_m(i))^(k2(i))))+alpha_L*log(bd_L_m(i)/be_L_m(i)))/(log(wf_L(i)/wf_H(i)))) ;
end
calibration_P2_c   = [beta(1) A(1) pi(1)];
disp('      beta          A            pi')
disp(calibration_P2_c )
disp('      R_m        R_m       R_m')
disp(R_m)
disp('     pb_m     pb_m      pb_m')
disp(pb_m)
disp('     pe_m     pe_m      pe_m')
disp(pe_m)

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

for i = 1:3
xi = [b_L(i), b_L(i), b_H(i), R(i), wf_L(i), wf_H(i)] ;
N = 1000000 ;
options = optimset('MaxFunEvals',N,'MaxIter',N, 'Algorithm', 'levenberg-marquardt');
% j == >  preference parameter 選取變數(for phi eta gamma)
[x,fval_m]=fsolve(@(x)myfun_2a_c(x,pb_m(i),pe_m(i),am(i),c_m(i),k0(i),k1(i),k2(i),N_L(i),N_H(i),rho,alpha_L,alpha_H,gamma,beta(i),A(i),pi(i)),xi,options);
    bd_L_m(i)    = x(1) ;
    be_L_m(i)    = x(2) ;
    b_L_m(i)     = x(4)*bd_L_m(i)+(1-x(4))*be_L_m(i) ;
    b_H_m(i)     = x(3) ;
    R_m(i)       = x(4) ;
    b_total_m(i) = N_L(i)*(b_L_m(i))+N_H(i)*b_H_m(i) ;
    wf_L_m(i)    = x(5) ;
    wf_H_m(i)    = x(6) ;
    tfd_L_m(i)   = (1-k0(i)-((gamma*pb_m(i)*bd_L_m(i)+(gamma*k2(i)+alpha_L)*wf_L_m(i)*k1(i)*bd_L_m(i)^k2(i))/(alpha_L*wf_L_m(i)))) ;
    tfe_L_m(i)   = (1-k0(i)-((gamma*(pb_m(i)+pe_m(i))*be_L_m(i)+(gamma*k2(i)+alpha_L)*wf_L_m(i)*k1(i)*be_L_m(i)^k2(i))/(alpha_L*wf_L_m(i)))) ;
    tf_L_m(i)    = R_m(i)*tfd_L_m(i)+(1-R_m(i))*tfe_L_m(i);
    tf_H_m(i)    = (1-k0(i)-((gamma*(pb_m(i)+pe_m(i))*b_H_m(i)+(gamma*k2(i)+alpha_H)*wf_H_m(i)*k1(i)*b_H_m(i)^k2(i))/(alpha_H*wf_H_m(i)))) ;
    Z_L_m(i)     = N_L(i)*tf_L_m(i);
    Z_H_m(i)     = N_H(i)*tf_H_m(i);
end
disp('      beta        beta       beta')
disp(beta)
disp('      R_m        R_m        R_m')
disp(R_m)