    %Master Thesis R06323052 陳思齊 
%calibration (wage deifinition is only wage(excluding owner))
clear
clc

%price parameter


%price parameter

%1 definition(受雇人員報酬)
am_1      = [1.0000 1.2835];  %from calibration 價格相關變數(受雇人員報酬)
wf_L_1    = [0.5235 0.6806];  %from calibration 價格相關變數(受雇人員報酬)
wf_H_1    = [1.0031 1.1490];  %from calibration 價格相關變數(受雇人員報酬)
am_1n     = [1.0000 1.0000*1.1357*1.1251];  %from calibration wage_growing(受雇人員報酬)
wf_L_1n   = [0.5235 0.5235*1.1478*1.1313];  %from calibration wage_growing(受雇人員報酬)(薪資以學歷區分各別算成長率)
wf_H_1n   = [1.0031 1.0031*1.1194*1.0993];  %from calibration wage_growing(受雇人員報酬)(薪資以學歷區分各別算成長率)
%original(including 1822) 
pb_1_o    = 0.0730 ; 
pe_1_o    = 0.0549 ;
c_t_1_ao  = 0.6485 ; % a: average(分層加權) ((生活開銷-子女部分)/成年人*2)(total)
c_t_1_to  = 0.7510 ; % t: total(全部平均)   ((生活開銷-子女部分)/成年人*2)(total)
c_t_1_aom = 0.7752 ; % a: average(分層加權) (僅考慮子女(<23歲)+夫妻家庭)(total)
c_t_1_tom = 0.8805 ; % t: total(全部平均)   (僅考慮子女(<23歲)+夫妻家庭)(total)

%2 definition(受雇人員報酬+產業主所得)
am_2      = [1.0000 1.2772];  %from calibration 價格相關變數(受雇人員報酬+產業主所得)
wf_L_2    = [0.4849 0.6391];  %from calibration 價格相關變數(受雇人員報酬+產業主所得)
wf_H_2    = [0.9242 1.0671];  %from calibration 價格相關變數(受雇人員報酬+產業主所得)
am_2n     = [1.0000 1.0000*1.1357*1.1251];  %from calibration wage_growing(受雇人員報酬+產業主所得)
wf_L_2n   = [0.4849 0.4849*1.1478*1.1313];  %from calibration wage_growing(受雇人員報酬+產業主所得)(薪資以學歷區分各別算成長率)
wf_H_2n   = [0.9242 0.9242*1.1194*1.0993];  %from calibration wage_growing(受雇人員報酬+產業主所得)(薪資以學歷區分各別算成長率)
%original(including 1822) 
pb_2_o    = 0.0676 ; 
pe_2_o    = 0.0508 ;
c_t_2_ao  = 0.6060 ; % a: average(分層加權) ((生活開銷-子女部分)/成年人*2)(total)
c_t_2_to  = 0.6892 ; % t: total(全部平均)   ((生活開銷-子女部分)/成年人*2)(total)
c_t_2_aom = 0.7261 ; % a: average(分層加權) (僅考慮子女(<23歲)+夫妻家庭)(total)
c_t_2_tom = 0.8108 ; % t: total(全部平均)   (僅考慮子女(<23歲)+夫妻家庭)(total)


%3 definition(經常性收入)
am_3      = [1.0000 1.2917] ;  %from calibration 價格相關變數(經常性收入)
wf_L_3    = [0.4612 0.6219] ;  %from calibration 價格相關變數(經常性收入)
wf_H_3    = [0.8671 1.0157] ;  %from calibration 價格相關變數(經常性收入)
am_3n     = [1.0000 1.0000*1.1357*1.1251] ;  %from calibration wage_growing(經常性收入)
wf_L_3n   = [0.4612 0.4612*1.1478*1.1313] ;  %from calibration wage_growing(經常性收入)(薪資以學歷區分各別算成長率)
wf_H_3n   = [0.8671 0.8671*1.1194*1.0993] ;  %from calibration wage_growing(經常性收入)(薪資以學歷區分各別算成長率)
%original(including 1822) 
pb_3_o    = 0.0724 ; 
pe_3_o    = 0.0501 ;
c_t_3_ao  = 0.5235 ; % a: average(分層加權) ((生活開銷-子女部分)/成年人*2)(total)
c_t_3_to  = 0.5844 ; % t: total(全部平均)   ((生活開銷-子女部分)/成年人*2)(total)
c_t_3_aom = 0.6254 ; % a: average(分層加權) (僅考慮子女(<23歲)+夫妻家庭)(total)
c_t_3_tom = 0.6867 ; % t: total(全部平均)   (僅考慮子女(<23歲)+夫妻家庭)(total)




%control part
am   =     am_3n ; 
wf_L =   wf_L_3n ; 
wf_H =   wf_H_3n ;
pb   =   pb_3_o ;     
pe   =   pe_3_o ;    
c_   = c_t_3_aom ;    


%home production
t_sleep            = [516 533]; %睡眠(含午睡) from Data_Time(new)
t_sleep_eat        = [588 596]; %睡眠(含午睡)+用餐(含吃宵夜) from Data_Time(new)
t_sleep_shower_eat = [644 662]; %睡眠(含午睡)+盥洗、沐浴、著裝或化妝+用餐(含吃宵夜)(必要時間) from Data_Time(new)
t_ness = t_sleep_eat ; 

t_free_work = (24-t_ness/60) ;  % 自由時間+工作時間
k0 = zeros(1,2) ; 
k1 = zeros(1,2) ; 
k2 = zeros(1,2) ; 
x0 = [187 164]; % from Data_Woman_new
x1 = [316 286]; % from Data_Woman_new
x2 = [375 332]; % from Data_Woman_new
x0 = x0/60./t_free_work;
x1 = x1/60./t_free_work;
x2 = x2/60./t_free_work;
for i = 1:2
  k0(i) = x0(i);
  k1(i) = x1(i)-k0(i);
  f=@(kk)(k0(i)+k1(i)*2^kk-x2(i));
  [kk, fval_k2]=fsolve(f,0.4845);
  k2(i) = kk;  
end


%preference parameter

%fertility related
b_L      = [2.253 1.748];  %from 各世代平均生育總數(66)
b_H      = [1.465 1.402];  %from 各世代平均生育總數(66)
b_total  = [2.213 1.781];  %from 各世代平均生育總數(66)
R        = [0.565 0.569];  %from PSFD (L->H percentage)

%production side parameter
tf_L = [143 190]; %from Data_Woman_work (min./day)
tf_H = [236 255]; %from Data_Woman_work (min./day)
tf_L = tf_L/60./t_free_work ; %標準化(換算成小時後，在換算成自由時間占比)
tf_H = tf_H/60./t_free_work ; %標準化(換算成小時後，在換算成自由時間占比)
n_L  = [1137265 1000711]; % from 各世代平均生育總數
n_H  = [ 144652  316335]; % from 各世代平均生育總數
N_L  = n_L./(n_L+n_H) ; %標準化N_L
N_H  = n_H./(n_L+n_H) ; %標準化N_H
Z_L  = N_L.*tf_L ; %L_type labor force
Z_H  = N_H.*tf_H ; %H_type labor force

% Assume
rho =  0.401 ; % from Krusell(2000) 1.669
%rho =  0.296 ; % from Katz and Murphy (1992) 1.42
%rho = -0.921 ; % from assumption(critical point)
%rho = -2.030 ; % from assumption(stable lower bound)



%calibration for gamma phi eta  & part of beta A pi
bd_L_m      = zeros(1,2);
be_L_m      = zeros(1,2);
b_L_m       = zeros(1,2);
R_m         = zeros(1,2);
b_total_m   = zeros(1,2);
tfd_L_m     = zeros(1,2);
tfe_L_m     = zeros(1,2);
tf_L_m      = zeros(1,2);
tf_H_m      = zeros(1,2);
Z_L_m       = zeros(1,2);
Z_H_m       = zeros(1,2);
beta        = zeros(1,2);
A           = zeros(1,2);
pi          = zeros(1,2);

for i = 1
x0_t = [b_L(i) b_L(i) b_H(i) R(i) 1.1418 1.3558 0.7686 1.1885 0.5454]; % bd_L be_L b_H R gamma phi eta A pi
N = 100000 ;
options = optimset('MaxFunEvals', N, 'MaxIter', N, 'Algorithm', 'levenberg-marquardt');
[x,fval_cali]=fsolve(@(x)myfun_2a(x,pb(i),pe(i),am(i),c_(i),k0(i),k1(i),k2(i),wf_L(i),wf_H(i),tf_L(i),tf_H(i),b_L(i),b_H(i),N_L(i),N_H(i),rho),x0_t,options);
    bd_L_m(i)    = x(1) ;
    be_L_m(i)    = x(2) ;
    R_m(i)       = x(3) ;
    b_L_m(i)     = (1-R_m(i))*bd_L_m(i)+R_m(i)*be_L_m(i) ;
    b_total_m(i) = N_L(i)*b_L_m(i)+N_H(i)*b_H(i) ;
    gamma        = x(4) ;
    alpha_L      = x(5) ;
    alpha_H      = x(6) ;
    A(i)         = x(7) ; 
    pi(i)        = x(8) ; 
    tfd_L_m(i)   = (1-k0(i)-((gamma*pb(i)*bd_L_m(i)+(gamma*k2(i)+alpha_L)*wf_L(i)*k1(i)*bd_L_m(i)^k2(i))/(alpha_L*wf_L(i)))) ;
    tfe_L_m(i)   = (1-k0(i)-((gamma*(pb(i)+pe(i))*be_L_m(i)+(gamma*k2(i)+alpha_L)*wf_L(i)*k1(i)*be_L_m(i)^k2(i))/(alpha_L*wf_L(i)))) ;
    tf_L_m(i)    = (1-R_m(i))*tfd_L_m(i)+R_m(i)*tfe_L_m(i) ;
    tf_H_m(i)    = (1-k0(i)-((gamma*(pb(i)+pe(i))*b_H(i)+(gamma*k2(i)+alpha_H)*wf_H(i)*k1(i)*b_H(i)^k2(i))/(alpha_H*wf_H(i)))) ;
    Z_L_m(i)     = N_L(i)*tf_L_m(i) ;
    Z_H_m(i)     = N_H(i)*tf_H_m(i) ;
    beta(i)      = (-((1+gamma)*log((pb(i)*bd_L_m(i)+wf_L(i)*k1(i)*k2(i)*(bd_L_m(i))^(k2(i)))/((pb(i)+pe(i))*be_L_m(i)+wf_L(i)*k1(i)*k2(i)*(be_L_m(i))^(k2(i))))+alpha_L*log(bd_L_m(i)/be_L_m(i)))/(log(wf_L(i)/wf_H(i)))) ;
calibration_P1    = [alpha_L alpha_H gamma];
disp('      alp_L      alp_H    gamma')
disp(calibration_P1)
end



%caliabration for part of beta A pi & part of pb pe c_
bd_L_m      = zeros(1,2);
be_L_m      = zeros(1,2);
b_L_m       = zeros(1,2);
R_m         = zeros(1,2);
b_total_m   = zeros(1,2);
tfd_L_m     = zeros(1,2);
tfe_L_m     = zeros(1,2);
tf_L_m      = zeros(1,2);
tf_H_m      = zeros(1,2);
Z_L_m       = zeros(1,2);
Z_H_m       = zeros(1,2);

beta   = zeros(1,2);
A      = zeros(1,2);
pi     = zeros(1,2);
pb_m   = zeros(1,2);
pe_m   = zeros(1,2);
c_m    = zeros(1,2);

for i = 1:2
x0_t = [b_L(i) b_L(i) 0.4 pb pe c_ 1.1885 0.5454 ];  
N = 100000 ;
options = optimset('MaxFunEvals', N, 'MaxIter', N, 'Algorithm', 'levenberg-marquardt');
[x,fval_cali]=fsolve(@(x)myfun_2a2c(x,alpha_L,alpha_H,am(i),gamma,k0(i),k1(i),k2(i),wf_L(i),wf_H(i),b_L(i),b_H(i),N_L(i),N_H(i),rho,tf_L(i),tf_H(i)),x0_t,options);
    bd_L_m(i)    = x(1) ;
    be_L_m(i)    = x(2) ;
    R_m(i)       = x(3) ;
    b_L_m(i)     = (1-R_m(i))*bd_L_m(i)+R_m(i)*be_L_m(i) ;
    b_total_m(i) = N_L(i)*(b_L_m(i))+N_H(i)*b_H(i) ;
    pb_m(i)      = x(4) ;
    pe_m(i)      = x(5) ;
    c_m(i)       = x(6) ;
    A(i)         = x(7) ; 
    pi(i)        = x(8) ; 
    tfd_L_m(i)   = (1-k0(i)-((gamma*pb_m(i)*bd_L_m(i)+(gamma*k2(i)+alpha_L)*wf_L(i)*k1(i)*bd_L_m(i)^k2(i))/(alpha_L*wf_L(i)))) ;
    tfe_L_m(i)   = (1-k0(i)-((gamma*(pb_m(i)+pe_m(i))*be_L_m(i)+(gamma*k2(i)+alpha_L)*wf_L(i)*k1(i)*be_L_m(i)^k2(i))/(alpha_L*wf_L(i)))) ;
    tf_L_m(i)    = (1-R_m(i))*tfd_L_m(i)+R_m(i)*tfe_L_m(i) ;
    tf_H_m(i)    = (1-k0(i)-((gamma*(pb_m(i)+pe_m(i))*b_H(i)+(gamma*k2(i)+alpha_H)*wf_H(i)*k1(i)*b_H(i)^k2(i))/(alpha_H*wf_H(i)))) ;
    Z_L_m(i)     = N_L(i)*tf_L_m(i) ;
    Z_H_m(i)     = N_H(i)*tf_H_m(i) ;
    beta(i)      = (-((1+gamma)*log((pb_m(i)*bd_L_m(i)+wf_L(i)*k1(i)*k2(i)*(bd_L_m(i))^(k2(i)))/((pb_m(i)+pe_m(i))*be_L_m(i)+wf_L(i)*k1(i)*k2(i)*(be_L_m(i))^(k2(i))))+alpha_L*log(bd_L_m(i)/be_L_m(i)))/(log(wf_L(i)/wf_H(i)))) ;
    calibration_P2    = [beta(i) A(i) pi(i) R_m(i)];
    disp('      beta          A             pi          R_m')
    disp(calibration_P2)
end
%%disp('        R_m           R_m')
%%fprintf('    %6.6f   %6.6f\n\n',R_m)
%%disp('        beta           beta')
%%fprintf('    %6.6f   %6.6f\n\n',beta)
%disp('        pb_m         pb_m')
%fprintf('    %6.6f   %6.6f\n\n',pb_m)
%disp('       pe_m         pe_m')
%fprintf('    %6.6f   %6.6f\n\n',pe_m)
%disp('         c_m           c_m')
%fprintf('    %6.6f   %6.6f\n\n',c_m)
disp('         b_L            b_L')
fprintf('    %6.6f   %6.6f\n\n',b_L)
disp('         b_H            b_H')
fprintf('    %6.6f   %6.6f\n\n',b_H)
disp('         pb')
fprintf('    %6.6f\n\n',pb)
disp('         pe')
fprintf('    %6.6f\n\n',pe)
disp('          c_')
fprintf('    %6.6f\n\n',c_)
disp('         b_t            b_t')
fprintf('    %6.6f   %6.6f\n\n\n\n\n',b_total)
disp('         tf_L            tf_L')
fprintf('    %6.6f   %6.6f\n\n\n\n\n',tf_L)
disp('          tf_H          tf_H ')
fprintf('    %6.6f   %6.6f\n\n\n\n\n',tf_H)
disp('          wf_L           wf_L')
fprintf('    %6.6f   %6.6f\n\n',wf_L)
disp('       wf_H          wf_H')
fprintf('    %6.6f   %6.6f\n\n',wf_H)
disp('         am             am')
fprintf('    %6.6f   %6.6f\n\n',am)
disp('          k0               k0')
fprintf('    %6.6f   %6.6f\n\n',k0)
disp('          k1               k1')
fprintf('    %6.6f   %6.6f\n\n',k1)
disp('          k2               k2')
fprintf('    %6.6f   %6.6f\n\n',k2)



%caliabration_check
bd_L_m      = zeros(1,2);
be_L_m      = zeros(1,2);
b_L_m       = zeros(1,2);
b_H_m       = zeros(1,2);
b_total_m   = zeros(1,2);
wf_L_m      = zeros(1,2);
wf_H_m      = zeros(1,2);
R_m         = zeros(1,2);
tfd_L_m     = zeros(1,2);
tfe_L_m     = zeros(1,2);
tf_L_m      = zeros(1,2);
tf_H_m      = zeros(1,2);
lfd_L_m     = zeros(1,2);
lfe_L_m     = zeros(1,2);
lf_L_m      = zeros(1,2);
lf_H_m      = zeros(1,2);
cd_L_m      = zeros(1,2);
ce_L_m      = zeros(1,2);
c_L_m       = zeros(1,2);
c_H_m       = zeros(1,2);
Z_L_m       = zeros(1,2);
Z_H_m       = zeros(1,2);
am_m        = zeros(1,2);
u_L         = zeros(1,2);
u_H         = zeros(1,2);
am_LL       = zeros(1,2);
am_LH       = zeros(1,2);
am_HH       = zeros(1,2);
%pb_m   =   pb_m*1.01;
%pe_m   =   pe_m*1.01;
%k0 = k0*0.99 ;
k1 = k1*0.99 ;
%am      =     am*1.01; 
%c_m     =    c_m*1.01; 
%A =  A*1.01;
%Q = 0.01 ;
for i = 1:2
    a = i ; 
    b = i ;
    c = i ;
    d = i ;
    e = i ;
xi = [b_L(i), b_L(i), b_H(i), R(i), wf_L(i), wf_H(i)] ;
N = 1000000 ;
options = optimset('MaxFunEvals',N,'MaxIter',N, 'Algorithm', 'levenberg-marquardt');
% j == >  preference parameter 選取變數(for phi eta gamma)
[x,fval_m]=fsolve(@(x)myfun_2a_c(x,pb_m(a),pe_m(b),am(c),c_m(d),k0(e),k1(e),k2(e),N_L(i),N_H(i),rho,alpha_L,alpha_H,gamma,beta(i),A(i),pi(i)),xi,options);
    %pb_m(i)      = pb_m(a)*(1-Q) ;
    %am_m(i)      = (am(c)-(N_L(i)*((1-x(4))*x(1)+x(4)*x(2))+N_H(i)*x(3))*pb_m(i)*Q) ; %p1
    %am_m(i)      = (am(c)-(N_L(i)*((1-x(4))*x(1)+x(4)*x(2)))*pb_m(i)*Q) ; %p2
    %am_m(i)      = (am(c)-(N_L(i)*((1-x(4))*x(1)+x(4)*x(2)))*pb_m(i)*Q/N_H(i)) ; %p3
    pe_m(i)      = pe_m(b) ;
    c_m(i)       = c_m(d) ;
    k0(i)        = k0(e) ;
    k1(i)        = k1(e) ;
    k2(i)        = k2(e) ;
    bd_L_m(i)    = x(1) ;
    be_L_m(i)    = x(2) ;
    b_L_m(i)     = ((1-x(4))*x(1)+x(4)*x(2)) ;
    b_H_m(i)     = x(3) ;
    R_m(i)       = x(4) ;
    b_total_m(i) = (N_L(i)*((1-x(4))*x(1)+x(4)*x(2))+N_H(i)*x(3)) ;
    wf_L_m(i)    = x(5) ;
    wf_H_m(i)    = x(6) ;
    tfd_L_m(i)   = (1-k0(i)-((gamma*pb_m(i)*bd_L_m(i)+(gamma*k2(i)+alpha_L)*wf_L(i)*k1(i)*bd_L_m(i)^k2(i))/(alpha_L*wf_L(i)))) ;
    tfe_L_m(i)   = (1-k0(i)-((gamma*(pb_m(i)+pe_m(i))*be_L_m(i)+(gamma*k2(i)+alpha_L)*wf_L(i)*k1(i)*be_L_m(i)^k2(i))/(alpha_L*wf_L(i)))) ;
    tf_L_m(i)    = (1-R_m(i))*tfd_L_m(i)+R_m(i)*tfe_L_m(i) ;
    tf_H_m(i)    = (1-k0(i)-((gamma*(pb_m(i)+pe_m(i))*b_H(i)+(gamma*k2(i)+alpha_H)*wf_H(i)*k1(i)*b_H(i)^k2(i))/(alpha_H*wf_H(i)))) ;
    lfd_L_m(i)   = 1-tfd_L_m(i)-k0(i)-k1(i)*bd_L_m(i)^k2(i) ;
    lfe_L_m(i)   = 1-tfe_L_m(i)-k0(i)-k1(i)*be_L_m(i)^k2(i) ;
    lf_L_m(i)    = (1-R_m(i))*lfd_L_m(i)+R_m(i)*lfe_L_m(i) ;
    lf_H_m(i)    = 1-tf_H_m(i)-k0(i)-k1(i)*b_H_m(i)^k2(i) ;
    cd_L_m(i)    = (pb_m(i)*bd_L_m(i)+wf_L_m(i)*k1(i)*k2(i)*bd_L_m(i)^k2(i))/(alpha_L*wf_L_m(i))+c_m(i) ;
    ce_L_m(i)    = ((pb_m(i)+pe_m(i))*be_L_m(i)+wf_L_m(i)*k1(i)*k2(i)*be_L_m(i)^k2(i))/(alpha_L*wf_L_m(i))+c_m(i) ;
    c_L_m(i)     = (1-R_m(i))*cd_L_m(i)+R_m(i)*ce_L_m(i) ;
    c_H_m(i)     = ((pb_m(i)+pe_m(i))*b_H_m(i)+wf_H_m(i)*k1(i)*k2(i)*b_H_m(i)^k2(i))/(alpha_H*wf_H_m(i))+c_m(i) ;
    Z_L_m(i)     = N_L(i)*tf_L_m(i) ;
    Z_H_m(i)     = N_H(i)*tf_H_m(i) ;
    u_L(i)       = log(cd_L_m(i)-c_m(i))+alpha_L*log(bd_L_m(i))+beta(i)*log(wf_L_m(i))+gamma*log(lfd_L_m(i)) ; %
    u_H(i)       = log(c_H_m(i)-c_m(i))+alpha_H*log(b_H_m(i))+beta(i)*log(wf_H_m(i))+gamma*log(lf_H_m(i)); %
    %am_LL(i)     = am(c)-(N_L(i)*((1-x(4))*x(1)+x(4)*x(2))+N_H(i)*x(3))*pb_m(i)*Q+x(1)*pb_m(i)*Q ; %p1_2
    %am_LH(i)     = am(c)-(N_L(i)*((1-x(4))*x(1)+x(4)*x(2))+N_H(i)*x(3))*pb_m(i)*Q+x(2)*pb_m(i)*Q ; %p1_2
    %am_HH(i)     = am(c)-(N_L(i)*((1-x(4))*x(1)+x(4)*x(2))+N_H(i)*x(3))*pb_m(i)*Q+x(3)*pb_m(i)*Q ; %p1_2
    %am_LL(i)     = am(c)-(N_L(i)*((1-x(4))*x(1)+x(4)*x(2)))*pb_m(i)*Q+x(1)*pb_m(i)*Q ; %p2_2
    %am_LH(i)     = am(c)-(N_L(i)*((1-x(4))*x(1)+x(4)*x(2)))*pb_m(i)*Q+x(2)*pb_m(i)*Q ; %p2_2
    %am_HH(i)     = am(c)-(N_L(i)*((1-x(4))*x(1)+x(4)*x(2)))*pb_m(i)*Q ; %p2_2
    %am_LL(i)     = am(c)+x(1)*pb_m(i)*Q ; %p3_2
    %am_LH(i)     = am(c)+x(2)*pb_m(i)*Q ; %p3_2
    %am_HH(i)     = am(c)-(N_L(i)*((1-x(4))*x(1)+x(4)*x(2)))*pb_m(i)*Q/N_H(i) ; %p3_2
end
disp('            R_m           R_m')
fprintf('    %6.6f   %6.6f\n\n',R_m)
disp('        bd_L            bd_L')
fprintf('    %6.6f   %6.6f\n\n',bd_L_m)
disp('        be_L           be_L')
fprintf('    %6.6f   %6.6f\n\n',be_L_m)
disp('       b_L_m        b_L_m')
fprintf('    %6.6f   %6.6f\n\n',b_L_m)
disp('      b_H_m      b_H_m')
fprintf('    %6.6f   %6.6f\n\n',b_H_m)
disp('       pb_m         pb_m')
fprintf('    %6.6f   %6.6f\n\n',pb_m)
disp('       pe_m         pe_m')
fprintf('    %6.6f   %6.6f\n\n',pe_m)
disp('         c_m           c_m')
fprintf('    %6.6f   %6.6f\n\n',c_m)
disp('       b_t_m        b_t_m')
fprintf('    %6.6f   %6.6f\n\n',b_total_m)
disp('       lfd_L            lfd_L')
fprintf('    %6.6f   %6.6f\n\n',lfd_L_m)
disp('       tfd_L            tfd_L')
fprintf('    %6.6f   %6.6f\n\n',tfd_L_m)
disp('       cd_L            cd_L')
fprintf('    %6.6f   %6.6f\n\n',cd_L_m)
disp('        lfe_L           lfe_L')
fprintf('    %6.6f   %6.6f\n\n',lfe_L_m)
disp('       tfe_L            tfe_L')
fprintf('    %6.6f   %6.6f\n\n',tfe_L_m)
disp('        ce_L            ce_L')
fprintf('    %6.6f   %6.6f\n\n',ce_L_m)
disp('      lf_L_m        lf_L_m')
fprintf('    %6.6f   %6.6f\n\n',lf_L_m)
disp('     tf_L_m        tf_L_m')
fprintf('    %6.6f   %6.6f\n\n',tf_L_m)
disp('       c_L_m        c_L_m')
fprintf('    %6.6f   %6.6f\n\n',c_L_m)
disp('      lf_H_m       lf_H_m')
fprintf('    %6.6f   %6.6f\n\n',lf_H_m)
disp('      tf_H_m      tf_H_m')
fprintf('    %6.6f   %6.6f\n\n',tf_H_m)
disp('      c_H_m        c_H_m')
fprintf('    %6.6f   %6.6f\n\n',c_H_m)
disp('      wf_L_m     wf_L_m')
fprintf('    %6.6f   %6.6f\n\n',wf_L_m)
disp('     wf_H_m     wf_H_m')
fprintf('    %6.6f   %6.6f\n\n',wf_H_m)
disp('         u_L            u_L')
fprintf('    %6.6f   %6.6f\n\n',u_L)
disp('        u_H          u_H')
fprintf('    %6.6f   %6.6f\n\n',u_H)
disp('      am_m        am_m')
fprintf('    %6.6f   %6.6f\n\n',am_m)
disp('      am_LL        am_LL')
fprintf('    %6.6f   %6.6f\n\n',am_LL)
disp('      am_LH       am_LH')
fprintf('    %6.6f   %6.6f\n\n',am_LH)
disp('      am_HH      am_HH')
fprintf('    %6.6f   %6.6f\n\n',am_HH)
%%
V_L         = zeros(1,2);
V_H         = zeros(1,2);
for i = 1:2
  %u_L_n = [-0.799514   -0.679298] ;  %rho =  0.401 p1
  %u_H_n = [-1.253773   -1.074555] ;  %rho =  0.401 p1
  %u_L_n = [-0.799065   -0.678238] ;  %rho =  0.401 p2
  %u_H_n = [-1.259221   -1.080921] ;  %rho =  0.401 p2
  %u_L_n = [-0.793854   -0.674262] ;  %rho =  0.401 p3
  %u_H_n = [-1.275048   -1.086167] ;  %rho =  0.401 p3
  %u_L_n = [-0.811617   -0.692652] ;  %rho =  0.401 p1_2
  %u_H_n = [-1.257854   -1.079507] ;  %rho =  0.401 p1_2
  %u_L_n = [-0.811168   -0.691590] ;  %rho =  0.401 p2_2
  %u_H_n = [-1.259209   -1.080906] ;  %rho =  0.401 p2_2
  %u_L_n = [-0.806001   -0.687661] ;  %rho =  0.401 p3_2
  %u_H_n = [-1.274918   -1.086094] ;  %rho =  0.401 p3_2
  %u_L_n = [-0.799556   -0.679373] ;  %rho = -2.030 p1
  %u_H_n = [-1.253790   -1.074589] ;  %rho = -2.030 p1
  %u_L_n = [-0.799107   -0.678313] ;  %rho = -2.030 p2
  %u_H_n = [-1.259239   -1.080955] ;  %rho = -2.030 p2
  %u_L_n = [-0.793855   -0.674265] ;  %rho = -2.030 p3
  %u_H_n = [-1.275202   -1.086301] ;  %rho = -2.030 p3
  %u_L_n = [-0.811659   -0.692725] ;  %rho = -2.030 p1_2
  %u_H_n = [-1.257871   -1.079540] ;  %rho = -2.030 p1_2
  %u_L_n = [-0.811210   -0.691663] ;  %rho = -2.030 p2_2
  %u_H_n = [-1.259226   -1.080939] ;  %rho = -2.030 p2_2
  u_L_n = [-0.806002   -0.687664] ;  %rho = -2.030 p3_2
  u_H_n = [-1.275069   -1.086225] ;  %rho = -2.030 p3_2
  f_L=@(VL)(log(cd_L_m(i)*(1+VL)-c_m(i))+alpha_L*log(bd_L_m(i))+beta(i)*log(wf_L_m(i))+gamma*log(lfd_L_m(i))-u_L_n(i));
  [VL, fval_VL]=fsolve(f_L,1);
  V_L(i) = VL*100 ;
  f_H=@(VH)(log(c_H_m(i)*(1+VH)-c_m(i))+alpha_H*log(b_H_m(i))+beta(i)*log(wf_H_m(i))+gamma*log(lf_H_m(i))-u_H_n(i));
  [VH, fval_VH]=fsolve(f_H,1);
  V_H(i) = VH*100 ;
end
disp('        V_L(%)       V_L(%)')
fprintf('    %6.6f   %6.6f\n\n',V_L)
disp('        V_H(%)      V_H(%)')
fprintf('    %6.6f   %6.6f\n\n',V_H)
