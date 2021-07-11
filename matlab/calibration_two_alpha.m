%Master Thesis R06323052 陳思齊 
%calibration (wage deifinition is only wage(excluding owner))
clear
clc

%price parameter

%1 definition(受雇人員報酬)
am_1      = [1.0000 1.1527 1.3158];  %from calibration 價格相關變數(受雇人員報酬)
wf_L_1    = [0.5182 0.5990 0.6938];  %from calibration 價格相關變數(受雇人員報酬)
wf_H_1    = [1.0034 1.1081 1.2063];  %from calibration 價格相關變數(受雇人員報酬)
am_1n     = [1.0000 1.0000*1.0425 1.0000*1.0425*1.0718];  %from calibration wage_growing(受雇人員報酬)
wf_L_1n   = [0.5182 0.5182*1.0327 0.5182*1.0327*1.0605];  %from calibration wage_growing(受雇人員報酬)(薪資以學歷區分各別算成長率)
wf_H_1n   = [1.0034 1.0034*1.0381 1.0034*1.0381*1.0687];  %from calibration wage_growing(受雇人員報酬)(薪資以學歷區分各別算成長率)
wf_L_1tn  = [0.5182 0.5182*1.0155 0.5182*1.0155*1.0278];  %from calibration wage_growing(受雇人員報酬)(薪資以學歷區分各別算成長率)
wf_H_1tn  = [1.0034 1.0034*1.0155 1.0034*1.0155*1.0278];  %from calibration wage_growing(受雇人員報酬)(以平均薪資計算成長)

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
am_2      = [1.0000 1.1461 1.3091];  %from calibration 價格相關變數(受雇人員報酬+產業主所得)
wf_L_2    = [0.4758 0.5625 0.6490];  %from calibration 價格相關變數(受雇人員報酬+產業主所得)
wf_H_2    = [0.9429 1.0197 1.1240];  %from calibration 價格相關變數(受雇人員報酬+產業主所得)
am_2n     = [1.0000 1.0000*1.0425 1.000*1.0425*1.0718];  %from calibration wage_growing(受雇人員報酬+產業主所得)
wf_L_2n   = [0.4758 0.4758*1.0327 0.4758*1.0327*1.0605];  %from calibration wage_growing(受雇人員報酬+產業主所得)(薪資以學歷區分各別算成長率)
wf_H_2n   = [0.9429 1.0034*1.0381 1.0034*1.0381*1.0687];  %from calibration wage_growing(受雇人員報酬+產業主所得)(薪資以學歷區分各別算成長率)
wf_L_2tn  = [0.4758 0.4758*1.0155 0.4758*1.0155*1.0278];  %from calibration wage_growing(受雇人員報酬+產業主所得)(薪資以學歷區分各別算成長率)
wf_H_2tn  = [0.9429 0.9429*1.0155 0.9429*1.0155*1.0278];  %from calibration wage_growing(受雇人員報酬+產業主所得)(以平均薪資計算成長)
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

%3 definition(經常性收入)
am_3      = [1.0000 1.1511 1.3137];  %from calibration 價格相關變數(經常性收入)
wf_L_3    = [0.4483 0.5351 0.6212];  %from calibration 價格相關變數(經常性收入)
wf_H_3    = [0.8748 0.9520 1.0560];  %from calibration 價格相關變數(經常性收入)
am_3n     = [1.0000 1.0000*1.0425 1.000*1.0425*1.0718];   %from calibration wage_growing(經常性收入)
wf_L_3n   = [0.4483 0.4483*1.0327 0.4483*1.0327*1.0605];  %from calibration wage_growing(經常性收入)(薪資以學歷區分各別算成長率)
wf_H_3n   = [0.8748 0.8748*1.0381 0.8748*1.0381*1.0687];  %from calibration wage_growing(經常性收入)(薪資以學歷區分各別算成長率)
wf_L_3tn  = [0.4483 0.4483*1.0155 0.4483*1.0155*1.0278];  %from calibration wage_growing(經常性收入)(薪資以學歷區分各別算成長率)
wf_H_3tn  = [0.8748 0.8748*1.0155 0.8748*1.0155*1.0278];  %from calibration wage_growing(經常性收入)(以平均薪資計算成長)
%original(including 1822) 
pb_3_o  = [0.0710 0.0891 0.1260]; 
pe_3_o  = [0.0520 0.0556 0.0555];
c_3_ao  = [0.5016 0.5795 0.6476]; % a: average(分層加權) ((生活開銷-子女部分)/成年人*2)
c_3_to  = [0.5575 0.6362 0.6953]; % t: total(全部平均)   ((生活開銷-子女部分)/成年人*2)
c_3_aom = [0.6064 0.6829 0.7492]; % a: average(分層加權) (僅考慮子女(<23歲)+夫妻家庭)
c_3_tom = [0.6626 0.7363 0.7873]; % t: total(全部平均)   (僅考慮子女(<23歲)+夫妻家庭)
%new(excluding 1822)
pb_3_n  = [0.0602 0.1082 0.1157];
pe_3_an = [0.0720 0.0801 0.0875]; % a: average(分層加權)
c_3_an  = [0.4998 0.5687 0.6390]; % a: average(分層加權)
pe_3_tn = [0.0767 0.0849 0.0918]; % t: total(全部平均)
c_3_tn  = [0.5512 0.6221 0.6861]; % t: total(全部平均)


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

for i = 1
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
calibration_P1    = [alpha_L alpha_H gamma];
disp('      alp_L      alp_H   gamma')
disp(calibration_P1)
calibration_P2    = [beta(i) A(i) pi(i) R_m(i)];
disp('      beta          A             pi          R_m')
disp(calibration_P2)
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

for i = 1:3
x0_t = [b_L(i) b_L(i) 0.4 pb(i) pe(i) c_(i) 1.1885 0.5454 ];  
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
    calibration_P2    = [beta(i) A(i) pi(i) R_m(i)];
    disp('      beta          A             pi          R_m')
    disp(calibration_P2)
end
disp('     pb_m     pb_m      pb_m')
disp(pb_m)
disp('     pe_m     pe_m      pe_m')
disp(pe_m)
disp('      c_m        c_m         c_m')
disp(c_m)


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
lfd_L_m     = zeros(1,3);
lfe_L_m     = zeros(1,3);
lf_L_m      = zeros(1,3);
lf_H_m      = zeros(1,3);
cd_L_m      = zeros(1,3);
ce_L_m      = zeros(1,3);
c_L_m       = zeros(1,3);
c_H_m       = zeros(1,3);
Z_L_m       = zeros(1,3);
Z_H_m       = zeros(1,3);

%pb_m = 1.01*pb_m  ;
%k1 = k1*0.985 ;
   Q = 0.0 ;
for i = 1:3
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
    pb_m(i)      = pb_m(a)*(1-Q) ;
    pe_m(i)      = pe_m(b) ;
    am(i)        = am(c)-(N_L(i)*(x(4)*x(1)+(1-x(4))*x(2)))*pb_m(i)*Q ;
    c_m(i)       = c_m(d) ;
    k0(i)        = k0(e) ;
    k1(i)        = k1(e) ;
    k2(i)        = k2(e) ;
    bd_L_m(i)    = x(1) ;
    be_L_m(i)    = x(2) ;
    b_L_m(i)     = (x(4)*x(1)+(1-x(4))*x(2)) ;
    b_H_m(i)     = x(3) ;
    R_m(i)       = x(4) ;
    b_total_m(i) = (N_L(i)*(x(4)*x(1)+(1-x(4))*x(2))+N_H(i)*x(3)) ;
    wf_L_m(i)    = x(5) ;
    wf_H_m(i)    = x(6) ;
    tfd_L_m(i)   = (1-k0(i)-((gamma*pb_m(i)*bd_L_m(i)+(gamma*k2(i)+alpha_L)*wf_L(i)*k1(i)*bd_L_m(i)^k2(i))/(alpha_L*wf_L(i)))) ;
    tfe_L_m(i)   = (1-k0(i)-((gamma*(pb_m(i)+pe_m(i))*be_L_m(i)+(gamma*k2(i)+alpha_L)*wf_L(i)*k1(i)*be_L_m(i)^k2(i))/(alpha_L*wf_L(i)))) ;
    tf_L_m(i)    = R_m(i)*tfd_L_m(i)+(1-R_m(i))*tfe_L_m(i) ;
    tf_H_m(i)    = (1-k0(i)-((gamma*(pb_m(i)+pe_m(i))*b_H(i)+(gamma*k2(i)+alpha_H)*wf_H(i)*k1(i)*b_H(i)^k2(i))/(alpha_H*wf_H(i)))) ;
    lfd_L_m(i)   = 1-tfd_L_m(i)-k0(i)-k1(i)*bd_L_m(i)^k2(i) ;
    lfe_L_m(i)   = 1-tfe_L_m(i)-k0(i)-k1(i)*be_L_m(i)^k2(i) ;
    lf_L_m(i)    = R_m(i)*lfd_L_m(i)+(1-R_m(i))*lfe_L_m(i) ;
    lf_H_m(i)    = 1-tf_H_m(i)-k0(i)-k1(i)*b_H_m(i)^k2(i) ;
    cd_L_m(i)    = (pb_m(i)*bd_L_m(i)+wf_L_m(i)*k1(i)*k2(i)*bd_L_m(i)^k2(i))/(alpha_L*wf_L_m(i)) ;
    ce_L_m(i)    = ((pb_m(i)+pe_m(i))*be_L_m(i)+wf_L_m(i)*k1(i)*k2(i)*be_L_m(i)^k2(i))/(alpha_L*wf_L_m(i)) ;
    c_L_m(i)     = R_m(i)*cd_L_m(i)+(1-R_m(i))*ce_L_m(i) ;
    c_H_m(i)    = ((pb_m(i)+pe_m(i))*b_H_m(i)+wf_H_m(i)*k1(i)*k2(i)*b_H_m(i)^k2(i))/(alpha_H*wf_H_m(i)) ;
    Z_L_m(i)     = N_L(i)*tf_L_m(i) ;
    Z_H_m(i)     = N_H(i)*tf_H_m(i) ;
end
disp('      R_m        R_m       R_m')
disp(R_m)
disp('       b_L          b_L         b_L')
disp(b_L)
%disp('   bd_L_m bd_L_m  bd_L_m')
%disp(bd_L_m)
%disp('   be_L_m  be_L_m  be_L_m')
%disp(be_L_m)
disp('    b_L_m    b_L_m    b_L_m')
disp(b_L_m)
disp('       b_H         b_H        b_H')
disp(b_H)
disp('   b_H_m    b_H_m   b_H_m')
disp(b_H_m)
%disp('       pb            pb           pb')
%disp(pb)
disp('     pb_m      pb_m     pb_m')
disp(pb_m)
%disp('        pe           pe           pe')
%disp(pe)
disp('     pe_m      pe_m     pe_m')
disp(pe_m)
%disp('        c_            c_            c_')
%disp(c_)
%%
disp('      c_m         c_m        c_m')
disp(c_m)
disp('    b_t_m     b_t_m     b_t_m')
disp(b_total_m)
disp('      lfd_L        lfd_L      lfd_L')
disp(lfd_L_m)
disp('     tfd_L        tfd_L       tfd_L')
disp(tfd_L_m)
disp('      cd_L        cd_L      cd_L')
disp(cd_L_m)
disp('      lfe_L        lfe_L      lfe_L')
disp(lfe_L_m)
disp('     tfe_L        tfe_L       tfe_L')
disp(tfe_L_m)
%disp('       tf_L         tf_L          tf_L')
disp('      ce_L        ce_L        ce_L')
disp(ce_L_m)
%disp(tf_L)
disp('    lf_L_m    lf_L_m     lf_L_m')
disp(lf_L_m)
disp('    tf_L_m    tf_L_m    tf_L_m')
disp(tf_L_m)
disp('    c_L_m     c_L_m    c_L_m')
disp(c_L_m)
%disp('       tf_H         tf_H          tf_H')
%disp(tf_H)
disp('   lf_H_m    lf_H_m   lf_H_m')
disp(lf_H_m)
disp('   tf_H_m    tf_H_m   tf_H_m')
disp(tf_H_m)
disp('   c_H_m    c_H_m    c_H_m')
disp(c_H_m)
disp('   wf_L_m  wf_L_m   wf_L_m')
disp(wf_L_m)
disp('   wf_H_m wf_H_m wf_H_m')
disp(wf_H_m)
disp('    am_m     am_m     am_m')
disp(am)