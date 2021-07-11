%Master Thesis R06323052 朝霍 

%calibration preference parameter

%calculate

%price parameter
pb_t = [0.0982 0.0995 0.1256 0.1141] ;  %from calibration 基闽跑计
pe_t = [0.0760 0.0611 0.0509 0.0326] ;  %from calibration 基闽跑计
R    = [0.411  0.402  0.401  0.436];    %from PSFD (L->L percentage)
pn_c = R.*pb_t+(1-R).*(pb_t+pe_t);      %for calibration

am   = ones(1,4);                       %from calibration 基闽跑计 
wf_t = [0.5406 0.5915 0.6223 0.6442] ;  %from calibration 基闽跑计 
wf_L = [0.4499 0.4749 0.4870 0.4885] ;  %from calibration 基闽跑计
wf_H = [0.9290 0.9165 0.8721 0.8223] ;  %from calibration 基闽跑计

c_   = [0.6359 0.6713 0.6632 0.7228] ;  %from calibration 基闽跑计

%fertitlity
n_t  = [2.309 2.154 1.923 1.698];       %from キАネ▅羆计(64)

%time related
k0_t = [0.2021 0.1875 0.1740 0.1615];   %from home_production.m
k1_t = [0.1344 0.1385 0.1333 0.1323];   %from home_production.m
k2_t = [0.5434 0.4839 0.4263 0.4207];   %from home_production.m
x_t  = k0_t+ k1_t .*n_t .^ k2_t;
tf_t = [144 168 187 201];               %from Human
tf_t = tf_t/60/16;              
lf_t = ones(1,4)-x_t-tf_t ; 

%cliabration
r    = (wf_t.*lf_t)./(wf_t.*tf_t+am-c_-pn_c.*n_t);
a    = (pn_c.*n_t+wf_t.*k1_t.*k2_t.*n_t.^k2_t)./(wf_t.*lf_t).*r;



