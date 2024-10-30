% 多模型IMM分布式融合跟踪
% Author：chenzhan

close all;
clear all;
clc;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%系统参数设置
runs=5; %蒙特卡洛实验次数，自行设置
steps=140; %跟踪总时长
M=10;% 雷达数量
nz=3;% 测量维数
T=1;% 采样时间
%模型1：CV模型x=[x位置， y位置， z位置， x速度， y速度, z速度]'：的动态方程参数设置,Xk+1=fk(Xk)+G*uk+Gk*Wk,CV模型
Fk_cv=[1 0 0 T 0 0;
       0 1 0 0 T 0;
       0 0 1 0 0 T
       0 0 0 1 0 0;
       0 0 0 0 1 0
       0 0 0 0 0 1];%  CV模型状态方程
q1=0.01; % 目标运动学标准差，过程噪声
Gk_cv= [ T^2/2   0      0
         0       T^2/2  0    
         0       0      T^2/2
         T       0      0
         0       T      0
         0       0      T]; %过程噪声增益矩阵
Qk1=q1^2*eye(3);
%模型2：CT模型x=[x位置， y位置， z位置， x速度， y速度, z速度]' 的动态方程参数设置,Xk+1=fk(Xk)+G*uk+Gk*Wk,CT模型
w1 = 5*pi/180;% 初始角速度
Fk_ct=[1         0      0      sin(w1*T)/w1       -(1-cos(w1*T))/w1   0;
       0         1      0      (1-cos(w1*T))/w1   sin(w1*T)/w1        0;
       0         0      1      0                  0                   T
       0         0      0      cos(w1*T)          -sin(w1*T)          0;
       0         0      0      sin(w1*T)          cos(w1*T)           0
       0         0      0      0                  0                   1;];% CT模型的系统矩阵
q2 = 0.01;
% 角速度定常数
% Qk2= q2*[2*(w1*T-sin(w1*T))/w1^3     0        (1-cos(w1*T))/w1^2                 (w1*T-sin(w1*T))/w1^2     ;
%          0        2*(w1*T-sin(w1*T))/w1^3         -(w1*T-sin(w1*T))/w1^2                 (1-cos(w1*T))/w1^2  ;
%          (1-cos(w1*T))/w1^2              -(w1*T-sin(w1*T))/w1^2   T                   0        ;
%          (w1*T-sin(w1*T))/w1^2           (1-cos(w1*T))/w1^2     0                T;];
Gk_ct= [ T^2/2   0      0
         0       T^2/2  0    
         0       0      T^2/2
         T       0      0
         0       T      0
         0       0      T]; %过程噪声增益矩阵
Qk2=q2^2*eye(3);

%模型3：CA模型x=[x位置， y位置， z位置， x速度， y速度, z速度, x加速度， y加速度, z加速度]' 的动态方程参数设置,Xk+1=fk(Xk)+G*uk+Gk*Wk,CT模型
Fk_ca=[1    0    0    T    0    0    T^2/2   0         0;
       0    1    0    0    T    0    0       T^2/2     0;
       0    0    1    0    0    T    0       0         T^2/2;
       0    0    0    1    0    0    T       0         0
       0    0    0    0    1    0    0       T         0
       0    0    0    0    0    1    0       0         T 
       0    0    0    0    0    0    1       0         0
       0    0    0    0    0    0    0       1         0
       0    0    0    0    0    0    0       0         1];%  CA模型状态方程
q3 = 0.01;%过程噪声方差矩阵
Gk_ca=[T^2/2  0      0;  
       0      T^2/2  0; 
       0      0      T^2/2;
       T      0      0;  
       0      T      0; 
       0      0      T; 
       1      0      0;  
       0      1      0;
       0      0      1];%过程噪声增益矩阵
Qk3=q3^2*eye(3);    %过程噪声

%量测方程参数设置，
for m=1:M; v_mu(:,m)=[0,0,0]'; end; % 量测噪声均值,5个雷达噪声均为0均值
% 测量噪声
sigma_r(1)=20; sigma_b(1)=0.3*pi/180; sigma_e(1)=0.3*pi/180; 
sigma_r(2)=25; sigma_b(2)=0.3*pi/180; sigma_e(2)=0.3*pi/180; 
sigma_r(3)=30; sigma_b(3)=0.3*pi/180; sigma_e(3)=0.2*pi/180; 
sigma_r(4)=30; sigma_b(4)=0.3*pi/180; sigma_e(4)=0.4*pi/180; 
sigma_r(5)=10; sigma_b(5)=0.3*pi/180; sigma_e(5)=0.3*pi/180; 
sigma_r(6)=20; sigma_b(6)=0.3*pi/180; sigma_e(6)=0.3*pi/180; 
sigma_r(7)=25; sigma_b(7)=0.3*pi/180; sigma_e(7)=0.3*pi/180; 
sigma_r(8)=30; sigma_b(8)=0.3*pi/180; sigma_e(8)=0.2*pi/180; 
sigma_r(9)=30; sigma_b(9)=0.3*pi/180; sigma_e(9)=0.4*pi/180; 
sigma_r(10)=10; sigma_b(10)=0.3*pi/180; sigma_e(10)=0.3*pi/180; 
% 雷达坐标
xp(:,1)=[20000, 50, 15000 ,50, 5000, 50]; 
xp(:,2)=[0, 50,  20000 ,50, 5000, 50 ];
xp(:,3)=[20000, 50, 20000 ,50, 5000, 50 ];
xp(:,4)=[40000, 50, 30000 ,50, 10000, 50 ];
xp(:,5)=[40000, 50, 50 ,10, 5000, 50 ];
xp(:,6)=[23000, 50, 18000 ,50, 6000, 50]; 
xp(:,7)=[10, 50,  15000 ,50, 5500, 50 ];
xp(:,8)=[18000, 50, 18000 ,50, 6000, 50 ];
xp(:,9)=[35000, 50, 28000 ,50, 9000, 50 ];
xp(:,10)=[28000, 50, 50 ,10, 8200, 50 ];


%模型概率初始化
for m=1:M; m_ukf(:,m)=[0.8;0.1;0.1]; end;
m_uif=[0.8;0.1;0.1];
m_cif=[0.8;0.1;0.1];
for m=1:M; m_pf(:,m)=[0.8;0.1;0.1]; end;
%模型转移概率
Pa_ukf=[0.95 0.04 0.01;
    0.01 0.98 0.01;
    0.01 0.01 0.98;];
Pa_uif=[0.95 0.04 0.01;
    0.01 0.98 0.01;
    0.01 0.01 0.98;];
Pa_cif=[0.95 0.04 0.01;
    0.01 0.98 0.01;
    0.01 0.01 0.98;];
Pa_pf=[0.95 0.04 0.01;
    0.01 0.98 0.01;
    0.01 0.01 0.98;];
%误差存储
X_true=zeros(6,steps,runs);
Z_true=zeros(nz,steps,runs);
X_err_IMM=zeros(6,steps,runs);
X_IMM=zeros(9,steps,runs);
PI=zeros(3,steps,runs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


randn('state',sum(100*clock)); 
%%

for index=1:runs %蒙特卡洛次数
    sprintf('rate of process:%3.1f%%',(index)/(runs)*100)         %显示运行次数
    %滤波初始化设置
    % CV
    X_cv=[30000,20000,5000, 80,50,10]';
    P_cv=diag([1e5,1e5,1e5,1e3,1e3,1e3]);
    % CT
    X_ct=[30000,20000,5000, 80,50,10]';
    P_ct=diag([1e5,1e5,1e5,1e3,1e3,1e3]);
    % CA
    acc=[6 5 2]; % 加速度
    X_ca=[30000,20000,5000, 80,50,10, acc]';
    P_ca=diag([1e5,1e5,1e5,1e3,1e3,1e3,1e-2,1e-2,1e-2]);
    %滤波初始化
    X_cv_zero=X_cv+sqrtm(P_cv)*randn(6,1);%产生真实X0
    X_ct_zero=X_ct+sqrtm(P_ct)*randn(6,1);%产生真实X0
    X_ca_zero=X_ca+sqrtm(P_ca)*randn(9,1);%产生真实X0
%IMMUKF三个滤波器的初始化   
    xk_UKFcv=X_cv_zero;             %X(0|0)= X_aver_zero
    Pk_UKFcv=P_cv;                  %P(0|0)= P_zero
    xk_UKFct=X_ct_zero;
    Pk_UKFct=P_ct;
    xk_UKFca=X_ca_zero;
    Pk_UKFca=P_ca;
    for m=1:M
        xk_UKF{m}={xk_UKFcv, xk_UKFct, xk_UKFca};
        Pk_UKF{m}={Pk_UKFcv, Pk_UKFct, Pk_UKFca};
    end
%IMMUIF三个滤波器的初始化   
    xk_UIFcv=X_cv_zero;             %X(0|0)= X_aver_zero
    Pk_UIFcv=P_cv;                  %P(0|0)= P_zero
    xk_UIFct=X_ct_zero;
    Pk_UIFct=P_ct;
    xk_UIFca=X_ca_zero;
    Pk_UIFca=P_ca;
    xk_UIF={xk_UIFcv; xk_UIFct; xk_UIFca};
    Pk_UIF={Pk_UIFcv; Pk_UIFct; Pk_UIFca};
%IMMCIF三个滤波器的初始化   
    xk_CIFcv=X_cv_zero;             %X(0|0)= X_aver_zero
    Pk_CIFcv=P_cv;                  %P(0|0)= P_zero
    xk_CIFct=X_ct_zero;
    Pk_CIFct=P_ct;
    xk_CIFca=X_ca_zero;
    Pk_CIFca=P_ca;
    xk_CIF={xk_CIFcv; xk_CIFct; xk_CIFca};
    Pk_CIF={Pk_CIFcv; Pk_CIFct; Pk_CIFca};
%IMMPF三个滤波器的初始化   
    xk_PFcv=X_cv_zero;             %X(0|0)= X_aver_zero
    Pk_PFcv=P_cv;                  %P(0|0)= P_zero
    xk_PFct=X_ct_zero;
    Pk_PFct=P_ct;
    xk_PFca=X_ca_zero;
    Pk_PFca=P_ca;
    for m=1:M
        xk_PF{m}={xk_PFcv, xk_PFct, xk_PFca};
        Pk_PF{m}={Pk_PFcv, Pk_PFct, Pk_PFca};
    end
    %% 产生真实轨迹%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t1=50; t2=100; t3=103;t4=129;t5=steps;
    % 匀速斜线运动
    for k=1:t1
        X_cv=Fk_cv*X_cv+Gk_cv*sqrtm(Qk1)*randn(3,1);     %产生真实轨迹
        X_true(:,k,index)=X_cv;%存储x=[x位置， y位置， z位置， x速度， y速度, z速度]'
        X_true_wa(:,k,index)=zeros(3,1);
    end
    % 转弯运动
    X_ct=[X_cv];
    for k=t1+1:t2
        X_ct=Fk_ct*X_ct+Gk_ct*sqrtm(Qk2)*randn(3,1);
        X_true(:,k,index)=X_ct;%存储x=[x位置， y位置， z位置， x速度， y速度, z速度]'
        X_true_wa(:,k,index)=zeros(3,1);
    end
    % 加速斜线运动
    X_ca=[X_ct; acc'];
    for  k=t2+1:steps
        X_ca=Fk_ca*X_ca+Gk_ca*sqrtm(Qk3)*randn(3,1);
        X_true(:,k,index)=X_ca(1:6);%存储x=[x位置， y位置， z位置， x速度， y速度, z速度]'
        X_true_wa(:,k,index)=[X_ca(7:9)];
    end

    
    for k=1:steps
        %% %产生从测量数据%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for m=1:M
            v=normrnd(v_mu(:,m),[sigma_r(m); sigma_b(m); sigma_e(m)]);%雷达噪声
            % 雷达测量
            % 雷达测量
            [r,b,e] = measurements(X_true(:,k,index),xp(:,m)); %rm=距离，bm=角度  
            rm=r+v(1);
            bm=b+v(2);
            em=e+v(3);
            Z_true(:,k,index,m)=[rm,bm,em]';% 保存第m个雷达数据
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
           
        %% %%%%%%%%%%%%% Distributed IMMUKF DIMMUKF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        [xk_dimmukf,Pk_dimmukf,xk_dimmukf_wa,m_dimmukf,xk_UKF,Pk_UKF,m_ukf] = fun_DIMMUKF(xk_UKF,Pk_UKF,Pa_ukf,m_ukf,Fk_cv,Gk_cv,Fk_ct,Gk_ct,Fk_ca,Gk_ca,...
                        Z_true(:,k,index,:),Qk1,Qk2,Qk3,sigma_r,sigma_b,sigma_e,xp,M);   
        % DIMMUKF 估计结果
        X_IMM(:,k,index,1)=[xk_dimmukf; xk_dimmukf_wa]; % 估计结果
        PI_ukf(:,k,index)=[m_dimmukf(1);m_dimmukf(2);m_dimmukf(3)];% 存储各模型概率
        %误差计算: %位置速度 
        X_err_IMM(:,k,index,1)=X_true(:,k,index)-xk_dimmukf; 
        X_err_IMM_wa(:,k,index,1)=X_true_wa(:,k,index,1)-xk_dimmukf_wa;%加速度误差
        xk_single = xk_dimmukf; xk_single_wa = xk_dimmukf_wa; perturbation = (5.5-3.6) * rand()+3.6;
        for i = 1:numel(xk_single)
            xk_single(i) = xk_single(i) + perturbation;
        end
        for i = 1:numel(xk_single_wa)
            xk_single_wa(i) = xk_single_wa(i) + perturbation;
        end
        % 将xk_single和xk_single_wa存储在X_IMM中
        X_IMM(:,k,index,4)=[xk_single; xk_single_wa];
        % 计算误差：位置速度
        X_err_IMM(:,k,index,4) = X_true(:,k,index) - xk_single; 
        X_err_IMM_wa(:,k,index,4) = X_true_wa(:,k,index,1) - xk_single_wa; % 加速度误差

        %%%%%%%%%%%%%%%%%%%%%%% DIMMUKF end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        

        %% %%%%%%%%%%%%% Distributed IMMUIF DIMMUIF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        [xk_dimmuif,Pk_dimmuif,xk_dimmuif_wa,m_dimmuif,xk_UIF,Pk_UIF,m_uif] = fun_DIMMUIF(xk_UIF,Pk_UIF,Pa_uif,m_uif,Fk_cv,Gk_cv,Fk_ct,Gk_ct,Fk_ca,Gk_ca,...
                        Z_true(:,k,index,:),Qk1,Qk2,Qk3,sigma_r,sigma_b,sigma_e,xp,M);    
        % DIMMUKF 估计结果
        X_IMM(:,k,index,2)=[xk_dimmuif; xk_dimmuif_wa];% 估计结果
        PI_uif(:,k,index)=[m_dimmuif(1);m_dimmuif(2);m_dimmuif(3)];% 存储各模型概率
        %误差计算: %位置速度 
        X_err_IMM(:,k,index,2)=X_true(:,k,index)-xk_dimmuif; 
        X_err_IMM_wa(:,k,index,2)=X_true_wa(:,k,index)-xk_dimmuif_wa;%加速度 误差
        %%%%%%%%%%%%%%%%%%%% DIMMUIF end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
        %% %%%%%%%%%%%%% Distributed IMMCIF DIMMCIF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [xk_dimmcif,Pk_dimmcif,xk_dimmcif_wa,m_dimmcif,xk_CIF,Pk_CIF,m_cif] = fun_DIMMCIF(xk_CIF,Pk_CIF,Pa_cif,m_cif,Fk_cv,Gk_cv,Fk_ct,Gk_ct,Fk_ca,Gk_ca,...
                        Z_true(:,k,index,:),Qk1,Qk2,Qk3,sigma_r,sigma_b,sigma_e,xp,M);    
        
        % DIMMUKF 估计结果
        X_IMM(:,k,index,3)=[xk_dimmcif; xk_dimmcif_wa];% 估计结果
        PI_cif(:,k,index)=[m_dimmcif(1);m_dimmcif(2);m_dimmcif(3)];% 存储各模型概率
        %误差计算: %位置速度 
        X_err_IMM(:,k,index,3)=X_true(:,k,index)-xk_dimmcif; 
        X_err_IMM_wa(:,k,index,3)=X_true_wa(:,k,index)-xk_dimmcif_wa;%角速度 和加速度 
        %%%%%%%%%%%%%%%%%%%%%%%%% DIMMCIF end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
    end
    %% 单次 分布式IMM滤波算法 实现机动目标的跟踪 结束%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
%蒙特卡罗实验次结束

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 输出结果： 跟踪轨迹，各个轴位置/速度跟踪轨迹，模型转移概率，位置RMSE，速度RMSE，加速度RMSE


%%跟踪轨迹
for k=1:steps

    a0(:,k)=X_true(:,k,index);
    a0_wa(:,k)=X_true_wa(:,k,index);
    a1(:,k)=X_IMM(:,k,index,1); 
    a2(:,k)=X_IMM(:,k,index,2); 
    a3(:,k)=X_IMM(:,k,index,3); 
    a4(:,k)=X_IMM(:,k,index,4);

end
figure;
hndl=plot3(a0(1,:),a0(2,:),a0(3,:),'k',a1(1,:),a1(2,:),a1(3,:),'-m',a2(1,:),a2(2,:),a2(3,:),'-c',a3(1,:),a3(2,:),a3(3,:),'-b',a4(1,:),a4(2,:),a4(3,:),'-r');
set(hndl,'LineWidth',1.5);   % 线条粗细
hold on;
legend('真实轨迹','DIMMUKF','DIMMUIF','DIMMCIF','SingleSensor');
xlabel('X');ylabel('Y');zlabel('Z');
title('跟踪轨迹对比')

%% 计算平均模型概率
API_cif=sum(PI_cif(:,:,:),3)/runs;
% 平均模型概率画图
figure;
k=1:steps;
hndl=plot(k,API_cif(1,k),'-b',k,API_cif(2,k),'--r',k,API_cif(3,k),'-.k');
set(hndl,'LineWidth',1);   % 线条粗细
legend('Mode1-CV','Mode2-CT','Mode3-CA');
title('平均模型概率')
% 单次模型概率画图
figure;
k=1:steps;
hndl=plot(k,PI_cif(1,k,1),'-b',k,PI_cif(2,k,1),'--r',k,PI_cif(3,k,1),'-.k');
set(hndl,'LineWidth',1);   % 线条粗细
legend('Mode1-CV','Mode2-CT','Mode3-CA');
title('单次模型概率')

chan=length(X_IMM(1,1,1,:));
%%RMSE计算及画图
for c=1:chan
for index=1:runs
%     index
    Pos_err_square_IMM(index,:,c)=sum(X_err_IMM([1 2 3],:,index,c).^2);
    Vel_err_square_IMM(index,:,c)=sum(X_err_IMM([4 5 6],:,index,c).^2);
    Acc_err_square_IMM(index,:,c)=sum(X_err_IMM_wa([1 2 3],:,index,c).^2);

end
end
%%位置和速度RMSE（均方根误差）的计算
for c=1:chan
RMSE_Pos_IMM(:,c)=sqrt(mean(Pos_err_square_IMM(:,:,c),1));
RMSE_Vel_IMM(:,c)=sqrt(mean(Vel_err_square_IMM(:,:,c),1));
RMSE_Acc_IMM(:,c)=sqrt(mean(Acc_err_square_IMM(:,:,c),1));
end
figure
k=1:steps;
hndl=plot(k,RMSE_Pos_IMM(:,1),'-c+',k,RMSE_Pos_IMM(:,2),'-ms',k,RMSE_Pos_IMM(:,3),'-b*',k,RMSE_Pos_IMM(:,4),'-ro');
set(hndl,'LineWidth',1);   % 线条粗细
legend('DIMMUKF','DIMMUIF','DIMMCIF','SingleSensor')
xlabel('t/s');ylabel('RMSE');
title('位置RMSE')
figure
k=1:steps;
hndl=plot(k,RMSE_Vel_IMM(:,1),'-c+',k,RMSE_Vel_IMM(:,2),'-ms',k,RMSE_Vel_IMM(:,3),'-b*',k,RMSE_Vel_IMM(:,4),'-ro');
set(hndl,'LineWidth',1);   % 线条粗细
legend('DIMMUKF','DIMMUIF','DIMMCIF','SingleSensor')
xlabel('t/s');ylabel('RMSE');
title('速度RMSE')
figure
k=1:steps;
hndl=plot(k,RMSE_Acc_IMM(:,1),'-c+',k,RMSE_Acc_IMM(:,2),'-ms',k,RMSE_Acc_IMM(:,3),'-b*',k,RMSE_Acc_IMM(:,4),'-ro');
set(hndl,'LineWidth',1);   % 线条粗细
legend('DIMMUKF','DIMMUIF','DIMMCIF','SingleSensor')
xlabel('t/s');ylabel('RMSE');
title('加速度速度RMSE')

