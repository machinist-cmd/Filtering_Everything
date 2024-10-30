function [X_immukf,P_immukf,X_immukf_wa,xk_UKF,Pk_UKF,m_ukf] = fun_2IMMUKF(xk_UKF,Pk_UKF,Pa_ukf,m_ukf,...
    Fk_cv,Gk_cv,Fk_ct,Gk_ct,Fk_ca,Gk_ca,Z_true,Qk1,Qk2,Qk3,sigma_r,sigma_b,sigma_e,xp)
%% IMM-UKF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% 测量方差设置
lambda_m=1;
sigma_r=sigma_r*lambda_m;
sigma_b=sigma_b*lambda_m;
% 过程方差设置
lambda_p=1;
% % Qk1=Qk1*lambda_p;
% % Qk2=Qk2*lambda_p;
% % Qk3=Qk3*lambda_p;
Qk1(logical(eye(size(Qk1))))=diag(Qk1)*lambda_p;
Qk2(logical(eye(size(Qk2))))=diag(Qk2)*lambda_p;
Qk3(logical(eye(size(Qk3))))=diag(Qk3)*lambda_p;

nz=size(Z_true(:,:),1);% 每个雷达量测维数
% 估计变量分配
xk_UKFcv=xk_UKF{1}{1}; Pk_UKFcv=Pk_UKF{1}{1}; % cv模型局部滤波器估计
xk_UKFct=xk_UKF{1}{2}; Pk_UKFct=Pk_UKF{1}{2}; % ct模型局部滤波器估计
xk_UKFca=xk_UKF{1}{3}; Pk_UKFca=Pk_UKF{1}{3}; % ca模型局部滤波器估计


%% IMMUKF开始
c_bar1=Pa_ukf(1,1)*m_ukf(1)+Pa_ukf(2,1)*m_ukf(2)+Pa_ukf(3,1)*m_ukf(3);
c_bar2=Pa_ukf(1,2)*m_ukf(1)+Pa_ukf(2,2)*m_ukf(2)+Pa_ukf(3,2)*m_ukf(3);
c_bar3=Pa_ukf(1,3)*m_ukf(1)+Pa_ukf(2,3)*m_ukf(2)+Pa_ukf(3,3)*m_ukf(3);

mu(1,1)=Pa_ukf(1,1)*m_ukf(1)/c_bar1;
mu(2,1)=Pa_ukf(2,1)*m_ukf(2)/c_bar1;
mu(3,1)=Pa_ukf(3,1)*m_ukf(3)/c_bar1;

mu(1,2)=Pa_ukf(1,2)*m_ukf(1)/c_bar2;
mu(2,2)=Pa_ukf(2,2)*m_ukf(2)/c_bar2;
mu(3,2)=Pa_ukf(3,2)*m_ukf(3)/c_bar2;

mu(1,3)=Pa_ukf(1,3)*m_ukf(1)/c_bar3;
mu(2,3)=Pa_ukf(2,3)*m_ukf(2)/c_bar3;
mu(3,3)=Pa_ukf(3,3)*m_ukf(3)/c_bar3;

% 分配各个模型的状态量中共同存在的状态，进行交互
xk_UKF1=xk_UKFcv; xk_UKF2=xk_UKFct; xk_UKF3=xk_UKFca(1:6);
Pk_UKF1=Pk_UKFcv; Pk_UKF2=Pk_UKFct; Pk_UKF3=Pk_UKFca([1:6],[1:6]);

% 输入交互：计算交互后目标于个模型的状态估计和协方差矩阵
X_update_hat1=xk_UKF1*mu(1,1)+xk_UKF2*mu(2,1)+xk_UKF3*mu(3,1);
X_update_hat2=xk_UKF1*mu(1,2)+xk_UKF2*mu(2,2)+xk_UKF3*mu(3,2);
X_update_hat3=xk_UKF1*mu(1,3)+xk_UKF2*mu(2,3)+xk_UKF3*mu(3,3);
P_update_hat1=(Pk_UKF1+(xk_UKF1-X_update_hat1)*(xk_UKF1-X_update_hat1)')*mu(1,1)+(Pk_UKF2+(xk_UKF2-X_update_hat1)*(xk_UKF2-X_update_hat1)')*mu(2,1)+(Pk_UKF3+(xk_UKF3-X_update_hat1)*(xk_UKF3-X_update_hat1)')*mu(3,1);
P_update_hat2=(Pk_UKF1+(xk_UKF1-X_update_hat2)*(xk_UKF1-X_update_hat2)')*mu(1,2)+(Pk_UKF2+(xk_UKF2-X_update_hat2)*(xk_UKF2-X_update_hat2)')*mu(2,2)+(Pk_UKF3+(xk_UKF3-X_update_hat2)*(xk_UKF3-X_update_hat2)')*mu(3,2);
P_update_hat3=(Pk_UKF1+(xk_UKF1-X_update_hat3)*(xk_UKF1-X_update_hat3)')*mu(1,3)+(Pk_UKF2+(xk_UKF2-X_update_hat3)*(xk_UKF2-X_update_hat3)')*mu(2,3)+(Pk_UKF3+(xk_UKF3-X_update_hat3)*(xk_UKF3-X_update_hat3)')*mu(3,3);
%%%%% 各个局部滤波器 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%filer1:CV
xk_UKFcv=X_update_hat1; Pk_UKFcv=P_update_hat1;
[xk_UKFcv,Pk_UKFcv,A_UKF1] = fun_UKF(xk_UKFcv,Pk_UKFcv,Fk_cv,Gk_cv,Z_true,Qk1,sigma_r,sigma_b,sigma_e,xp,6,nz);

%filer2:CT
xk_UKFct=X_update_hat2; Pk_UKFct=P_update_hat2;
[xk_UKFct,Pk_UKFct,A_UKF2] = fun_UKF(xk_UKFct,Pk_UKFct,Fk_ct,Gk_ct,Z_true,Qk2,sigma_r,sigma_b,sigma_e,xp,6,nz);

%filer3：CA
xk_UKFca(1:6)=X_update_hat3; Pk_UKFca([1:6],[1:6])=P_update_hat3;
[xk_UKFca,Pk_UKFca,A_UKF3] = fun_UKF(xk_UKFca,Pk_UKFca,Fk_ca,Gk_ca,Z_true,Qk3,sigma_r,sigma_b,sigma_e,xp,9,nz);
%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%分配三个滤波器的估计结果，进行IMM交互输出最终的结果
xk_UKF1=xk_UKFcv; xk_UKF2=xk_UKFct; xk_UKF3=xk_UKFca(1:6);
Pk_UKF1=Pk_UKFcv; Pk_UKF2=Pk_UKFct; Pk_UKF3=Pk_UKFca([1:6],[1:6]);

%似然函数，计算个模型概率
c=A_UKF1*c_bar1+A_UKF2*c_bar2+A_UKF3*c_bar3;
m_ukf(1)=A_UKF1*c_bar1/c;
m_ukf(2)=A_UKF2*c_bar2/c;
m_ukf(3)=A_UKF3*c_bar3/c;

%交互多模型滤波结果
X_immukf=xk_UKF1*m_ukf(1)+xk_UKF2*m_ukf(2)+xk_UKF3*m_ukf(3);
P_immukf=m_ukf(1)*(Pk_UKF1+(X_immukf-xk_UKF1)*(X_immukf-xk_UKF1)')+m_ukf(2)*(Pk_UKF2+(X_immukf-xk_UKF2)*(X_immukf-xk_UKF2)')+m_ukf(3)*(Pk_UKF3+(X_immukf-xk_UKF3)*(X_immukf-xk_UKF3)');

%分配CA模型的加速度估计: CV和CT模型局部滤波器认为目标非加速度运动，即加速度估计为0
xk_UKFcv_acc=zeros(3,1); xk_UKFct_acc=zeros(3,1); xk_UKFca_acc=xk_UKFca(7:9);
X_immukf_wa=xk_UKFcv_acc*m_ukf(1) + xk_UKFct_acc*m_ukf(2) + xk_UKFca_acc*m_ukf(3);

% 重新分配各个局部滤波器估计，进行下一次循环
xk_UKF={xk_UKFcv, xk_UKFct, xk_UKFca};
Pk_UKF={Pk_UKFcv, Pk_UKFct, Pk_UKFca};



end

