function [X_immuif,P_immuif,X_immuif_wa,xk_UIF,Pk_UIF,m_uif] = fun_2IMMUIF(xk_UIF,Pk_UIF,Pa_uif,m_uif,...
    Fk_cv,Gk_cv,Fk_ct,Gk_ct,Fk_ca,Gk_ca,Z_true,Qk1,Qk2,Qk3,sigma_r,sigma_b,sigma_e,xp,M)
%% IMM-UIF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
xk_UIFcv=xk_UIF{1}; Pk_UIFcv=Pk_UIF{1}; % cv模型局部滤波器估计
xk_UIFct=xk_UIF{2}; Pk_UIFct=Pk_UIF{2}; % ct模型局部滤波器估计
xk_UIFca=xk_UIF{3}; Pk_UIFca=Pk_UIF{3}; % ca模型局部滤波器估计


%% IMMUIF开始
c_bar1=Pa_uif(1,1)*m_uif(1)+Pa_uif(2,1)*m_uif(2)+Pa_uif(3,1)*m_uif(3);
c_bar2=Pa_uif(1,2)*m_uif(1)+Pa_uif(2,2)*m_uif(2)+Pa_uif(3,2)*m_uif(3);
c_bar3=Pa_uif(1,3)*m_uif(1)+Pa_uif(2,3)*m_uif(2)+Pa_uif(3,3)*m_uif(3);

mu(1,1)=Pa_uif(1,1)*m_uif(1)/c_bar1;
mu(2,1)=Pa_uif(2,1)*m_uif(2)/c_bar1;
mu(3,1)=Pa_uif(3,1)*m_uif(3)/c_bar1;

mu(1,2)=Pa_uif(1,2)*m_uif(1)/c_bar2;
mu(2,2)=Pa_uif(2,2)*m_uif(2)/c_bar2;
mu(3,2)=Pa_uif(3,2)*m_uif(3)/c_bar2;

mu(1,3)=Pa_uif(1,3)*m_uif(1)/c_bar3;
mu(2,3)=Pa_uif(2,3)*m_uif(2)/c_bar3;
mu(3,3)=Pa_uif(3,3)*m_uif(3)/c_bar3;

% 分配各个模型的状态量中共同存在的状态，进行交互
xk_UIF1=xk_UIFcv; xk_UIF2=xk_UIFct; xk_UIF3=xk_UIFca(1:6);
Pk_UIF1=Pk_UIFcv; Pk_UIF2=Pk_UIFct; Pk_UIF3=Pk_UIFca([1:6],[1:6]);

% 输入交互：计算交互后目标于个模型的状态估计和协方差矩阵
X_update_hat1=xk_UIF1*mu(1,1)+xk_UIF2*mu(2,1)+xk_UIF3*mu(3,1);
X_update_hat2=xk_UIF1*mu(1,2)+xk_UIF2*mu(2,2)+xk_UIF3*mu(3,2);
X_update_hat3=xk_UIF1*mu(1,3)+xk_UIF2*mu(2,3)+xk_UIF3*mu(3,3);
P_update_hat1=(Pk_UIF1+(xk_UIF1-X_update_hat1)*(xk_UIF1-X_update_hat1)')*mu(1,1)+(Pk_UIF2+(xk_UIF2-X_update_hat1)*(xk_UIF2-X_update_hat1)')*mu(2,1)+(Pk_UIF3+(xk_UIF3-X_update_hat1)*(xk_UIF3-X_update_hat1)')*mu(3,1);
P_update_hat2=(Pk_UIF1+(xk_UIF1-X_update_hat2)*(xk_UIF1-X_update_hat2)')*mu(1,2)+(Pk_UIF2+(xk_UIF2-X_update_hat2)*(xk_UIF2-X_update_hat2)')*mu(2,2)+(Pk_UIF3+(xk_UIF3-X_update_hat2)*(xk_UIF3-X_update_hat2)')*mu(3,2);
P_update_hat3=(Pk_UIF1+(xk_UIF1-X_update_hat3)*(xk_UIF1-X_update_hat3)')*mu(1,3)+(Pk_UIF2+(xk_UIF2-X_update_hat3)*(xk_UIF2-X_update_hat3)')*mu(2,3)+(Pk_UIF3+(xk_UIF3-X_update_hat3)*(xk_UIF3-X_update_hat3)')*mu(3,3);

%% 
%%% 各个局部滤波器： 分布式无迹信息滤波DUIF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%filer1:CV
xk_UIFcv=X_update_hat1; Pk_UIFcv=P_update_hat1;
[xk_UIFcv,Pk_UIFcv,A_UIF1] = fun_UIF(xk_UIFcv,Pk_UIFcv,Fk_cv,Gk_cv,Z_true,Qk1,sigma_r,sigma_b,sigma_e,xp,6,nz,M);

%filer2:CT
xk_UIFct=X_update_hat2; Pk_UIFct=P_update_hat2;
[xk_UIFct,Pk_UIFct,A_UIF2] = fun_UIF(xk_UIFct,Pk_UIFct,Fk_ct,Gk_ct,Z_true,Qk2,sigma_r,sigma_b,sigma_e,xp,6,nz,M);

%filer3：CA
xk_UIFca(1:6)=X_update_hat3; Pk_UIFca([1:6],[1:6])=P_update_hat3;
[xk_UIFca,Pk_UIFca,A_UIF3] = fun_UIF(xk_UIFca,Pk_UIFca,Fk_ca,Gk_ca,Z_true,Qk3,sigma_r,sigma_b,sigma_e,xp,9,nz,M);
%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%分配三个滤波器的估计结果，进行IMM交互输出最终的结果
xk_UIF1=xk_UIFcv; xk_UIF2=xk_UIFct; xk_UIF3=xk_UIFca(1:6);
Pk_UIF1=Pk_UIFcv; Pk_UIF2=Pk_UIFct; Pk_UIF3=Pk_UIFca([1:6],[1:6]);

%似然函数，计算个模型概率
c=A_UIF1*c_bar1+A_UIF2*c_bar2+A_UIF3*c_bar3;
m_uif(1)=A_UIF1*c_bar1/c;
m_uif(2)=A_UIF2*c_bar2/c;
m_uif(3)=A_UIF3*c_bar3/c;

%交互多模型滤波结果
X_immuif=xk_UIF1*m_uif(1)+xk_UIF2*m_uif(2)+xk_UIF3*m_uif(3);
P_immuif=m_uif(1)*(Pk_UIF1+(X_immuif-xk_UIF1)*(X_immuif-xk_UIF1)')+m_uif(2)*(Pk_UIF2+(X_immuif-xk_UIF2)*(X_immuif-xk_UIF2)')+m_uif(3)*(Pk_UIF3+(X_immuif-xk_UIF3)*(X_immuif-xk_UIF3)');

%分配CA模型的加速度估计: CV和CT模型局部滤波器认为目标非加速度运动，即加速度估计为0
xk_UIFcv_acc=zeros(3,1); xk_UIFct_acc=zeros(3,1); xk_UIFca_acc=xk_UIFca(7:9);
X_immuif_wa=xk_UIFcv_acc*m_uif(1) + xk_UIFct_acc*m_uif(2) + xk_UIFca_acc*m_uif(3);

% 重新分配各个局部滤波器估计，进行下一次循环
xk_UIF={xk_UIFcv, xk_UIFct, xk_UIFca};
Pk_UIF={Pk_UIFcv, Pk_UIFct, Pk_UIFca};



end

