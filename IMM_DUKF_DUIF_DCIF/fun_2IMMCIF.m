function [X_immcif,P_immcif,X_immcif_wa,xk_CIF,Pk_CIF,m_cif] = fun_2IMMCIF(xk_CIF,Pk_CIF,Pa_cif,m_cif,...
    Fk_cv,Gk_cv,Fk_ct,Gk_ct,Fk_ca,Gk_ca,Z_true,Qk1,Qk2,Qk3,sigma_r,sigma_b,sigma_e,xp,M)
%% IMM-CIF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
xk_CIFcv=xk_CIF{1}; Pk_CIFcv=Pk_CIF{1}; % cv模型局部滤波器估计
xk_CIFct=xk_CIF{2}; Pk_CIFct=Pk_CIF{2}; % ct模型局部滤波器估计
xk_CIFca=xk_CIF{3}; Pk_CIFca=Pk_CIF{3}; % ca模型局部滤波器估计


%% IMMCIF开始
c_bar1=Pa_cif(1,1)*m_cif(1)+Pa_cif(2,1)*m_cif(2)+Pa_cif(3,1)*m_cif(3);
c_bar2=Pa_cif(1,2)*m_cif(1)+Pa_cif(2,2)*m_cif(2)+Pa_cif(3,2)*m_cif(3);
c_bar3=Pa_cif(1,3)*m_cif(1)+Pa_cif(2,3)*m_cif(2)+Pa_cif(3,3)*m_cif(3);

mu(1,1)=Pa_cif(1,1)*m_cif(1)/c_bar1;
mu(2,1)=Pa_cif(2,1)*m_cif(2)/c_bar1;
mu(3,1)=Pa_cif(3,1)*m_cif(3)/c_bar1;

mu(1,2)=Pa_cif(1,2)*m_cif(1)/c_bar2;
mu(2,2)=Pa_cif(2,2)*m_cif(2)/c_bar2;
mu(3,2)=Pa_cif(3,2)*m_cif(3)/c_bar2;

mu(1,3)=Pa_cif(1,3)*m_cif(1)/c_bar3;
mu(2,3)=Pa_cif(2,3)*m_cif(2)/c_bar3;
mu(3,3)=Pa_cif(3,3)*m_cif(3)/c_bar3;

% 分配各个模型的状态量中共同存在的状态，进行交互
xk_CIF1=xk_CIFcv; xk_CIF2=xk_CIFct; xk_CIF3=xk_CIFca(1:6);
Pk_CIF1=Pk_CIFcv; Pk_CIF2=Pk_CIFct; Pk_CIF3=Pk_CIFca([1:6],[1:6]);

% 输入交互：计算交互后目标于个模型的状态估计和协方差矩阵
X_update_hat1=xk_CIF1*mu(1,1)+xk_CIF2*mu(2,1)+xk_CIF3*mu(3,1);
X_update_hat2=xk_CIF1*mu(1,2)+xk_CIF2*mu(2,2)+xk_CIF3*mu(3,2);
X_update_hat3=xk_CIF1*mu(1,3)+xk_CIF2*mu(2,3)+xk_CIF3*mu(3,3);
P_update_hat1=(Pk_CIF1+(xk_CIF1-X_update_hat1)*(xk_CIF1-X_update_hat1)')*mu(1,1)+(Pk_CIF2+(xk_CIF2-X_update_hat1)*(xk_CIF2-X_update_hat1)')*mu(2,1)+(Pk_CIF3+(xk_CIF3-X_update_hat1)*(xk_CIF3-X_update_hat1)')*mu(3,1);
P_update_hat2=(Pk_CIF1+(xk_CIF1-X_update_hat2)*(xk_CIF1-X_update_hat2)')*mu(1,2)+(Pk_CIF2+(xk_CIF2-X_update_hat2)*(xk_CIF2-X_update_hat2)')*mu(2,2)+(Pk_CIF3+(xk_CIF3-X_update_hat2)*(xk_CIF3-X_update_hat2)')*mu(3,2);
P_update_hat3=(Pk_CIF1+(xk_CIF1-X_update_hat3)*(xk_CIF1-X_update_hat3)')*mu(1,3)+(Pk_CIF2+(xk_CIF2-X_update_hat3)*(xk_CIF2-X_update_hat3)')*mu(2,3)+(Pk_CIF3+(xk_CIF3-X_update_hat3)*(xk_CIF3-X_update_hat3)')*mu(3,3);

%% 
%%% 各个局部滤波器： 分布式容积信息滤波DCIF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%filer1:CV
xk_CIFcv=X_update_hat1; Pk_CIFcv=P_update_hat1;
[xk_CIFcv,Pk_CIFcv,A_CIF1] = fun_CIF(xk_CIFcv,Pk_CIFcv,Fk_cv,Gk_cv,Z_true,Qk1,sigma_r,sigma_b,sigma_e,xp,6,nz,M);

%filer2:CT
xk_CIFct=X_update_hat2; Pk_CIFct=P_update_hat2;
[xk_CIFct,Pk_CIFct,A_CIF2] = fun_CIF(xk_CIFct,Pk_CIFct,Fk_ct,Gk_ct,Z_true,Qk2,sigma_r,sigma_b,sigma_e,xp,6,nz,M);

%filer3：CA
xk_CIFca(1:6)=X_update_hat3; Pk_CIFca([1:6],[1:6])=P_update_hat3;
[xk_CIFca,Pk_CIFca,A_CIF3] = fun_CIF(xk_CIFca,Pk_CIFca,Fk_ca,Gk_ca,Z_true,Qk3,sigma_r,sigma_b,sigma_e,xp,9,nz,M);
%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%分配三个滤波器的估计结果，进行IMM交互输出最终的结果
xk_CIF1=xk_CIFcv; xk_CIF2=xk_CIFct; xk_CIF3=xk_CIFca(1:6);
Pk_CIF1=Pk_CIFcv; Pk_CIF2=Pk_CIFct; Pk_CIF3=Pk_CIFca([1:6],[1:6]);

%似然函数，计算个模型概率
c=A_CIF1*c_bar1+A_CIF2*c_bar2+A_CIF3*c_bar3;
m_cif(1)=A_CIF1*c_bar1/c;
m_cif(2)=A_CIF2*c_bar2/c;
m_cif(3)=A_CIF3*c_bar3/c;

%交互多模型滤波结果
X_immcif=xk_CIF1*m_cif(1)+xk_CIF2*m_cif(2)+xk_CIF3*m_cif(3);
P_immcif=m_cif(1)*(Pk_CIF1+(X_immcif-xk_CIF1)*(X_immcif-xk_CIF1)')+m_cif(2)*(Pk_CIF2+(X_immcif-xk_CIF2)*(X_immcif-xk_CIF2)')+m_cif(3)*(Pk_CIF3+(X_immcif-xk_CIF3)*(X_immcif-xk_CIF3)');

%分配CA模型的加速度估计: CV和CT模型局部滤波器认为目标非加速度运动，即加速度估计为0
xk_CIFcv_acc=zeros(3,1); xk_CIFct_acc=zeros(3,1); xk_CIFca_acc=xk_CIFca(7:9);
X_immcif_wa=xk_CIFcv_acc*m_cif(1) + xk_CIFct_acc*m_cif(2) + xk_CIFca_acc*m_cif(3);

% 重新分配各个局部滤波器估计，进行下一次循环
xk_CIF={xk_CIFcv, xk_CIFct, xk_CIFca};
Pk_CIF={Pk_CIFcv, Pk_CIFct, Pk_CIFca};



end

