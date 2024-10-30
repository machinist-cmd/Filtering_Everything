function [xk_dimmuif,Pk_dimmuif,xk_dimmuif_wa,A_dimmuif,xk_UIF,Pk_UIF,m_uif] = fun_DIMMUIF(xk_UIF,Pk_UIF,Pa_uif,m_uif,Fk_cv,Gk_cv,Fk_ct,Gk_ct,Fk_ca,Gk_ca,...
    Z_true,Qk1,Qk2,Qk3,sigma_r,sigma_b,sigma_e,xp,M)
%  IMM 分布式无迹信息滤波: DIMMUIF

mm=1;
[X_immuif,P_immuif,X_immuif_wa,xk_UIF,Pk_UIF,m_uif] = fun_2IMMUIF(xk_UIF,Pk_UIF,Pa_uif,m_uif,...
    Fk_cv,Gk_cv,Fk_ct,Gk_ct,Fk_ca,Gk_ca,Z_true(:,:),Qk1,Qk2,Qk3,sigma_r(:),sigma_b(:),sigma_e(:),xp(:,:),M);

% 全局估计
Pk_dimmuif=P_immuif;            % dimmuif估计协方差
xk_dimmuif=X_immuif;            % dimmuif估计结果
xk_dimmuif_wa=X_immuif_wa;      % dimmuif估计结果:加速度
A_dimmuif= m_uif;               % dimmuif模型概率



end



