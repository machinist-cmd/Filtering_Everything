function [xk_dimmcif,Pk_dimmcif,xk_dimmcif_wa,A_dimmcif,xk_CIF,Pk_CIF,m_cif] = fun_DIMMCIF(xk_CIF,Pk_CIF,Pa_cif,m_cif,Fk_cv,Gk_cv,Fk_ct,Gk_ct,Fk_ca,Gk_ca,...
    Z_true,Qk1,Qk2,Qk3,sigma_r,sigma_b,sigma_e,xp,M)
%  IMM 分布式无迹信息滤波: DIMMCIF

mm=1;
[X_immcif,P_immcif,X_immcif_wa,xk_CIF,Pk_CIF,m_cif] = fun_2IMMCIF(xk_CIF,Pk_CIF,Pa_cif,m_cif,...
    Fk_cv,Gk_cv,Fk_ct,Gk_ct,Fk_ca,Gk_ca,Z_true(:,:),Qk1,Qk2,Qk3,sigma_r(:),sigma_b(:),sigma_e(:),xp(:,:),M);

% global estimate 全局估计
Pk_dimmcif=P_immcif;        % dimmcif估计协方差
xk_dimmcif=X_immcif;        % dimmcif估计结果
xk_dimmcif_wa=X_immcif_wa;  % dimmcif估计结果:加速度
A_dimmcif= m_cif;           % dimmcif模型概率



end



