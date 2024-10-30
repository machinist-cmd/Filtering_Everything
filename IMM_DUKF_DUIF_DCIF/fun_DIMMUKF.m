function [xk_dimmukf,Pk_dimmukf,xk_dimmukf_wa,A_dimmukf,xk_UKF,Pk_UKF,m_ukf] = fun_DIMMUKF(xk_UKF,Pk_UKF,Pa_ukf,m_ukf,Fk_cv,Gk_cv,Fk_ct,Gk_ct,Fk_ca,Gk_ca,...
    Z_true,Qk1,Qk2,Qk3,sigma_r,sigma_b,sigma_e,xp,M)
%  IMM 分布式无迹卡尔曼滤波融合: DIMMUKF

%%  基于imm-ukf算法，得到每个传感器的局部估计
for m=1:M
% % [xk_m,Pk_m,A_m] = fun_2IMMUKF(xk,Pk,Fk,Gk,Zm(:,m),Qk,sigma_r(m),sigma_b(m), xp(:,m),nx,nz);
[X_immukfm,P_immukfm,X_immukf_wa,xk_UKFm,Pk_UKFm,m_ukfm] = fun_2IMMUKF(xk_UKF(m),Pk_UKF(m),Pa_ukf,m_ukf(:,m),...
    Fk_cv,Gk_cv,Fk_ct,Gk_ct,Fk_ca,Gk_ca,Z_true(:,m),Qk1,Qk2,Qk3,sigma_r(m),sigma_b(m),sigma_e(m),xp(:,m));

xk_UKF{m}=xk_UKFm; Pk_UKF{m}=Pk_UKFm; 

xkm(:,m)=X_immukfm; Pkm(:,:,m)=P_immukfm; xkm_wa(:,m)=X_immukf_wa; m_ukf(:,m)=m_ukfm;
end

%% 基于分布式融合，对所有局部估计进行融合，得到全局估计
invPk_dukf=0;invxk_dukf=0;
for m=1:M
    invPk_dukf=inv(Pkm(:,:,m) ) + invPk_dukf;
    invxk_dukf=inv(Pkm(:,:,m) )*xkm(:,m) + invxk_dukf;
    % % Pk_dukf=inv( inv(Pk_1) +  inv(Pk_2) + inv(Pk_3) + inv(Pk_4) + inv(Pk_5) );
% % xk_dukf=Pk_dukf *  ( inv(Pk_1)*xk_1 +  inv(Pk_2)*xk_2 + inv(Pk_3)*xk_3 + inv(Pk_4)*xk_4 + inv(Pk_5)*xk_5   );
end
% global estimate
Pk_dimmukf=inv(invPk_dukf);
xk_dimmukf=Pk_dimmukf *  invxk_dukf;
xk_dimmukf_wa=( sum(xkm_wa, 2)   )./M;
A_dimmukf= ( sum(m_ukf, 2)   )./M;
A_dimmukf=A_dimmukf./sum(A_dimmukf);

% % A_dimmukf=m_ukf(:,5);
% % A_dukf=Am(m);



end