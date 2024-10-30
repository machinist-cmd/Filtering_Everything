function [xk_a,Pk_a,A]  = fun_UIF(xk_a,Pk_a,Fk,Gk,Z_true,Qk,sigma_r,sigma_b,sigma_e, xp, nx,nz,M)
% 分布式 无迹信息滤波： UIF
%% 
zk=Z_true;% 雷达量测:rm bm em
n=nx;%状态维数，CV6,CT为6 CA为9
nz=nz;% 测量维数
Pk=Pk_a([1:nx],[1:nx]);
xk=xk_a([1:nx]);
% 提出奇异Pk
d=eig(Pk);
if all(d > 0) == 1
    Pk=Pk;
else 
    Pk=diag(diag(Pk)) ;
end
%%
% xkk=Fk*xk;
% Pkk=Fk*Pk*Fk'+Gk*Qk*Gk';
%UT transformation
alpha=0.9;
kk=0;
beta=2; 
lambda=alpha^2*(n+kk)-n;
Wm=[lambda/(lambda+n),  (0.5/(lambda+n))+zeros(1,2*n)];%权值确定
Wc=[lambda/(lambda+n)+1-alpha^2+beta,   (0.5/(lambda+n))+zeros(1,2*n)];
%产生xk的Sigma点
SPk=sqrt(n+lambda)*(chol(Pk))';
Xsigma0=xk;
for i=1:n 
    Xsigma1(:,i)=xk+SPk(:,i);
    Xsigma2(:,i)=xk-SPk(:,i); 
end
Xsigma=[Xsigma0,Xsigma1,Xsigma2];
%产生xkk的Gama点
for i=1:2*n+1
% %     Xgama(:,i)=ProssEq(Xsigma(:,i));
    Xgama(:,i)=Fk*Xsigma(:,i);
end
xkk=Xgama*Wm';
Pkk=zeros(n,n);for i=1:2*n+1; Pkk=Pkk+Wc(i)*((Xgama(:,i)-xkk)*(Xgama(:,i)-xkk)');end; 
Pkk=Pkk + Gk*Qk*Gk';

P_pre_IF_inv = inv(Pkk);                            %%预测p-inv
y_pre_IF = P_pre_IF_inv*xkk;                         %%预测yk+1|k


% 产生xkk的Sigma点
SPkk=sqrt(n+lambda)*(chol(Pkk))';
Zsigma0=xkk;
for i=1:n 
    Zsigma1(:,i)=xkk+SPkk(:,i);
    Zsigma2(:,i)=xkk-SPkk(:,i);
end
Zsigma=[Zsigma0,Zsigma1,Zsigma2];

%M个传感器
for j=1:M
% zkk
for i=1:2*n+1
    [z1,z2,z3] = measurements(Zsigma(:,i), xp(:,j)); 
    Zgama(:,i)=[z1,z2,z3]';
end

zkk=Zgama*Wm';
Rk(:,:,j)=diag([sigma_r(j)^2, sigma_b(j)^2, sigma_e(j)^2]);
Sk=zeros(nz,nz);for i=1:2*n+1; Sk=Sk+Wc(i)*((Zgama(:,i)-zkk)*(Zgama(:,i)-zkk)');end; 
Sk=Sk+Rk(:,:,j);
% % Ck=zeros(nx,nz);for i=1:2*n+1; Ck=Ck+Wc(i)*((Zsigma(:,i)-xkk)*(Zgama(:,i)-zkk)');end;
Ck=zeros(nx,nz);for i=1:2*n+1; Ck=Ck+Wc(i)*((Xgama(:,i)-xkk)*(Zgama(:,i)-zkk)');end;


% 计算分布式信息融合的各个参数：新息
Pxzkk(:,:,j)=Ck;

Ip(:,:,j)=(P_pre_IF_inv*Pxzkk(:,:,j))*inv(Rk(:,:,j))*(P_pre_IF_inv*Pxzkk(:,:,j))';  %% 新信息P

Ix(:,j)=(P_pre_IF_inv*Pxzkk(:,:,j))*inv(Rk(:,:,j))*(zk(:,j)-zkk+Pxzkk(:,:,j)'*y_pre_IF); %% 新信息X

A(j)=(2*pi)^(-1/2)*(det(Sk))^(-1/2)*exp((-1/2)*(zk(:,j)-zkk)'*inv(Sk)*(zk(:,j)-zkk)); %% 模型概率
end

% 分布式新息滤波
Y_update_IF  = y_pre_IF ;
P_update_IF_inv = P_pre_IF_inv ;

%M个传感器，新息叠加
for j =1:M
    Y_update_IF  = Y_update_IF +  Ix(:,j); 
    P_update_IF_inv = P_update_IF_inv + Ip(:,:,j); 
end

A=sum(A)/M;
% A=A(j);

Pk=inv(P_update_IF_inv);   %% p_update
xk = Pk*Y_update_IF;      %% x_update
% Pk
% xk
Pk_a([1:nx],[1:nx])=Pk;
xk_a([1:nx])=xk;


end