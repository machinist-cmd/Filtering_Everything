function [xk_a,Pk_a,A]  = fun_CIF(xk_a,Pk_a,Fk,Gk,Z_true,Qk,sigma_r,sigma_b,sigma_e, xp, nx,nz,M)
% 分布式 容积信息滤波： CIF
%% 
zk=Z_true;% 雷达量测:rm bm
n=nx;%状态维数，CV4,CT为5 CA为6
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
mx=2*n;
%%%%%%%%%%%%%%%%%%%%%%%
Ix=[eye(nx), -eye(nx)];
xi=zeros(nx, mx);
for i=1:mx
    xi(:,i)=sqrt(mx/2)*Ix(:,i);
end
SPk=chol(Pk)';
%。产生容积点，通过目标函数传播容积点
for i=1:mx
    Xsigma(:,i)=xk+SPk*xi(:,i);
    Xgama(:,i)=Fk*Xsigma(:,i);
end
xkk=sum(Xgama,2)/mx;
Pkk_i=zeros(nx,nx,mx);
for i=1:mx
    Pkk_i(:,:,i)=(Xgama(:,i)-xkk)*(Xgama(:,i)-xkk)';
end
Pkk=sum(Pkk_i, 3)/mx + Gk*Qk*Gk';

% 
P_pre_IF_inv = inv(Pkk);                            %%预测p-inv
y_pre_IF = P_pre_IF_inv*xkk;                         %%预测yk+1|k


%产生xkk的容积点
SPkk=chol(Pkk)';
for i=1:mx 
    Zsigma(:,i)=xkk+SPkk*xi(:,i);
end


%M个传感器
for j=1:M
% zkk
for i=1:mx
    [z1,z2,z3] = measurements(Zsigma(:,i), xp(:,j)); 
    Zgama(:,i)=[z1,z2,z3]';
end

zkk=sum(Zgama,2)/mx;
Rk(:,:,j)=diag([sigma_r(j)^2, sigma_b(j)^2, sigma_e(j)^2]);
Sk_i=zeros(nz,nz,mx);
for i=1:mx
    Sk_i(:,:,i)=(Zgama(:,i)-zkk)*(Zgama(:,i)-zkk)';
end
Sk=sum(Sk_i, 3)/mx + Rk(:,:,j);
Ck_i=zeros(nx,nz,mx);for i=1:mx; Ck_i(:,:,i)=(Zsigma(:,i)-xkk)*(Zgama(:,i)-zkk)';end;
Ck=sum(Ck_i, 3)/mx;

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