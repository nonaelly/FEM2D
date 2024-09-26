% function  main()
clear;
tic
%% read the element and boundary condation 
fprintf(1,'read the modal\n')
node = load('NLIST.DAT');
sumNode = size(node,1);
elem = load('ELIST.DAT');
sumElem = size(elem,1);
fixNode = load('fixNode.dat');
nodeForce = load('nodeForce.dat'); % 节点力
nodeForce(1,:) = [];

%% material

% -------------------------------------------------------------------------
mat = [200000,0.3;]; % 弹性模量、泊松比, 多种材料需要换行
ndim = 2;isStress = 1; % 平面应力还是平面应变，为1时为平面应力

h = 1;   % 应用于二维问题，默认是1

matID  = elem(:,2); % elem的第二列为材料编号
EX = mat(matID,1); % 每个单元的杨氏模量
mu = mat(matID,2); % 每个单元的泊松比

elem(:,1:2) = [];
node(:,1)   = [];
node = node(:,1:ndim);

mnode = size(elem,2);  % 单元类型

if mnode == 4
    reduce = 0;  % 是否采用减缩积分，为1表示就采用减缩积分，0表示全积分。
else
    reduce = 1;  
end

% ---------------------------转化压强---------------------------------------
if exist('press.dat','file')
    fprintf(1,'trans the pressure\n')
    nodeForceNew = transPres(node,elem,ndim,mnode);
    nodeForce = [nodeForce;nodeForceNew];
end

%% Compute viscoelastic StiffMatrix K1 K2
fprintf(1,'sort the global K1, K2 for vis\n')
if ndim == 2
    [GK_u,GK_v,GK1_a,GK2_a,GK3_a,B,D] = globalKVis2D(EX,mu,h,elem,node,isStress,reduce,mnode);
end

%% Solve viscoelastic process

t0=0;   % start time
tEnd=300;  % total time
dt=0.1; % time step

tauG = 1/2.3979;  % relaxation times
tauK = 1/2.3979;  % relaxation times
visG = 0.99; % nomalized shear modulus
visK = 0;    % nomalized volume modulus

K0 = EX(1)./(3*(1-2*mu(1)));
G0 = EX(1)./(2*(1+mu(1)));

[qG_jk, qG_jk1] = deal(zeros(length(u), length(visG)));
[qK_jk, qK_jk1] = deal(zeros(length(u), length(visK)));

t = 0;
while t < tEnd
    G_dt = G0*(1 - visG.*(1-exp(-dt./tauG))) ;
    K_dt = K0*(1 - visK.*(1-exp(-dt./tauK)));

    GK_a = GK1_a + G0*(1 - G_dt)*GK2_a + K0*(1 - K_dt)*GK3_a;
    [GK, feM] = boundary_chang_one(GK_u,GK_v,GK_a,fixNode,nodeForce,sumNode,ndim);
    feM1 = - 2*GK2_a*(G0*sum(visG.*qG_jk, 2) + 1/2*G0*(1 - G_dt)*uk1);
    feM2 = - 2*GK3_a*(K0*sum(visK.*qK_jk, 2) + 1/2*K0*(1 - K_dt)*uk1);
    %     force = feT + feM + feM1 + feM2 + feM3;
    force = feM + feM1 + feM2;

    u = GK\force;

    uk2 = uk1;
    uk1 = u;
    uStar_k2 = (uk2 + uk1)/2;

    qG_jk1 = qG_jk;
    qK_jk1 = qK_jk;
    qG_jk = exp(-dt./tauG).*((1-exp(-dt./tauG)).*uStar_k2+qG_jk1);
    qK_jk = exp(-dt./tauK).*((1-exp(-dt./tauK)).*uStar_k2+qK_jk1);

    t = t + dt;
end












