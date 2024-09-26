% function  main()
close all
clear;
tic
%% read the element and boundary condation
fprintf(1,'generate the modal\n')

a = 5;
b = 10;
numR = 5;
numTheta = 10;
[node, elem, nodeBou, elemBou] = cylinderMesh(a, b, numR, numTheta);
sumNode = size(node,1);
sumElem = size(elem,1);

% set boundary condition
pressure = 1;
BC = {'p', 'dy', 'f', 'dx'};
fixNode = setBCCylinder(BC, nodeBou);
nodeForce = [];

%% material

% -------------------------------------------------------------------------
G0 = 1e4;
K0 = 2e4;
E0 = 9*G0*K0/(3*K0+G0);
nu0 = (3*K0-2*G0)/(2*(3*K0+G0));

mat = [E0, nu0;]; % 弹性模量、泊松比, 多种材料需要换行
ndim = 2;isStress = 0; % 平面应力还是平面应变，为1时为平面应力

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
if sum(ismember(BC,'p'))
    fprintf(1,'trans the pressure\n')
    nodeForceNew = transPresV2(node,elem,ndim,mnode,BC,pressure,elemBou);
    nodeForce = [nodeForce;nodeForceNew];
end

%% Solve elastic
% ----------------------------形成总体刚度阵--------------------------------
fprintf(1,'sort the global K\n')
if ndim == 2
    [GK_u,GK_v,GK_a,B,D] = globalK2D(EX,mu,h,elem,node,isStress,reduce,mnode);
end
% ----------------------------第一类边界条件-----------------------------------
fprintf(1,'boundary condition\n')

% 对角元改1法施加第一类边界条件
[GK,force] = boundary_chang_one(GK_u,GK_v,GK_a,fixNode,nodeForce,sumNode,ndim);

% ----------------------------求解方程--------------------------------------
fprintf(1,'solve the equation\n')
u = GK\force;
% 将计算所得结果改成非稀疏形式
[nodeNdfID,~,rst] = find(u);
u = zeros(ndim*sumNode,1);
u(nodeNdfID) = rst;
u = reshape(u,ndim,[]);
if mnode == 4 && reduce == 0
    [s,mises] = getStress(B,D,u,elem,ndim,sumNode,sumElem,mnode);
elseif mnode == 8 && reduce == 1
    [s,mises] = getStress(B,D,u,elem,ndim,sumNode,sumElem,mnode);
else
    fprintf(1,'不计算应力\n')
end
u = u';

% 显示结果
% ux
figure
axis equal;
showContour(node,elem,u(:,1));
axis off;
% uy
figure
axis equal;
showContour(node,elem,u(:,2));
axis off;
% von-Mises stress
figure
axis equal;
showContour(node,elem,mises);
axis off;
disp('solution is done')

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

[qG_jk, qG_jk1] = deal(zeros(numel(u), length(visG)));
[qK_jk, qK_jk1] = deal(zeros(numel(u), length(visK)));
[uk1, uk2] = deal(zeros(numel(u), 1));
t = t0;
while t < tEnd
    G_dt = (1 - visG.*(1-exp(-dt./tauG))) ;
    K_dt = (1 - visK.*(1-exp(-dt./tauK)));

    if t == t0
        GK_a = GK1_a;
    else
        GK_a = GK1_a + G0*(1 - G_dt)*GK2_a + K0*(1 - K_dt)*GK3_a;
    end
    [GK, feM] = boundary_chang_one(GK_u,GK_v,GK_a,fixNode,nodeForce,sumNode,ndim);

    GK2 = sparse(GK_u,GK_v,GK2_a,sumNode*ndim,sumNode*ndim);
    GK3 = sparse(GK_u,GK_v,GK3_a,sumNode*ndim,sumNode*ndim);

    feM1 = - 2*GK2*(G0*sum(visG.*qG_jk, 2) + 1/2*G0*(1 - G_dt)*uk1);
    feM2 = - GK3*(K0*sum(visK.*qK_jk, 2) + 1/2*K0*(1 - K_dt)*uk1);
    %     force = feT + feM + feM1 + feM2 + feM3;
    force = feM + feM1 + feM2;

    u = GK\force;
    % 将计算所得结果改成非稀疏形式
    [nodeNdfID,~,rst] = find(u);
    u = zeros(ndim*sumNode,1);
    u(nodeNdfID) = rst;

    uStar_k2 = (uk2 + uk1)/2;
    uk2 = uk1;
    uk1 = u;
    
    qG_jk1 = qG_jk;
    qK_jk1 = qK_jk;
    qG_jk = exp(-dt./tauG).*((1-exp(-dt./tauG)).*uStar_k2+qG_jk1);
    qK_jk = exp(-dt./tauK).*((1-exp(-dt./tauK)).*uStar_k2+qK_jk1);

    t = t + dt;
end

uAna = CylinderAnalysis([a, 0], a, b, p, (1-visG)*G0, visG*G0, mu, t0:dt:tEnd, tauG);










