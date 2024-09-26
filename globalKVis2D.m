function [u,v,a1,a2,a3,B,D] = globalKVis2D(elemEX,elemMu,h,elem,node,isStress,reduce,mnode)
% Stiffness matrix for viscoelastic problem
% K1 : elastic

% 生成总体刚度矩阵
% 适用于二维问题
% 总体刚度阵GK采用行-列-值(u-v-a)的稀疏矩阵形式存储
% sumNode = size(node,1);
sumElem = size(elem,1);
% 判断是否输出应力
if mnode == 4
    if reduce == 0
        B = zeros(3*4*sumElem,8);
    else
        B = 0;
    end
elseif mnode == 8
    if reduce == 1
        B = zeros(3*8*sumElem,16);
    else
        B = 0;
    end
end
D = zeros(sumElem*3,3);

u = zeros(sumElem*(mnode*2)*(mnode*2),1); % 稀疏矩阵的行-列-值
v = zeros(sumElem*(mnode*2)*(mnode*2),1);
a1 = zeros(sumElem*(mnode*2)*(mnode*2),1);
a2 = zeros(sumElem*(mnode*2)*(mnode*2),1);
a3 = zeros(sumElem*(mnode*2)*(mnode*2),1);

ndim = 2;
K0 = elemEX./(3*(1-2*elemMu));
G0 = elemEX./(2*(1+elemMu));

% GK = zeros(sumNode*2);
if mnode == 4
    for n = 1:sumElem
        EX = elemEX(n);
        mu = elemMu(n);
        x1 = node(elem(n,1),1);
        y1 = node(elem(n,1),2);
        x2 = node(elem(n,2),1);
        y2 = node(elem(n,2),2);
        x3 = node(elem(n,3),1);
        y3 = node(elem(n,3),2);
        x4 = node(elem(n,4),1);
        y4 = node(elem(n,4),2);
        [~, D2, D3] = PlaneData(ndim, G0(n), K0(n));
        % elastic       
        [K1,elemB1,elemD1] = elemK2D4(EX,mu,h,x1,y1,x2,y2,x3,y3,x4,y4,isStress,reduce);  % 生成单元刚度阵
        % visco : G 
        % K2 = 1/2*B'*D2*B
        K2 = elemK2D4Vis(h,x1,y1,x2,y2,x3,y3,x4,y4,reduce,D2/2);  % 生成单元刚度阵
        % visco : K 
        % K3 = B'*D3*B
        K3 = elemK2D4Vis(h,x1,y1,x2,y2,x3,y3,x4,y4,reduce,D3);  % 生成单元刚度阵
        if reduce == 0
            B(12*n-11:12*n,:) = elemB1;
            D(3*n-2:3*n,:) = elemD1;
        end
        nodeID = elem(n,:);
        ndfID([1,3,5,7]) = nodeID*2-1;
        ndfID([2,4,6,8]) = nodeID*2;
        ddd = repmat(ndfID,8,1); % ddd 为临时变量，无实际意义
        u(64*(n-1)+1:64*n) = ddd(:);
        v(64*(n-1)+1:64*n) = repmat(ndfID,1,8);
        a1(64*(n-1)+1:64*n) = K1(:);
        a2(64*(n-1)+1:64*n) = K2(:);
        a3(64*(n-1)+1:64*n) = K3(:);
    end
% GK = sparse(u,v,a,sumNode*2,sumNode*2);
else % mnode == 8
    x = zeros(8,1);
    y = zeros(8,1);
    for n = 1:sumElem

        EX = elemEX(n);
        mu = elemMu(n);
        for m = 1:8
            x(m) = node(elem(n,m),1);
            y(m) = node(elem(n,m),2);
        end
        [K1,elemB1,elemD1] = elemK2D8(EX,mu,h,x,y,isStress,reduce);  % 生成单元刚度阵
        if reduce == 1
            B(12*n-11:12*n,:) = elemB1;
            D(3*n-2:3*n,:) = elemD1;
        end
        nodeID = elem(n,:);
        ndfID([1,3,5,7,9,11,13,15]) = nodeID*2-1;
        ndfID([2,4,6,8,10,12,14,16]) = nodeID*2;
        ddd = repmat(ndfID,16,1); % ddd 为临时变量，无实际意义
        u(256*(n-1)+1:256*n) = ddd(:);
        v(256*(n-1)+1:256*n) = repmat(ndfID,1,16);
        a1(256*(n-1)+1:256*n) = K1(:);
    end
end










