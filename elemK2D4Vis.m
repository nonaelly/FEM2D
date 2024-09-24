function [K2] = elemK2D4Vis(h,x1,y1,x2,y2,x3,y3,x4,y4,reduce,D)
% 计算单元刚度矩阵
% 四节点四边形单元
% 采用二点高斯积分，全积分高斯点位置+0.577350269189626,-0.577350269189626，权系数为1
% 若采取减缩积分，积分点位置为0，权系数为2
% 可以采用更高阶的高斯积分
% mnode == 4
% 4----3
% |    |
% |    |
% 1----2

K2 = zeros(8,8);
if reduce == 0
   elemB = zeros(12,8);
else
    elemB = 0;
end
 
if reduce == 0 % 采用完全积分
    nip = 2; % 单维方向的积分点数
    x = [-0.577350269189626,0.577350269189626]; % 积分点
    w = [1,1]; % 权系数
else
    nip = 1;
    x = 0; % 积分点
    w = 2; % 权系数
end


for m = 1:nip
    for n = 1:nip
        [J,B] = elemB2D4(x1,y1,x2,y2,x3,y3,x4,y4,x(n),x(m));
        if reduce == 0
            elemB(6*m+3*n-8:6*m+3*n-6,:) = B;
        end
        K2 = K2 + w(m)*w(n)*B'*D*B*det(J)*h;
    end
end


