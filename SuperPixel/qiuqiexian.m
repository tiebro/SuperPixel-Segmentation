function [x1,y1,x2,y2] = qiuqiexian(r,c);
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
r1 = unique(r);
f = length(r1);
for i = 1:f
    a1 = r1(i,1);%行的各值
    [x1,y1]=find(r==a1);%x1表示行中等于a1的值的集合
    e = length(x1);
    for j = 1:e
        b1 = x1(j,1);
        d1 = c(b1,1);
        D(j) = d1;
    end
    x(i) = max(D);
    y(i) = a1;
    D(:) = 0;
end


x = reshape(x,1*f,1);%变成列向量
y = reshape(y,1*f,1);
xmean = mean(x);
xstd = std(x);
ymean = mean(y);
ystd = std(y);

for i = 1:f
    if ((x(i,1) > xmean-xstd) && (x(i,1) < xmean+xstd)) && ((y(i,1) > ymean-ystd) && (y(i,1) < ymean+ystd))
        xz(i) = x(i,1);
        yz(i) = y(i,1);
    else
        xz(i) = 0;
        yz(i) = 0;
    end
end

[~,g] = find(xz~=0);
f1 = length(g);
for i = 1:f1 %
    xzz(1,i) = xz(1,g(1,i));
    yzz(1,i) = yz(1,g(1,i));
end

[~,k]=size(xzz);
o = min(xzz);
p = max(xzz);

% for n=1:9
n = 1;
    ANSS=polyfit(xzz,yzz,n);  %用polyfit拟合曲线
    for i=1:n+1           %answer矩阵存储每次求得的方程系数，按列存储
       answer(i,n)=ANSS(i);
   end
    x0=o:1:p;
    y0=ANSS(1)*x0.^n    ; %根据求得的系数初始化并构造多项式方程
    for num=2:1:n+1     
        y0=y0+ANSS(num)*x0.^(n+1-num);
    end
    subplot(1,1,n)
    plot(xzz,yzz,'*')
    hold on
    plot(x0,y0)
    
    k = ANSS(1,1);
    b = ANSS(1,2);
%     y = k*x+b;
    x1 = o;y1 = k*x1+b;
    x2 = p;y2 = k*x2+b;

% end
% suptitle('不同次数方程曲线拟合结果，从1到9阶')


end

