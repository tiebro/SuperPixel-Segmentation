function [x1,y1,x2,y2] = qiuqiexian(r,c);
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
r1 = unique(r);
f = length(r1);
for i = 1:f
    a1 = r1(i,1);%�еĸ�ֵ
    [x1,y1]=find(r==a1);%x1��ʾ���е���a1��ֵ�ļ���
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


x = reshape(x,1*f,1);%���������
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
    ANSS=polyfit(xzz,yzz,n);  %��polyfit�������
    for i=1:n+1           %answer����洢ÿ����õķ���ϵ�������д洢
       answer(i,n)=ANSS(i);
   end
    x0=o:1:p;
    y0=ANSS(1)*x0.^n    ; %������õ�ϵ����ʼ�����������ʽ����
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
% suptitle('��ͬ��������������Ͻ������1��9��')


end

