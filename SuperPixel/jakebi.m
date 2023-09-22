%==基本参数
%矩阵阶数n
n = 3;
%计算误差
eps = 0.000005;
%最大迭代次数
Nmax = 10000;
 % 方法一 ： 自行设置矩阵
 A = [2,1,1;
     1,3,1;
     1,1,2];
 B = [3,2,3];
 % 方法二：随机生成矩阵
 % 系数矩阵(默认生成一个对角占优矩阵)

% x = zeros(1,n);
% y = zeros(1,n);
% count = 1;
% 
% while(1)
%      for i = 1:n
%         for j = 1:n
%             %j≠k时
%             if i ~= j
%                 y(i) = y(i) + A(i,j)*x(j);
%             end
%         end
%         y(i)=( B(i) - y(i))/A(i,i);
%     end
%     if max(abs(x-y)) < eps
%         %c为将结果回带到方程组中求解的结果向量，以确认结果无误
%         c = A * y';
%         fprintf('迭代结束，次数%d，最终结果：\n',count);
%         disp(y);
%         disp(c);
%         break;
%     else
%         fprintf('第%d次迭代结果：\n',count);
%         disp(y)
%     end
%     if count == Nmax
%         c = A * y';
%         fprintf('超过最大迭代次数，迭代结束，最终结果：\n');
%         disp(y);
%         disp(c);
%         break; 
%     end
%     count = count + 1;
%     x = y;
%     y(1: n) = 0;
% end

x = zeros(1,n);
y = zeros(1,n);
count = 1;
%中间变量
x1 = zeros(1,n);
x2 = zeros(1,n);

while(1)
    x1(1:n) = 0;
    x2(1:n) = 0;
     for i = 1:n
        for j = 1: i-1
            x1(i) =x1(i) + A(i,j) * y(j);
        end
        for k = i + 1 : n
           x2(i) = x2(i) + A(i,k) * x(k); 
        end
        y(i)=(B(i) - x1(i) - x2(i))/A(i,i);
    end
    if max(abs(x-y)) < eps
        c = A * y';
        fprintf('迭代结束，次数%d，最终结果：\n',count);
        disp(y);
        disp(c);
        break;
    else
        fprintf('第%d次迭代结果：\n',count);
        disp(y)
    end
    if count == Nmax
        c = A * y';
        fprintf('超过最大迭代次数，迭代结束，最终结果：\n');
        disp(y);
        disp(c);
        break; 
    end
    count = count + 1;
    x = y;
end
