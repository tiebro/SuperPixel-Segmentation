%==��������
%�������n
n = 3;
%�������
eps = 0.000005;
%����������
Nmax = 10000;
 % ����һ �� �������þ���
 A = [2,1,1;
     1,3,1;
     1,1,2];
 B = [3,2,3];
 % ��������������ɾ���
 % ϵ������(Ĭ������һ���Խ�ռ�ž���)

% x = zeros(1,n);
% y = zeros(1,n);
% count = 1;
% 
% while(1)
%      for i = 1:n
%         for j = 1:n
%             %j��kʱ
%             if i ~= j
%                 y(i) = y(i) + A(i,j)*x(j);
%             end
%         end
%         y(i)=( B(i) - y(i))/A(i,i);
%     end
%     if max(abs(x-y)) < eps
%         %cΪ������ش��������������Ľ����������ȷ�Ͻ������
%         c = A * y';
%         fprintf('��������������%d�����ս����\n',count);
%         disp(y);
%         disp(c);
%         break;
%     else
%         fprintf('��%d�ε��������\n',count);
%         disp(y)
%     end
%     if count == Nmax
%         c = A * y';
%         fprintf('�����������������������������ս����\n');
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
%�м����
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
        fprintf('��������������%d�����ս����\n',count);
        disp(y);
        disp(c);
        break;
    else
        fprintf('��%d�ε��������\n',count);
        disp(y)
    end
    if count == Nmax
        c = A * y';
        fprintf('�����������������������������ս����\n');
        disp(y);
        disp(c);
        break; 
    end
    count = count + 1;
    x = y;
end
