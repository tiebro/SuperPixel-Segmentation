% Segment based on area, Region Growing;
clear;clc;close all;
[fileName,pathName] = uigetfile('*.*','Please select an image');%�ļ���ѡ���ļ�
if(fileName)
    fileName = strcat(pathName,fileName);
    fileName = lower(fileName);%һ�µ�Сд��ĸ��ʽ
else 
    J = 0;%��¼�����������ָ�õ�������
    msgbox('Please select an image');
    return; %�˳�����
end
 
img = imread(fileName);
I = img;
if( ~( size(img,3)-3 ))
    I = rgb2gray(I);%ת��Ϊ��ͨ���Ҷ�ͼ
end
I = im2double(I); %ͼ��Ҷ�ֵ��һ����[0,1]֮��
% Ireshape = imresize(I,[600,800]);
% I = Ireshape(51:475,200:699);
gausFilter = fspecial('gaussian',[5 5],0.5);%fspecial����Ԥ������˲�����   'gaussian'��Ϊ��˹��ͨ�˲���
I = imfilter(I,gausFilter,'replicate');%replicateͼ���Сͨ��������߽��ֵ����չ
 
%���ӵ�Ľ���ʽѡ��
if( exist('x','var') == 0 && exist('y','var') == 0)
    figure,imshow(I,[]);
    hold on;
    [y,x] = getpts;%���ȡ��  �س�ȷ��
    x = round(x(1));%ѡ�����ӵ�
    y = round(y(1));
end
 
% if( nargin == 0)
%     reg_maxdist = 0.1;
%     %nargin��matlab�����д�г��õ�һ�����ɣ���Ҫ���ڼ��㵱ǰ�����������������
%     %����һ����Ը���nargin�ķ���ֵ��ȷ�����������������ȱʡֵ����ʵ���У����
%     %�û�����Ĳ�������Ϊ�㣬��ôĬ��Ϊ0.2
% end
reg_maxdist = 0.2;
J = zeros(size(I)); % �������ķ���ֵ����¼�����������õ�������
Isizes = size(I);
reg_mean = I(x,y);%��ʾ�ָ�õ������ڵ�ƽ��ֵ����ʼ��Ϊ���ӵ�ĻҶ�ֵ
reg_size = 1;%�ָ�ĵ������򣬳�ʼ��ֻ�����ӵ�һ��
neg_free = 10000; %��̬�����ڴ��ʱ��ÿ������������ռ��С
neg_list = zeros(neg_free,3);
%���������б�����Ԥ�ȷ������ڴ�������������ص������ֵ�ͻҶ�ֵ�Ŀռ䣬����
%���ͼ��Ƚϴ���Ҫ���neg_free��ʵ��matlab�ڴ�Ķ�̬����
neg_pos = 0;%���ڼ�¼neg_list�еĴ����������ص�ĸ���
pixdist = 0;
%��¼�������ص����ӵ��ָ������ľ�����
%��һ�δ��������ĸ��������ص�͵�ǰ���ӵ�ľ���
%�����ǰ����Ϊ��x,y����ôͨ��neigb���ǿ��Եõ����ĸ��������ص�λ��
neigb = [ -1 0;
          1  0;
          0 -1;
          0  1];
 %��ʼ�������������������д��������������ص���Ѿ��ָ�õ��������ص�ĻҶ�ֵ����
 %����reg_maxdis,������������
 
 while (pixdist < 0.06 && reg_size < numel(I)) %numel��Ԫ�ظ���
     %�����µ��������ص�neg_list��
     for j=1:4
         xn = x + neigb(j,1);
         yn = y + neigb(j,2);
         %������������Ƿ񳬹���ͼ��ı߽�
         ins = (xn>=1)&&(yn>=1)&&(xn<=Isizes(1))&&(yn<=Isizes(1));
         %�������������ͼ���ڲ���������δ�ָ�ã���ô������ӵ������б���
         if( ins && J(xn,yn)==0)
             neg_pos = neg_pos+1;
             neg_list(neg_pos,:) =[ xn, yn, I(xn,yn)];%�洢��Ӧ��ĻҶ�ֵ
             J(xn,yn) = 1;%��ע���������ص��Ѿ������ʹ� ������ζ�ţ����ڷָ�������
         end
     end
    %���������ڴ���ʲ����������µ��ڴ�ռ�
    if (neg_pos+10>neg_free)
        neg_free = neg_free + 100000;
        neg_list((neg_pos +1):neg_free,:) = 0;
    end
    %�����д����������ص���ѡ��һ�����ص㣬�õ�ĻҶ�ֵ���Ѿ��ָ������ҶȾ�ֵ��
    %��ľ���ֵʱ����������������С��
    dist = abs(neg_list(1:neg_pos,3)-reg_mean);
    [pixdist,index] = min(dist);%pixdist��¼ÿ�е���Сֵ��index��¼ÿ����Сֵ���к�
%     if dist <= reg_maxdist
    %����������µľ�ֵ
    reg_mean = (reg_mean * reg_size +neg_list(index,3))/(reg_size + 1);
    reg_size = reg_size + 1;
    %���ɵ����ӵ���Ϊ�Ѿ��ָ�õ��������ص�
    J(x,y)=2;%��־�����ص��Ѿ��Ƿָ�õ����ص�
    x = neg_list(index,1);
    y = neg_list(index,2);
    
%     pause(0);%��̬����
%     [x1,y1] = find(J==2);
%     plot(y1,x1,'r.');
    
    %���µ����ӵ�Ӵ����������������б����Ƴ�
    neg_list(index,:) = neg_list(neg_pos,:);
    neg_pos = neg_pos -1;
%     else
%         break
%     end
 end
 
 J = (J==2);%����֮ǰ���ָ�õ����ص���Ϊ2
 hold off;
 figure,imshow(J);
 J = bwmorph(J,'dilate');%����ն�dilate
 figure,imshow(J);
 figure,imshow(I+J);
 

hold on 
%������
[r,c] = size(J);
[x2,y2] = find(J==1);
A = unique(x2);
a = length(A);
D = [];
for i = 1:a
    b = A(i,1);%������
    value = 0;
    k = 0;
    for j = 1:c
        if J(b,j) == 1
            value = value + j;
            k = k+1;
        end
    end
    g = round(value/k);
    plot(g,b,'r.');
    D(i,1) = b;
    D(i,2) = g;
end
B = D;
B(:,[2]) = [];
B = B';
C = D;
C(:,[1]) = [];
C = C';

[~,e]=size(B);
max = max(B);
for n=1:3
    ANSS=polyfit(B,C,n);  %��polyfit�������
    for i=1:n+1           %answer����洢ÿ����õķ���ϵ�������д洢
       answer(i,n)=ANSS(i);
   end
    x0=0:1:max;
    y0=ANSS(1)*x0.^n    ; %������õ�ϵ����ʼ�����������ʽ����
    for num=2:1:n+1     
        y0=y0+ANSS(num)*x0.^(n+1-num);
    end
    subplot(1,3,n)
    plot(B,C,'*')
    hold on
    plot(x0,y0)
end
suptitle('��ͬ��������������Ͻ������1��3��')
figure,imshow(img);
hold on
plot(y0,x0,'r.')

%���׶���ʽ
syms x
y=answer(1,3)*x^3+answer(2,3)*x^2+answer(3,3)*x+answer(4,3);
y=diff(y);
x=86;
k1=eval(vpa(subs(y)));
x=661;
k2=eval(vpa(subs(y)));

if k1 ~= inf && k2 ~= inf
    q = abs(atan((k1-k2)/(1+k1*k2))/pi*180);
elseif k1 == inf && k2 ~= inf
    q = abs((pi/2-atan(k2))/pi*180);
elseif k1 ~= inf && k2 == inf
    q = abs((pi/2+atan(k1))/pi*180);
elseif k1 == inf && k2 == inf
    q = 0;
end

if q > 90
    q = 180-q;
end



