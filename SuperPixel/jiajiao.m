function [q]=jiajiao(K,i)
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

k1 = K(1,i);
k2 = K(1,i+1);
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
end

