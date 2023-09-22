function [r3,c3,r4,c4,k] = zhaoqiexiandian(rectx,recty)

x1 = rectx(1,1); y1 = recty(1,1);
x2 = rectx(2,1); y2 = recty(2,1);
x3 = rectx(3,1); y3 = recty(3,1);
x4 = rectx(4,1); y4 = recty(4,1);

a1 = y1-y3; b1 =x3-x1; c1 = x1*y3-x3*y1;
a2 = y2-y4; b2 =x4-x2; c2 = x2*y4-x4*y2;

x = sym('x');
y = sym('y');
[x,y] = solve(a1*x+b1*y+c1==0,a2*x+b2*y+c2==0,x,y);

if x1 > x && y1 < y
    x = x1;y = y1;
elseif x2 > x && y2 < y
    x = x2;y = y2;
elseif x3 > x && y3 < y
    x = x3;y = y3;
elseif x4 > x && y4 < y
    x = x4;y = y4;
end
r1 = x;c1 = y;

if c1 == recty(1,1) && c1 ~= recty(2,1)
    x2 = rectx(2,1); y2 = recty(2,1);
elseif c1 == recty(2,1) && c1 ~= recty(3,1)
    x2 = rectx(3,1); y2 = recty(3,1);
elseif c1 == recty(3,1) && c1 ~= recty(4,1)
    x2 = rectx(4,1); y2 = recty(4,1);
elseif c1 == recty(4,1) && c1 ~= recty(5,1)
    x2 = rectx(5,1); y2 = recty(5,1);
end
r2 = x2; c2 = y2;

if r1 == r2
    n = 4*abs(c2-c1);
    c3 = c1-n;
    c4 = c2+n;
    r3 = r1;
    r4 = r2;
    k = inf;
elseif r1 ~= r2
    m = 2*abs(r2-r1);
    n = 4*abs(c2-c1);
    k = (c2-c1)/(r2-r1);
    b = c2-k*r2;
    %y变化
    c3 = c1-n;
    r3 = (c3-b)/k;
    c4 = c2+n;
    r4 = (c4-b)/k;

    %x变化
    % r3 = r1-m;
    % c3 = k*r3+b;
    % r4 = r2+m;
    % c4 = k*r4+b;
end

end