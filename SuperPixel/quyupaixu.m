function [labelpic]=quyupaixu(labelpic)
[u1,v1] = find(labelpic==1);
[u2,v2] = find(labelpic==2);
[u3,v3] = find(labelpic==3);
[u4,v4] = find(labelpic==4);
[u5,v5] = find(labelpic==5);
[u6,v6] = find(labelpic==6);
m1 = v1(1,1);
m2 = v2(1,1);
m3 = v3(1,1);
m4 = v4(1,1);
m5 = v5(1,1);
m6 = v6(1,1);
n1 = u1(1,1);
n2 = u2(1,1);
n3 = u3(1,1);
n4 = u4(1,1);
n5 = u5(1,1);
n6 = u6(1,1);
A = [n1,n2,n3,n4,n5,n6];
V = sort(A);
f1 = length(v1);
f2 = length(v2);
f3 = length(v3);
f4 = length(v4);
f5 = length(v5);
f6 = length(v6);

[x,y] = find(V==n1);
for i = 1:f1
labelpic(u1(i,1),v1(i,1))=y;
end

[x,y] = find(V==n2);
for i = 1:f2
labelpic(u2(i,1),v2(i,1))=y;
end

[x,y] = find(V==n3);
for i = 1:f3
labelpic(u3(i,1),v3(i,1))=y;
end

[x,y] = find(V==n4);
for i = 1:f4
labelpic(u4(i,1),v4(i,1))=y;
end

[x,y] = find(V==n5);
for i = 1:f5
labelpic(u5(i,1),v5(i,1))=y;
end

[x,y] = find(V==n6);
for i = 1:f6
labelpic(u6(i,1),v6(i,1))=y;
end

end