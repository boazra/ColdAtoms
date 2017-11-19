n = 10000;
v1 = rand(n,1);
v2 = rand(n,1);
v3 = rand(n,1);
V = [v1,v2,v3];
a1 = rand(1)*1000;
a2 = rand(1)*1000;
a3 = rand(1)*1000;
b = a1*v1 + a2*v2 + a3*v3;
x = linsolve(V,b);
[a1,a2,a3]'-x

