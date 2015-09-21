1;
A=trajcomp_load("../samples/prague1.dat");
full = trajcomp_get(A);
B=trajcomp_douglas_peucker(A,25);
simple=trajcomp_get(B);
clf;
hold on;
plot(simple(:,1),simple(:,2),"r");
plot(full(:,1),full(:,2),"b");
