% check stuff
% urandom = [1*rand(Z+K1+2,1)];
% urandom2 =[1*rand(Z+K1+2,1)];
delta = 1e-8;


Rand1 = [2*rand(K1,1)];
Rand2 = [2*rand(K1,1)];


res      = getDCL(Rand1,Rand2);
respertT = getDCL(Rand1+delta,Rand2);
respertC = getDCL(Rand1,Rand2+delta);
 
d_dTnum = (respertT-res)./delta;
d_dCnum = (respertC-res)./delta;


d_dT = getDCdTL(Rand1,Rand2);
d_dC = getDCdCL(Rand1,Rand2);


max(abs(d_dT-d_dTnum))
max(abs(d_dC-d_dCnum))


figure(33)

temp = [0:0.01:10]; % temp
plot(temp,getDCL(temp,temp))
hold on

plot(temp,getDCdTL(temp,temp))

hold off

