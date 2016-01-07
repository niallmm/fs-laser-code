% check stuff
% urandom = [1*rand(Z+K1+2,1)];
% urandom2 =[1*rand(Z+K1+2,1)];
delta = 1e-12;


Rand1 = [4*rand(K1,1)];
state = ones(1,K1) + 1;
%Rand2 = [2*rand(K1,1)];


res      = getDTS(Rand1,state);
respertT = getDTS(Rand1+delta,state);
%respertC = getDTL(Rand1,Rand2+delta);
 
d_dTnum = (respertT-res)./delta;
%d_dCnum = (respertC-res)./delta;


d_dT = getDTSdT(Rand1,state);
%d_dC = getDTdCL(Rand1,Rand2);


max(abs(d_dT-d_dTnum))
%max(abs(d_dC-d_dCnum))


figure(33)

temp = [0.3:0.01:10]; % temp
plot(temp,getDTS(temp,ones(1,length(temp))))
hold on

plot(temp,getDTSdT(temp,ones(1,length(temp))))

%plot(temp,getDTdTL(temp,temp))

hold off

