Linfty = 10;
maxRelError = 1e-5;
a = 1.15;
b = 0.0169;
c = 0.068;

DT = 82.4;
Tmelt = 1;
mu = 36;
hdot = mu*(2*a-Tmelt);

t_tilde = 0.02/hdot;
U_tilde = hdot;

f = @(xint)a*(exp(-b*xint)+exp(-c*xint)).* ...
                getGreensFct(U_tilde.*t_tilde, xint, t_tilde, 0.0,DT);
   
            
spint = quadgk(f,0,Linfty,'RelTol',1e-9);

xint = linspace(0,30,1000);
integrand = f(xint);
figure
plot(xint, integrand)