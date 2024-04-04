function xm = findS(e,thetaCM)

# inital guesses: Mendenhall & Weller NIMB 58 (1991) 11, eqs. 23-25
d = thetaCM/pi;
gam = 1-d;
if d<10*eps,
  e0 = 2*d*e;
else
  e0 = (1-gam*gam)*e;
end

x0 = minApproach(e0, 1.e-8)

if e >= 1.,
  x1 = 0.7*gam*x0;
  x2 = 1.0/(2.0*e*tan(thetaCM/2.0));
else
  x1 = 0.9*gam*x0;
  x2 = 1.4*gam*x0;
end

if x2>1e4, x2=1e4 end

dtheta = @(x) thetaCM - screened_coulomb_theta(e,x);

while dtheta(x1)>=0.0, x1 *= 0.1; end
while dtheta(x2)<=0.0, x2 *= 1.001; end

if dtheta(x1)>=0.0 || dtheta(x2)<=0.0,
  error("no bracketing!");
end


xm = 0.5*(x1+x2);
fm = dtheta(xm);
k = 0;
d = 1;
while abs(d)>eps,
  if fm<0, x1 = xm; else x2 = xm; end
  q = 0.5*(x1+x2);
  d = (q - xm)/xm;
  xm = 0.5*(x1+x2);
  fm = dtheta(xm);
  k++;
end
