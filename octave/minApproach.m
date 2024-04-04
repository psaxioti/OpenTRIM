function xm = minApproach(e,s)

ffunc = @(x) 1 - screening_zbl(x)./x/e - s*s./x.^2;

x0 = 1/2./e + sqrt(1/4./e.^2 + s.*s);
x2 = x0;
x1 = x2/10;

while ffunc(x2)<=0, x2 = x2*1.001 end
while ffunc(x1)>=0, x1 = x1*0.1 end

xm = 0.5*(x1+x2);
fm = ffunc(xm);
d = 1;
while abs(d) > eps && abs(fm) > eps,
  if fm<0, x1 = xm; else x2 = xm; end
  d = (0.5*(x1+x2)-xm)/xm;
  xm = 0.5*(x1+x2);
  fm = ffunc(xm);
end
