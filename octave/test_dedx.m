clear

Z1 = 26;
Z2 = 10:20:90;

[e, d] = dedx(Z1,Z2(1));

D = zeros(length(Z2),length(d));

for i=1:length(Z2)
  [e, D(i,:)] = dedx(Z1,Z2(i));
end

figure 1
clf
loglog(e,D)
