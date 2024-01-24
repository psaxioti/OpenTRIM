clear



folder = '/home/george/Dev/ionsim/iradinapp/test/iradina++/bench_FeFC/';

A = load([folder 'FeFC.h5']);
x = A.X;
x(end) = [];
x += (x(2)-x(1))/2;
N = double(A.Nions);
disp(["Nions = " num2str(A.Nions)])
disp(["Npkas = " num2str(A.Npkas)])

figure 1

subplot(2,2,1)
plot(x, double(A.Vacancies)/N);
title('Vacancies')

subplot(2,2,2)
plot(x, double(A.Implantations)/N);
title('Implanted Ions')

subplot(2,2,3)
plot(x, A.IonizationEnergy/N);
title('Ionization Energy')

subplot(2,2,4)
plot(x, A.PhononEnergy/N);
title('Phonon Energy')



return

sum(A.PhononEnergy)

sum(A.IonizationEnergy)

sum(sum(A.PhononEnergy))

sum(sum(A.IonizationEnergy))

sum(sum(A.IonizationEnergy+A.PhononEnergy))/N

v = sum(A.Vacancies(1:50,:))/N

sqrt(sum(A.Vacancies(1:50,:)))/N

sum(v)/v(1)

sum(A.PhononEnergy(1:50,:))/N/v(1)

V1 = double(sum(A.Vacancies,2))./double(A.Vacancies(:,1));

Tdam = A.PhononEnergy(:,2)./A.Vacancies(:,1);

i=1:100;
plot(i,Tdam(i,:))
