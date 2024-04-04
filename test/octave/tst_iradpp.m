clear

folder = '../test/iradina++/bench_FeFC/';

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
plot(x, double(A.Implantations(:,1))/N);
title('Implanted ions')

subplot(2,2,3)
plot(x, A.IonizationEnergy/N);
title('Ionization Energy')

subplot(2,2,4)
plot(x, A.PhononEnergy/N);
title('Phonon Energy')

