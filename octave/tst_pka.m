clear



folder = '../test/iradina++/bench_FeFC/';

A = load([folder 'FeFC.pka.h5']);

ionId = A.pka(1,:);
atomId = A.pka(2,:);
cellId = A.pka(3,:);
Er = A.pka(4,:);
Tdam = A.pka(5,:);
V = A.pka(6,:);
R = A.pka(7,:);
I = A.pka(8,:);

disp(["Nions = " num2str(ionId(end))])
disp(["Npkas = " num2str(length(ionId))])

Ed = 40;
L = 2.5*Ed;
Emax = 5*L;

i = find(Tdam < Emax);

figure 1
clf

subplot(2,2,1)
x = [Ed L Emax];
y = [1 1 Emax/L];
loglog(Tdam(i),V(i),'.',x,y)

subplot(2,2,2)
plot(Tdam(i),R(i),'.',x,y)

subplot(2,2,3)
x = [0 Emax];
plot(Er(i),Tdam(i),'.',x,x,x,[1 1]*Ed)

subplot(2,2,4)
semilogy(cellId(i),Er(i),'.')


