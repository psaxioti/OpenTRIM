clear

nb =1

qc = srim2mat(['../test/srim/qc/' num2str(nb)]);
fc = srim2mat(['../test/srim/fc2/' num2str(nb)]);

dx = fc.x(2)-fc.x(1);

function y = NRT(T,Ed)
  y = T/Ed;
  i = find(y<1); y(i)=0;
  i = find(y>=1 & y<2.5); y(i) = 1;
  i = find(y>=2.5); y(i) = y(i)/2.5;
end

Eb = 3; Ed=40;
disp( "Vac/Ion   QC\tFC\tFC-NRT")
disp(sprintf("SRIM      %.0f \t%.0f \t%.0f",...
    sum(qc.Vr+qc.Vi)*dx,...
    sum(sum(fc.Vr))*dx,...
    sum(NRT(fc.Vr*Eb + fc.EPr,Ed))*dx))

ipp = load(['../test/iradina++/b' num2str(nb) '/b' num2str(nb) '.h5']);
N =ipp.Nions;
Ri = (ipp.Implantations(:,1)+ipp.Replacements(:,1))/N/dx/1e-8;
V = sum(ipp.Vacancies(:,2:end),2)/N/dx;
disp(sprintf("IRADINA++ %.0f \t%.0f \t%.0f",...
    sum(sum(ipp.Vnrt_LSS(:,2:end),2))/N,...
    sum(V)*dx, ...
    sum(sum(ipp.Vnrt(:,2:end),2))/N))

x = ipp.cell_xyz(:,1);

figure 1
clf

subplot(2,1,1)
plot(x,qc.Vi+qc.Vr,...
     x,sum(fc.Vr,2),...
     x,sum(ipp.Vnrt_LSS(:,2:end),2)/N/dx,'--',...
     x,sum(ipp.Vacancies(:,2:end),2)/N/dx,'--', ...
     x,sum(ipp.Vnrt(:,2:end),2)/N/dx,'--')
title('Vacancy production')
xlabel('nm')
ylabel('Vacancies/ion-A')
legend('SRIM QC','SRIM FC',...
       'IRADINA++ QC','IRADINA++ FC','IRADINA++ FC-NRT')


subplot(2,1,2)
Eb = 3;
plot(x,qc.Vi*Eb + qc.Vr*Eb + qc.EPr, ... # YR Lin D3***
     x,fc.Vr*Eb + fc.EPr, ...  # YR Lin D3**
     x,sum(ipp.Tdam_LSS(:,2:end),2)/N/dx,'--',...
     x,sum(ipp.PhononEnergy(:,2:end),2)/N/dx,'--', ...
     x,sum(ipp.Tdam(:,2:end),2)/N/dx,'--')
title('Tdam')
xlabel('nm')
ylabel('eV/ion-A')
legend('SRIM QC','SRIM FC',...
       'IRADINA++ QC','IRADINA++ FC','IRADINA++ FC-NRT')

figure 2
clf
x = ipp.cell_xyz(:,1);
plot(x,ipp.Tdam(:,2),x,ipp.PhononEnergy(:,2))
legend('Tdam','Fe')
disp(num2str([sum(ipp.Tdam(:,2)) sum(ipp.PhononEnergy(:,2))]/N))
