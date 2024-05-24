function Tdam = dam_ener(ipp, S, cellsize, E0, S1)
Tdref = sum(sum(ipp.tally.damage.Tdam(:,2:end)));
% Tdam D1
str = 'Td/Tdam D1';
Tdam.Tdam1(1) = (sum(sum(S.ERr(1:end))) - sum(S.EIr))*cellsize*10;
Tdam.Tdam1(2) = NaN;
Tdam.Tdam1(3) = Tdref;
if ~isempty(S1)
    Tdam.Tdam1(4) = (sum(sum(S1.ERr(1:end))) - sum(S1.EIr))*cellsize*10;  % srim monolayer
endif
Tdam.difff1 = Tdam.Tdam1./Tdref;
disp(sprintf("%s\t%-10.5e\t%-10.5e\t%-10.5e\t%-10.5e",str,Tdam.difff1))

% Tdam D2
str = 'Td/Tdam D2';
Tdam.Tdam2(1) = E0 - sum(S.EIi + S.EIr + S.EPi)*cellsize*10; %srim
Tdam.Tdam2(2) = NaN; %iradina
Tdam.Tdam2(3) = sum(sum(ipp.tally.damage.Tdam(:,2:end))); %iradina++
if ~isempty(S1)
    Tdam.Tdam2(4) = E0 - sum(S1.EIi + S1.EIr + S1.EPi)*cellsize*10;  % srim monolayer
endif
Tdam.difff2 = Tdam.Tdam2./Tdref;
disp(sprintf("%s\t%-10.5e\t%-10.5e\t%-10.5e\t%-10.5e",str,Tdam.difff2))

% Tdam D3
str = 'Td/Tdam D2';
Tdam.Tdam3(1) = sum(S.EPr - S.EPi)*cellsize*10;
Tdam.Tdam3(2) = NaN;
Tdam.Tdam3(3) = Tdref;
if ~isempty(S1)
    Tdam.Tdam3(4) = sum(S1.EPr - S1.EPi)*cellsize*10;  % srim monolayer
endif
Tdam.difff3 = Tdam.Tdam3./Tdref;
disp(sprintf("%s\t%-10.5e\t%-10.5e\t%-10.5e\t%-10.5e",str,Tdam.difff3))

% Tdam D1*
str = 'Td/Tdam D1*';
Tdam.Tdam4(1) = (sum(sum(S.ERr(:,1:end))) - sum(S.EIr) + sum(S.EPi))*cellsize*10;
Tdam.Tdam4(2) = NaN;
Tdam.Tdam4(3) = Tdref;
if ~isempty(S1)
    Tdam.Tdam4(4) = (sum(sum(S1.ERr(:,1:end))) - sum(S1.EIr) + sum(S1.EPi))*cellsize*10;  % srim monolayer
endif
Tdam.difff4 = Tdam.Tdam4./Tdref;
disp(sprintf("%s\t%-10.5e\t%-10.5e\t%-10.5e\t%-10.5e",str,Tdam.difff4))

% Tdam D2**
str = 'Td/Tdam D2**';
Tdam.Tdam5(1) = E0 - sum(S.EIi + S.EIr)*cellsize*10;
Tdam.Tdam5(2) = NaN;
Tdam.Tdam5(3) = Tdref;
if ~isempty(S1)
    Tdam.Tdam5(4) = E0 - sum(S1.EIi + S1.EIr)*cellsize*10;  % srim monolayer
endif
Tdam.difff5 = Tdam.Tdam5./Tdref;
disp(sprintf("%s\t%-10.5e\t%-10.5e\t%-10.5e\t%-10.5e",str,Tdam.difff5))

smb = '^vs*d';

figure 6
clf
%plot sum Tdam
tst ={'D1', 'D2', 'D3', 'D1*', 'D2**'};
clr=get(gca,'colororder');
hold on
for i = 1:5,
   var_name = sprintf('Tdam.Tdam%d\t', i);
    % Access the variable using eval
    y = eval(var_name);
    h(i) =  plot(i,y(1)/1e3,['--o;' tst{i} ';'],...
                 'color',clr(i,:),'markerfacecolor',clr(i,:));
    if ~isempty(S1)
      h1(i) =  plot(i,y(4)/1e3,['--v;' tst{i} ';'],...
                 'color',clr(i,:),'markerfacecolor',clr(i,:));
    endif
end
set(gca,'xticklabel',{'D1','D2','D3','D1*','D2**'})
plot([1 5],[Tdref/1e3 Tdref/1e3],'--r')
hold off
box on
legend(h,'location','southeast','interpreter','latex')
xlabel('$T_{dam}$ Method','interpreter','latex')
ylabel(['$T_{dam}^{Srim} (keV) $'], ...
        'interpreter','latex')


figure 7
clf
%plot difference with Tdref
tst ={'D1', 'D2', 'D3', 'D1*', 'D2**'};
clr=get(gca,'colororder');
hold on
for i = 1:5,
   var_name = sprintf('Tdam.difff%d\t', i);
    % Access the variable using eval
    y = eval(var_name);
    h(i) =  plot(i,y(1),['--o;' tst{i} ';'],...
                 'color',clr(i,:),'markerfacecolor',clr(i,:));
    if ~isempty(S1)
      h1(i) =  plot(i,y(4),['--v;' tst{i} ';'],...
                 'color',clr(i,:),'markerfacecolor',clr(i,:));
    endif
end
set(gca,'xticklabel',{'D1','D2','D3','D1*','D2**'})
plot([1 5],[1 1],'--r')
hold off
box on
legend(h,'location','southeast','interpreter','latex')
xlabel('$T_{dam}$ Method','interpreter','latex')
ylabel(['$T_{dam}^{Srim} / T_{dam}^{irad++}$'], ...
        'interpreter','latex')


end
