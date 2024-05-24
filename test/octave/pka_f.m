function pka = pka_f(ipp, S, cellsize, E0,pkai, S1),
  str = ' # PKAs ';
  pka(1) = sum(S.Vi)*cellsize*10; % srim pka/ion
  pka(2) = sum(pkai);
  pka(3) = sum(sum(ipp.tally.defects.PKAs(:,2:end)));
  if ~isempty(S1)
    pka(4) = sum(S1.Vi)*cellsize*10;  % srim monolayer
  endif
  disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e",str,pka))
end

