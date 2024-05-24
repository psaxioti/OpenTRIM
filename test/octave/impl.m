function X = impl(S,ions_total, ipp, cellsize, S1)
  X(1) = sum(S.Ri.*cellsize*10e-8);  % srim
  X(2) = sum(ions_total); % iradina
  X(3) = sum(ipp.tally.defects.Implantations(:,1)+ipp.tally.defects.Replacements(:,1)); %iradina++
  if ~isempty(S1)
    X(4) = sum(S1.Ri.*cellsize*10e-8);  % srim monolayer
  endif
str = 'Impl. Ions';
disp(sprintf("%s\t%-10.4e\t%-10.4e\t%-10.4e\t%-10.4e",str,X))
end

