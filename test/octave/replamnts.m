function Replac = replamnts(nb, ipp, S, cellsize, replmnts, S1)
  defects = ipp.tally.defects;
  Replac.Repl(1) = sum(S.RC*cellsize*10);% srim
  Replac.Repl(2) = sum(replmnts);% iradina
if (nb == 3 )|| (nb == 4) || (nb == 5),
  Replac.Repl(3) = sum(defects.Replacements(:,1)+defects.Replacements(:,2)+defects.Replacements(:,3)); %iradina++
else
  Replac.Repl(3) = sum(defects.Replacements(:,1)+defects.Replacements(:,2)); %iradina++
end

if ~isempty(S1)
    Replac.Repl(4) = sum(S1.RC*cellsize*10);  % srim monolayer
endif

str = 'repl/ion';
disp(sprintf("%s\t%-10.4f\t%-10.4f\t%-10.4f\t%-10.4f\t%-10.4f",str,Replac.Repl))
end
