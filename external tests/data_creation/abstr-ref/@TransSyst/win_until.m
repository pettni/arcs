function [V, Cv, cont] = win_until(ts, B, P, quant1)
  % Compute the winning set of
  %  B U P
  % under (quant1, forall)-controllability
  %
  % Returns a sorted set
  %
  % Exanding algo
  
  persistent test_count

  V = uint32([]);
  Vlist = {};
  Klist = {};
  while true
    if nargout > 2
      % Normal pre
      [preV, preK] = ts.pre(V, [], quant1, false);
      Vlist{end+1} = preV;
      Klist{end+1} = preK;
      Vt = union(P, intersect(B, preV));

      % PG pre
      [preVinv, CpreVinv, preKinv] = ts.pre_pg(Vt, B, quant1);
      if ~isempty(setdiff(preVinv, Vt))
        Vlist{end+1} = preVinv;
        Klist{end+1} = preKinv;
      end
      Vt = union(Vt, preVinv);
    elseif nargout > 1
      % Normal pre
      preV = ts.pre(V, [], quant1, false);
      Vt = union(P, intersect(B, preV));

      % PG pre
      [preVinv, CpreVinv] = ts.pre_pg(Vt, B, quant1);
      Vt = union(Vt, preVinv);
    else
      preV = ts.pre(V, [], quant1, false);
      Vt = union(P, intersect(B, preV));
      preVinv = ts.pre_pg(V, B, quant1);
      Vt = union(Vt, preVinv);
    end

    Vt = reshape(Vt, 1, length(Vt));

    if length(V) == length(Vt)
      break
    end

    V = Vt;
  end

  if nargout > 1
    Cv = union(CpreVinv, setdiff(ts.pre(V, [], quant1, false), V));
  end

  if nargout > 2
    cont = Controller(Vlist, Klist, 'reach', 'win_until');
  end
  
  persistent file_perm
  global run_setting
  if (strcmp(run_setting, 'write'))
      if isempty(file_perm)
          file_perm = 'w';
          test_count = 1;
      end
      
      test_count = test_count + 1;
      fileID = fopen('win_until_test.txt', file_perm);
      disp('Writing to file');
      fprintf(fileID, 'test\n');
      fprintf(fileID, 'B ');
      fprintf(fileID, '%d ', [length(B), B]);
      fprintf(fileID, '\n');
      fprintf(fileID, 'Z ');
      fprintf(fileID, '%d ', [length(P), P]);
      fprintf(fileID, '\n');
      if quant1
          fprintf(fileID, 'q e\n');
      else
          fprintf(fileID, 'q a\n');
      end
      fprintf(fileID, 'ans ');
      fprintf(fileID, '%d ', [length(V), V]);
      fprintf(fileID, '\n');
      fprintf(fileID, 'mode %d\n', nargout);
      if nargout > 1
          fprintf(fileID, 'cand ');
          fprintf(fileID, '%d ', [length(Cv), Cv']);
          fprintf(fileID, '\n');
      end
      fclose(fileID);
      file_perm = 'a';
      disp('closed file');
  end
end
