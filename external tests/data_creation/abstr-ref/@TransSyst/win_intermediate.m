function [V, Cv, cont] = win_intermediate(ts, A, B, P, C_list, quant1)
  % Compute winning set of
  %  []A && ( (B U P) || [] (B &&_i <>C_i) )
  % under (quant1, forall)-controllability
  %
  % Note: A must be sorted
  % Returns a sorted set
  % 
  % Contracting algorithm

  V = uint32(1:ts.n_s);
  Klist = {};

  iter = 1;
  while true
    Vt = V;
    preV = ts.pre(V, [], quant1, false);
    
    if nargout > 1
      Cv = [];
    end

    Qs = intersect(B, preV);
    for i=1:length(C_list)
      Qi = union(P, intersect(Qs, C_list{i}));
      Qi = reshape(Qi, 1, length(Qi));
      if nargout > 2 && iter == 1
        [Vti, Cvi, Ki] = ts.win_until_and_always(A, B, Qi, quant1);
        Klist{i} = Ki;
        Cv = union(Cv, Cvi);
      elseif nargout > 1 && iter == 1
        [Vti, Cvi] = ts.win_until_and_always(A, B, Qi, quant1);
        Cv = union(Cv, Cvi);
        disp('ran cand');
      else
        Vti = ts.win_until_and_always(A, B, Qi, quant1);
      end
      Vt = intersect(Vt, Vti);
    end
    
    if nargout > 1 && iter == 1
      V1 = Vt;
      C_rec = Cv;
    end

    if length(V) == length(Vt)
      break
    end
    V = Vt;
    iter = iter + 1;
  end

  if nargout > 1
    % Contracting: C_rec U (V_1\V_last)
    Cv = union(C_rec, setdiff(V1, V));
  end

  if nargout > 2
    Vlist = {V};
    for i = 1:length(C_list)
      Vlist{end+1} = intersect(B, C_list{i});
    end
    cont = Controller(Vlist, Klist, 'recurrence', 'win_intermediate');
  end
  
  persistent file_perm
  persistent test_count
  global run_setting
  if (strcmp(run_setting, 'write'))
      if isempty(file_perm)
          file_perm = 'w';
          test_count = 1;
      end
      
      test_count = test_count + 1;
      fileID = fopen('win_intermediate_test.txt', file_perm);
      disp('Writing to file');
      fprintf(fileID, 'test\n');
      fprintf(fileID, 'A ');
      fprintf(fileID, '%d ', [length(A), A]);
      fprintf(fileID, '\n');
      fprintf(fileID, 'B ');
      fprintf(fileID, '%d ', [length(B), B]);
      fprintf(fileID, '\n');
      fprintf(fileID, 'Z ');
      fprintf(fileID, '%d ', [length(P), P]);
      fprintf(fileID, '\n');
      fprintf(fileID, 'C %d', length(C_list));
      fprintf(fileID, '\n');
      for i = 1:length(C_list)
          fprintf(fileID, '%d ', [length(C_list{i}), C_list{i}]);
          fprintf(fileID, '\n');
      end
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
          if (size(Cv, 1) > 1)
              Cv = Cv';
          end
          fprintf(fileID, '%d ', [length(Cv), Cv]);
          fprintf(fileID, '\n');
      end
      fclose(fileID);
      file_perm = 'a';
      disp('closed file');
      
  end
end
