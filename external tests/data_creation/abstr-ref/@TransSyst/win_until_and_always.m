function [V, Cv, cont] = win_until_and_always(ts, A, B, P, quant1)
  % Compute the winning set of
  %   []A && B U P
  % under (quant1, forall)-controllability
  %
  % Note: A must be sorted
  % Returns a sorted set
  
    persistent test_count
    
    if (test_count == 12)
        disp('hello');
    end
  
  if ~isempty(setdiff(1:ts.n_s, A))
    % Need to worry about []A
    if nargout > 2
      [Vinv, Cvinv, cont_inv] = ts.win_intermediate(uint32(1:ts.n_s), A, [], {uint32(1:ts.n_s)}, quant1);
      [V, Cvu, cont] = ts.win_until(intersect(B, Vinv), intersect(P, Vinv), quant1);
      cont.from = 'win_until_and_always';
      Cv = union(Cvinv, Cvu);
    elseif nargout > 1
      [Vinv, Cvinv] = ts.win_intermediate(uint32(1:ts.n_s), A, [], {uint32(1:ts.n_s)}, quant1);
      [V, Cvu] = ts.win_until(intersect(B, Vinv), intersect(P, Vinv), quant1);
      Cv = union(Cvinv, Cvu);
    else
      Vinv = ts.win_intermediate(uint32(1:ts.n_s), A, [], {uint32(1:ts.n_s)}, quant1);
      V = ts.win_until(intersect(B, Vinv), intersect(P, Vinv), quant1);
    end
  else
    % No need to worry about []A
    if nargout > 2
      [V, Cv, cont] = ts.win_until(B, P, quant1);
      cont.from = 'win_until';
    elseif nargout > 1
      [V, Cv] = ts.win_until(B, P, quant1);
    else
      V = ts.win_until(B, P, quant1);
    end
  end
  
  persistent file_perm
  global run_setting
  if (strcmp(run_setting, 'write'))
      if isempty(file_perm)
          file_perm = 'w';
          test_count = 1;
      end
      
      test_count = test_count + 1;
      fileID = fopen('win_until_and_always_test.txt', file_perm);
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
