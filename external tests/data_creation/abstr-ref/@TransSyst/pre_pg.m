function [W, Cw, cont] = pre_pg(ts, V, B, quant1)
  % pre_pg: pre(V) under (quant1, forall) while remaining in B using progress groups
  % 
  % Returns a sorted set
  W = uint32(V);
  Vlist = {V};
  Cw = [];
  Klist = {Controller(W, containers.Map(), 'simple')};

  persistent test_count

  if (test_count == 4)
      disp('here');
  end
  for i=1:length(ts.pg_U)
    % Progress groups
    if nargout > 2
      [preVinv, Cw, preKinv] = ts.pginv(ts.pg_U{i}, ts.pg_G{i}, W, B, quant1);
      if length(preVinv) > 0
        Vlist{end+1} = preVinv;
        Klist{end+1} = preKinv;
      end
    elseif nargout > 1
      [preVinv, Cw] = ts.pginv(ts.pg_U{i}, ts.pg_G{i}, W, B, quant1);
    else
      preVinv = ts.pginv(ts.pg_U{i}, ts.pg_G{i}, W, B, quant1);
    end
    W = union(W, preVinv);
    W = reshape(W, 1, length(W));
  end

  if nargout > 2
    cont = Controller(Vlist, Klist, 'reach', 'pre_pg');
  end
  
  persistent file_perm
  global run_setting
  if isempty(test_count)
    test_count = 1;
  end
  test_count = test_count + 1;
  if (strcmp(run_setting, 'write'))
      if isempty(file_perm)
          file_perm = 'w';
          test_count = 1;
      end

      
      fileID = fopen('pre_pg_test.txt',file_perm);
      disp('Writing to file');
      fprintf(fileID, 'test\n');
      fprintf(fileID, 'V ');
      fprintf(fileID, '%d ', [length(V), V]);
      fprintf(fileID, '\n');
      fprintf(fileID, 'B ');
      fprintf(fileID, '%d ', [length(B), B]);
      fprintf(fileID, '\n');
      if quant1
          fprintf(fileID, 'q e\n');
      else
          fprintf(fileID, 'q a\n');
      end
      fprintf(fileID, 'ans ');
      fprintf(fileID, '%d ', [length(W), W]);
      fprintf(fileID, '\n');
      fprintf(fileID, 'mode %d\n', nargout);
      if nargout > 1
          fprintf(fileID, 'cand ');
          fprintf(fileID, '%d ', [length(Cw), Cw]);
          fprintf(fileID, '\n');
      end
      fclose(fileID);
      file_perm = 'a';
      disp('closed file');
  end
end