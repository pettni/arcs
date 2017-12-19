function [W, Cw, cont] = pginv(ts, U, G, Z, B, quant1)
  % Compute U-controlled set
  % contained in G \cap B \setdiff Z that can be
  % used to force a transition to Z using the progress
  % group (U,G) under (quant1, forall)-controllability
  %  
  % Returns a sorted set
  % 
  % Contracting algorithm

  persistent test_count
  if isempty(test_count)
     test_count = 1;
  end
  persistent file_perm
  global run_setting
  
  W1 = setdiff(intersect(uint32(G), uint32(B)), uint32(Z));
  W = W1;

  if (test_count == 10)
     disp('here');
  end
  test_count = test_count + 1;
  if ts.b_disable_pg || ...  % pgs disabled
     (~quant1 && ~isempty(setdiff(1:ts.n_a, U))) || ... % uncontrolled actions
     isempty(intersect(ts.pre(Z, U, quant1, 1), W))  % no reach to Z
    W = [];
    Cw = [];
    cont = Controller(W, containers.Map(), 'simple');
    
    if (strcmp(run_setting, 'write') && isempty(intersect(ts.pre(Z, U, quant1, 1), W)))
        if isempty(file_perm)
            file_perm = 'w';
            test_count = 0;
        end
        
        disp('Writing to file');
        fileID = fopen('inv_test.txt',file_perm);
        fprintf(fileID, 'test\n');
        fprintf(fileID, 'U ');
        fprintf(fileID, '%d ', [length(U), U']);
        fprintf(fileID, '\n');
        fprintf(fileID, 'G ');
        fprintf(fileID, '%d ', [length(G), G']);
        fprintf(fileID, '\n');
        fprintf(fileID, 'Z ');
        fprintf(fileID, '%d ', [length(Z), Z]);
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
    return
  end

  while true
    % union is sorted
    preW = ts.pre(union(W, Z), U, quant1, false);
    Wt = intersect(W, preW);
    if length(W) == length(Wt)
      break
    end
    W = Wt;
  end

  if nargout > 1
    % Contracting: first set minus final set
    Cw = setdiff(W1, W);
  end

  if nargout > 2
    [~, cont] = ts.pre(union(W, Z), U, quant1, false);
    cont.restrict_to(W);
    cont.from = 'pginv';
  end
  
  if (strcmp(run_setting, 'write'))
      if isempty(file_perm)
          file_perm = 'w';
          test_count = 0;
      end

      disp('Writing to file');
      fileID = fopen('inv_test.txt',file_perm);
      fprintf(fileID, 'test\n');
      fprintf(fileID, 'U ');
      fprintf(fileID, '%d ', [length(U), U']);
      fprintf(fileID, '\n');
      fprintf(fileID, 'G ');
      fprintf(fileID, '%d ', [length(G), G']);
      fprintf(fileID, '\n');
      fprintf(fileID, 'Z ');
      fprintf(fileID, '%d ', [length(Z), Z]);
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
      fprintf(fileID, '%d ', [length(W), W']);
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
