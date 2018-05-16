function [V, cont] = pre(ts, X, U, quant1, quant2)
  % Compute pre(X) under (quant1, quant2)-controllability
  % and action set U. If U == [] all actions are considered
  %
  % Note: X must be sorted!
  % Returns a sorted set
  %
  % quant = false: forall  (non controllable)
  % quant = true: exists   (controllable)

  if strcmp(ts.sys_setting, TransSyst.bdd_set)
    if (nargout == 1)
      V = ts.bdd_sys.pre(X, U, quant1, quant2);
    elseif nargout == 2
      [V, cont] = ts.bdd_sys.pre(X,U,quant1,quant2);
    end
    return;
  end
  
  if nargout > 1
    Kmap = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
  end

  if isempty(U)
    U = uint32(1:ts.n_a);
  end
  
  if(isempty(ts.trans_array))
      log_idx = false(1, ts.n_s);
      if ts.fast_enabled
        for i=1:length(X)
          for j = ts.fast_pre_all{X(i)}
            log_idx(j) = 1;
          end
        end
      else
        for i=1:ts.num_trans()
          if ismember(ts.state2(i), X) && ismember(ts.action(i), U)
            log_idx(ts.state1(i)) = 1;
          end
        end
      end

      for q = 1:ts.n_s
        if ~log_idx(q)
            continue
        end
        act_list = false(1, length(U));   % outcome per action 
        for i = 1:length(U)
          a = U(i);
          if ts.fast_enabled
            aPost = ts.fast_post{(a-1) * ts.n_s + q};
          else
            aPost = ts.post(q, a);
          end
          if quant2
            act_list(i) = any(builtin('_ismemberhelper',aPost, X));
          else
            act_list(i) = all(builtin('_ismemberhelper',aPost, X)) && ~isempty(aPost);
          end
        end
        if quant1
          if ~any(act_list)
            log_idx(q) = 0;
          elseif nargout > 1
            Kmap(q) = U(act_list); 
          end
        else
          if ~all(act_list)
            log_idx(q) = 0;
          elseif nargout > 1
            Kmap(q) = U;
          end
        end
      end
      
      V = zeros(1, sum(log_idx), 'uint32');
      V(:) = find(log_idx)';
  else
      if(ischar(quant1) || ischar(quant2))
          if(strcmp(quant1,'exists'))
            quant1 = true;
          elseif(strcmp(quant1,'forall'))
              quant1 = false;
          end

          if(strcmp(quant2,'exists'))
              quant2 = true;
          elseif(strcmp(quant2,'forall'))
              quant2 = false;
          end
      end
      
      if(quant1 == true && quant2 == false)
           V_list = false(ts.n_s,ts.n_a);
           for i = 1:length(U)
               u = U(i);
               sum1 = sum(ts.trans_array{u}(:,X),2);
               sum2 = sum(ts.trans_array{u},2);
               idx_V = sum1 == sum2 & sum1~=0;
               V_list(idx_V,u)=true;
           end
           V = uint32(find(sum(V_list,2)>=1))';
           if(nargout > 1)
              for i = 1:length(V) 
                  Kmap(V(i)) = uint32(find(V_list(V(i),:)~=0));
              end
           end
      elseif(quant1 == false && quant2 == true)
           idx_V = true(ts.n_s,1);
           for i = 1:length(U)
               u = U(i);
               sum1 = sum(ts.trans_array{u}(:,X),2);
               idx_V = idx_V & sum1~=0;
           end
           V = uint32(find(idx_V~=0))';
      elseif(quant1 == true && quant2 == true)
          V_list = false(ts.n_s,ts.n_a);
          for i = 1:length(U)
               u = U(i);
               sum1 = sum(ts.trans_array{u}(:,X),2);
               V_list(:,u) = sum1~=0;
          end
          V = uint32(find(sum(V_list,2)>=1))';
          if(nargout > 1)
              for i = 1:length(V) 
                  Kmap(V(i)) = uint32(find(V_list(V(i),:)~=0));
              end
          end
      elseif(quant1 == false && quant2 == false)
          idx_V = true(ts.n_s,1);
           for i = 1:length(U)
               u = U(i);
               sum1 = sum(ts.trans_array{u}(:,X),2);
               sum2 = sum(ts.trans_array{u},2);
               idx_V = idx_V & sum1==sum2 & sum1~=0;
           end
           V = uint32(find(idx_V~=0))';
      end
  end
   
  if nargout > 1
    cont = Controller(V, Kmap, 'simple', 'pre');
  end
end
