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

  if nargout > 1
    Kmap = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
  end

  if isempty(U)
    U = uint32(1:ts.n_a);
  end

  ts.trans_array_enable();
  
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

  if nargout > 1
    cont = Controller(V, Kmap, 'simple', 'pre');
  end

end
