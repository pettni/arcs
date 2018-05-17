function ret = post(ts, q, a)
  % post(ts, q, a): compute the post set of state q and action a
  % 
  % Returns a sorted set
  
  if strcmp(ts.sys_setting, TransSyst.bdd_set)
      ret = ts.post(ts, q, a);
      return;
  end
  
  ret = uint32([]);
  for i = 1:ts.num_trans()
    if q == ts.state1(i) && a == ts.action(i)
      ret(end+1) = ts.state2(i);
    end
  end
  ret = unique(ret);

end
