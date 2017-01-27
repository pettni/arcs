function ret = intersects(rec1, rec2)
  % return true if rec1 and rec2 intersect
  if length(rec1)>1 || length(rec2)>1
    ret = false;
    for i=1:length(rec1)
      for j=1:length(rec2)
        if intersects(rec1(i), rec2(j))
          ret = true;
          return;
        end
      end
    end
    return
  end

  ret = all(max(rec1.xmin, rec2.xmin) <= min(rec1.xmax, rec2.xmax));

end