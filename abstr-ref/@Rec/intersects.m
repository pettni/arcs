function ret = intersects(rec1, rec2)
  % return true if rec1 and rec2 intersect
  if length(rec1)>1
    ret = false;
    for rec = rec1
      if intersects(rec, rec2)
        ret = true;
        return;
      end
    end
    return;
  end
  if length(rec2)>1
    ret = intersects(rec2, rec1);
    return;
  end
  ret = all(rec1.xmax>=rec2.xmin) || all(rec1.xmin<=rec2.xmax);
end