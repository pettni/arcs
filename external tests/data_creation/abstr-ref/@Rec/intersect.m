function irec = intersect(rec1, rec2)
  % Computes the intersection between two Rec's rec1 and rec2

  if length(rec1)>1
    irec = [];
    for rec = rec1
      isect = intersect(rec, rec2);
      if isect.isFullDim
        irec = [irec isect];
      end
    end
    return;
  end
  if length(rec2)>1
    irec = intersect(rec2, rec1);
    return;
  end

  irec = Rec([max(rec1.xmin, rec2.xmin); ...
              min(rec1.xmax, rec2.xmax)], ...
              union(rec1.ap, rec2.ap));
end