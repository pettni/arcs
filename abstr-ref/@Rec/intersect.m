function irec = intersect(rec1, rec2)
  % Computes the intersection between two Rec's rec1 and rec2

  if length(rec1)>1
    irec = [];
    for rec = rec1
      isect = intersect(rec, rec2)
      if ~isect.isEmptySet
        irec = [irec isect];
      end
    end
    return;
  end
  if length(rec2)>1
    irec = intersect(rec2, rec1);
    return;
  end

  imin = max(rec1.xmin, rec2.xmin);
  imax = min(rec1.xmax, rec2.xmax);

  % if ~isempty(setxor(rec1.ap, rec2.ap)) && (~isempty(rec1.ap) || ~isempty(rec2.ap))
  %   warning('intersecting sets with different APs. Keeping both..')
  % end
  ap = union(rec1.ap, rec2.ap);
  irec = Rec([imin; imax], ap);
end