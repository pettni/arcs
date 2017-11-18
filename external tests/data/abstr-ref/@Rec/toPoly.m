function poly = toPoly(rec)
  % Convert Rec rec to Polyhedron object
  poly = Polyhedron('A', [eye(rec.dim); -eye(rec.dim)], 'b', [rec.xmax'; -rec.xmin' ]);
end