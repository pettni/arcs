load poly

clf
part.plot
legend off
title('FTS without disturbance')

matlab2tikz('poly.tex','interpretTickLabelsAsTex',true, ...
         'width','\figurewidth', 'height', '\figureheight', ...
         'parseStrings',false, 'showInfo', false)

load poly_pg

clf
part.plot
title('AFTS without disturbance')

matlab2tikz('poly_pg.tex','interpretTickLabelsAsTex',true, ...
         'width','\figurewidth', 'height', '\figureheight', ...
         'parseStrings',false, 'showInfo', false)

load poly_dist

clf
part.plot
legend off
title('FTS with disturbance')

matlab2tikz('poly_dist.tex','interpretTickLabelsAsTex',true, ...
         'width','\figurewidth', 'height', '\figureheight', ...
         'parseStrings',false, 'showInfo', false)


load poly_pg_dist

clf
part.plot
legend off
title('AFTS with disturbance')

matlab2tikz('poly_pg_dist.tex','interpretTickLabelsAsTex',true, ...
         'width','\figurewidth', 'height', '\figureheight', ...
         'parseStrings',false, 'showInfo', false)

