function plot(part, alpha, add_numbers)
  % plot(part, alpha, add_numbers): 
  % plot a 2d partition, regions with the same AP are plotted in the same color.

  if nargin<2
    alpha = 1;
  end
  if nargin<3
    add_numbers = 0;
  end

  hh = ishold; 
  if (~hh) 
    clf; 
  end
  hold on;

  aps = part.get_all_aps;
  colors = winter(length(aps));
  h_list = [];
  legend_list = {};

  recs_noap = part.get_cells_with_ap([]);
  if ~isempty(recs_noap)
    h = plot(part.cell_list(recs_noap), [0.5 0.5 0.5], alpha, 0);
  end

  % Sort aps by number of cells to plot "smallest" aps on top
  count_list = zeros(1,length(aps));
  for i=1:length(aps)
    count_list(i) = length(part.get_cells_with_ap(aps{i}));
  end

  [~, I] = sort(count_list, 2, 'descend');

  for i=I
    aplist = part.get_cells_with_ap(aps{i});
    h = plot(part.cell_list(aplist), colors(i,:), alpha, 0);
    h_list = [h_list h];
    legend_list{end+1} = [aps{i}];
  end
  if length(h_list) > 0
    legend(h_list,legend_list, 'location', 'NorthEastOutside')
  end

  if add_numbers
    for i = 1:length(part)
      xc = part.cell_list(i).getMidpoint();
      if part.dim == 2
        text(xc(1), xc(2), int2str(i), 'color', 'black');
      end
      if part.dim == 3
        text(xc(1), xc(2), xc(3), int2str(i), 'color', 'black');
      end
    end
  end 
  
  if (~hh) hold off; end
end