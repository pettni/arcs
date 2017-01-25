function plot(part, alpha, text)
  % Plot a partition, regions with the same AP are plotted
  % in the same color.

  if nargin<2
    alpha = 1;
  end
  if nargin<3
    text = 0;
  end
  hh = ishold; if (~hh) clf; end
  hold on;
  if ~isempty(part.get_cells_with_ap([]))
    plot(part.cell_list(part.get_cells_with_ap([])), -1, alpha, 1);
  end
  aps = part.get_all_aps;
  colors = winter(length(aps));
  h_list = [];
  legend_list = {};
  for i=1:length(aps)
    aplist = part.get_cells_with_ap(aps(i));
    h = plot(part.cell_list(aplist), colors(i,:), alpha, 0);
    h_list = [h_list h];
    legend_list{end+1} = ['AP ' int2str(aps(i))];
  end
  if text
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
  legend(h_list,legend_list)
  if (~hh) hold off; end
end