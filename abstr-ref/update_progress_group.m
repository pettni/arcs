function calG = update_progress_group(calG, part, chInd)
    for act_ind = 1:length(calG)
        G = calG{act_ind};
        for gind = 1:length(G)
            if not(isempty(find(G{gind}== chInd,1)))
                G{gind} = [G{gind}, length(part)];
            end
        end
        calG{act_ind} = G;
    end
end