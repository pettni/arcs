function calG = compute_progress_group(act_set, part)

    calG = cell(1,length(act_set));
    for act_ind = 1:length(act_set)
        calG{act_ind} = {};
        if isTransientLin(part.domain,act_set{act_ind})
            calG{act_ind} = [calG{act_ind}, [1:length(part)]];
        else
            for st_ind = 1:length(part)
                if isTransientLin(part(st_ind), act_set{act_ind})
                    calG{act_ind} = [calG{act_ind}, [st_ind]];
                end
            end
        end
    end
end