function calG = compute_progress_group_nld(act_set, part,vars, deg, dBound)

    calG = cell(1,length(act_set));
    for act_ind = 1:length(act_set)
        calG{act_ind} = {};
        if isTransientNLinD(part.domain,act_set{act_ind},vars, deg, dBound)
            calG{act_ind} = [calG{act_ind}, [1:length(part)]];
        else
            for st_ind = 1:length(part)
                if isTransientNLinD(part(st_ind), act_set{act_ind},vars, deg, dBound)
                    calG{act_ind} = [calG{act_ind}, [st_ind]];
                end
            end
        end
    end
end