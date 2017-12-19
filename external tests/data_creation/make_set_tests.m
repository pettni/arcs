function make_set_tests(data_file, setting)
    clear pginv
    clear pre_pg
    clear pre
    clear win_until
    clear win_until_and_always
    clear win_intermediate
    clear win_primal
    
    global run_setting
    run_setting = setting;
    load(data_file);
    part = s.part;
    
    Win = [];
    B = part.get_cells_with_ap({'SET'});
    if (strcmp(run_setting, 'time'))
        time = cputime;
    end
    [Win, Cwin] = part.ts.win_primal([], B, [], 'exists', 'forall', Win);
    if (strcmp(run_setting, 'time'))
        time = cputime - time;
        disp(['win_primal (win set and cand set) runtime : ', num2str(time), ' s']);
    end
end