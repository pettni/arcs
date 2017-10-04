
number_of_systems = 7;
encs = {'log', 'split'};
results = zeros(number_of_systems, length(encs));
for s = 1:number_of_systems
  for enc_ind = 1:length(encs) 
    data_string = ['system_data/BDD_statetest', num2str(s), '.mat'];
    system_string = ['system_data/BDD_statetest', num2str(s), '_', encs{enc_ind}, '.txt'];
    part_data = load(data_string);
    fileID = fopen(system_string, 'r');
    
    % transition system
    state_args = [];
    fscanf(fileID, '%s', 1)
    state_args = fscanf(fileID, '%d', 2)';
    
    state_encs = {};
    for i = 1:(state_args(1)+1)
      id = fscanf(fileID, '%d', 1);
      enc_string = fscanf(fileID, '%s', 1);
      enc = zeros(1, length(enc_string));
      for j = 1:length(enc_string)
        enc(j) = str2num(enc_string(j));
      end
      state_encs{id} = enc;
    end
    
    trans_args = [];
    fscanf(fileID, '%s', 1)
    trans_args = fscanf(fileID, '%d', 2)';
    
    system = TransSyst(state_args(1)+1, trans_args(1), TransSyst.bdd_set, 'split');
    
    sparse_system = TransSyst(state_args(1)+1, trans_args(1), 'sparse');
    for i = 1:length(state_encs)
      system.bdd_sys.set_state_enc(i, uint32(state_encs{i}));
    end
    
    for i = 1:trans_args(2)
      state1 = fscanf(fileID, '%d', 1);
      a_enc_string = fscanf(fileID, '%s', 1);
      a_enc = zeros(1, length(a_enc_string));
      state2 = fscanf(fileID, '%d', 1);
      for j = 1:length(a_enc)
        a_enc(j) = str2num(a_enc_string(j));
      end
      action = bi2de(a_enc) + 1;
      
      system.add_transition(state1, state2, action);
      sparse_system.add_transition(state1, state2, action);
    end
    sparse_system.create_fast();
    
    % read progress groups
    fscanf(fileID, '%s', 1);
    pg_num = fscanf(fileID, '%d', 1);
    for i = 1:pg_num
      fscanf(fileID, '%s', 1);
      group_args = fscanf(fileID, '%d', 2);
      actions = fscanf(fileID, '%d', group_args(1));
      states = fscanf(fileID, '%d', group_args(2));
      system.add_progress_group(uint32(actions), double(states));
      sparse_system.add_progress_group(actions, states);
    end
    fclose(fileID);
    
    % win_primal
    
    B = part_data.part.get_cells_with_ap({'SET'});
    
    quant1 = 'exists';
    quant2 = 'forall';
    
    system.bdd_sys.reorder(5000);
    
    time = system.bdd_sys.win_primal_time([], B, [], quant1, quant2);
    results(s, enc_ind) = time;
  end

end