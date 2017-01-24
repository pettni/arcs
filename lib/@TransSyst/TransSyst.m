classdef TransSyst<handle

	properties (SetAccess=protected)
		n_s;
		n_a;
		% Transitions
		state1;
		state2; 
		action;
		% Progress groups
		pg_U = {};
		pg_G = {};

		% Fastaccess_first
		fast_forward = containers.Map();
		fast_backward = containers.Map();
		fast_enabled = false;
	end

	methods
		function ts = TransSyst(n_s, n_a)
			% Create a TransSyst with n_s states and n_a actions
			ts.n_s = n_s;
			ts.n_a = n_a;
			ts.state1 = [];
			ts.state2 = [];
			ts.action = [];
		end

		function create_fast(ts)
			% Compute and store quick-access backward and forward
			% transition maps
			ts.fast_forward = containers.Map(num2cell(1:ts.n_s * ts.n_a), ...
											cell(1, ts.n_s * ts.n_a));
			ts.fast_backward = containers.Map(num2cell(1:ts.n_s * ts.n_a), ...
											cell(1, ts.n_s * ts.n_a));
			for i=1:ts.num_trans()
				s1 = ts.state1(i);
				s2 = ts.state2(i);
				a = ts.action(i);
				ts.fast_forward((a-1)*ts.n_s + s1) = ...
					[ts.fast_forward((a-1)*ts.n_s + s1) s2];
				ts.fast_backward((a-1)*ts.n_s + s2) = ...
					[ts.fast_backward((a-1)*ts.n_s + s2) s1];
			end
			ts.fast_enabled = true;
		end


		function add_transition(ts, s1, s2, a)
			% Add transition from s1 to s2 under action a
			ts.state1(end+1) = s1;
			ts.state2(end+1) = s2;
			ts.action(end+1) = a;
			ts.fast_enabled = false;
		end

		function split_state(ts, s1)
			% Split state number n:
			% - remove all transitions pertaining to s1
			for i=ts.num_trans():-1:1
				if ts.state1(i) == s1 || ts.state2(i) == s1
					ts.state1(i) = [];
					ts.state2(i) = [];
					ts.action(i) = [];
				end
			end

			% - move last (outside) state forward
			for i=1:ts.num_trans()
				if ts.state2(i) == ts.n_s
					ts.state2(i) = ts.n_s+1;
				end
			end

			% - update progress groups
			for i = 1:length(ts.pg_G)
				if ismember(s1, ts.pg_G{i})
					ts.pg_G{i} = union(ts.pg_G{i}, ts.n_s);
				end
			end

			% - increase state counter
			ts.n_s = ts.n_s + 1;
			ts.fast_enabled = false;
		end

		function add_progress_group(ts, U, G)
			% Add progress group G under modes U
			ts.pg_U{end+1} = U;
			ts.pg_G{end+1} = G;
		end

		function ret = num_trans(ts)
			% Number of transitions
			ret = length(ts.state1);
		end

		function ret = post(ts, q, a)
			% Compute the post set of state q and action a
			ret = [];
			if ts.fast_enabled
				ret = ts.fast_forward((a-1) * ts.n_s + q);
			else
				for i = 1:ts.num_trans()
					if q == ts.state1(i) && a == ts.action(i)
						ret(end+1) = ts.state2(i);
					end
				end
			end
			ret = unique(ret);
		end

		function [ret, K] = pre(ts, X, U, quant1, quant2)
			% Compute pre(X) under (quant1, quant2)-controllability
			% and action set U
            K = containers.Map('KeyType', 'uint64', 'ValueType', 'any');

			all_pre = [];

			if ts.fast_enabled
				for i=1:length(X)
					for j=1:length(U)
						all_pre = [all_pre ts.fast_backward((U(j)-1) * ts.n_s + X(i))];
					end
				end
			else
				for i=1:ts.num_trans()
					if ismember(ts.state2(i), X) && ismember(ts.action(i), U)
						all_pre(end+1) = ts.state1(i);
					end
				end
			end

			all_pre = unique(all_pre); 

			if strcmp(quant1, 'exists') && strcmp(quant2, 'exists')
				ret = all_pre;
				return
			end

			if strcmp(quant1, 'exists') && strcmp(quant2, 'forall')
				ret = [];
				for i = 1:length(all_pre)
					q = all_pre(i);
					act_list = zeros(1, ts.n_a);
					for a = 1:ts.n_a
						aPost = ts.post(q, a);
						act_list(a) = all(ismember(aPost, X)) && length(aPost) > 0;
					end
					if any(act_list)
						ret(end+1) = q;
                        K(q) = find(act_list); 
					end
				end
                assert(all(ismember(ret, cell2mat(K.keys))))
				return
			end

			if strcmp(quant1, 'forall') && strcmp(quant2, 'forall')
				ret = [];
				for i = 1:length(all_pre)
                    q = all_pre(i);
					act_list = zeros(1, ts.n_a);
					for a = 1:ts.n_a
						aPost = ts.post(q, a);
						act_list(a) = all(ismember(aPost, X)) && length(aPost) > 0;
					end
					if all(act_list)
						ret(end+1) = q;
                        K(q) = 1:ts.n_a;
					end
				end
				return
			end

			if strcmp(quant1, 'forall') && strcmp(quant2, 'exists')
				ret = [];
                for i = 1:length(all_pre)
                    q = all_pre(i);
					act_list = zeros(1, ts.n_a);
					for a = 1:ts.n_a
						act_list(a) = any(ismember(ts.post(q, a), X));
					end
					if all(act_list)
						ret(end+1) = q;
					end
				end
				return
			end
		end

		function [W, K] = pginv(ts, U, G, Z, B, quant1)
			% Compute U-controlled set
			% contained in G \cap B \setdiff Z that can be
			% used to force a transition to Z using the progress
			% group (U,G) under (quant1, forall)-controllability

			if strcmp(quant1, 'forall') && ~isempty(setdiff(1:ts.n_a, U))
				W = [];
                K = containers.Map('KeyType', 'uint64', 'ValueType', 'any');
				return
			end

			W = setdiff(intersect(G, B), Z);

			while true
                [preW, K] = ts.pre(union(W, Z), U, quant1, 'forall');
				Wt = intersect(W, preW);
    			if length(W) == length(Wt)
					break
				end
				W = Wt;
			end

            assert(all(ismember(W, cell2mat(K.keys))))
		end

		function [V, K] = win_until(ts, B, P, quant1)
			% Compute the winning set of
      		%  B U P
    		% under (quant1, forall)-controllability
    		V = [];
            K = containers.Map('KeyType', 'uint64', 'ValueType', 'any');
    		while true
                [preV, preK] = ts.pre(V, 1:ts.n_a, quant1, 'forall');
                K = [K; preK];
    			Vt = union(P, intersect(B, preV));
				Vt = reshape(Vt, 1, length(Vt));
    			for i=1:length(ts.pg_U)
    				% Progress groups
                    [preVinv, preKinv] = ts.pginv(ts.pg_U{i}, ts.pg_G{i}, V, B, quant1);
                    K = [K; preKinv];
    				Vt = union(Vt, preVinv);
    				Vt = reshape(Vt, 1, length(Vt));
    			end
    			if length(V) == length(Vt)
					break
				end
				V = Vt;
    		end
            if ~ all(ismember(setdiff(V, P), cell2mat(K.keys)))
                V
                P
                K.keys
            end
            assert(all(ismember(setdiff(V, P), cell2mat(K.keys))))
    	end

    	% Should be replaced when dual algos implemented!
    	function [V, K] = win_always(ts, B, quant1)
    		% Compute the winning set of
    		%  [] B
    		% under (quant1, forall)-controllability
    		V = B;
    		while true
                [preV, K] = ts.pre(V, 1:ts.n_a, quant1, 'forall');
    			Vt = intersect(V, preV);
    			if length(V) == length(Vt)
    				break
    			end
    			V = Vt;
    		end
    	end

    	function [ret, K] = win_until_and_always(ts, A, B, P, quant1)
    		% Compute the winning set of
		    %   []A && B U P
		    % under (quant1, forall)-controllability
		    if ~isempty(setdiff(1:ts.n_s, A))
		    	[Vinv, Kinv] = ts.win_always(A, quant1);
			 	[ret, K] = ts.win_until(intersect(B, Vinv), intersect(P, Vinv), quant1);
                for i=1:length(ret)
                    if ~ismember(ret(i), P)
                        K(ret(i)) = intersect(K(ret(i)), Kinv(ret(i)));
                    end
                end
		    else
		    	[ret, K] = ts.win_until(B, P, quant1);
		    end
	    end

    	function [V, Klist] = win_intermediate(ts, A, B, P, C_list, quant1)
    		% Compute winning set of
     		%  []A && ( (B U P) || [] (B &&_i <>C_i) )
    		% under (quant1, forall)-controllability
    		V = 1:ts.n_s;
            Klist = {};
    		while true
    			Vt = V;
    			[preV, ~] = ts.pre(V, 1:ts.n_a, quant1, 'forall');
    			for i=1:length(C_list)
    				Qi = union(P, intersect(intersect(B, C_list{i}), preV));
    				Qi = reshape(Qi, 1, length(Qi));
                    [Vti, Ki] = ts.win_until_and_always(A, B, Qi, quant1);
    				Vt = intersect(Vt, Vti);
                    Klist{i} = Ki;
    			end
    			
    			if length(V) == length(Vt)
    				break
    			end
    			V = Vt;
    		end
    	end

    	function [Vlist, Klist] = win_primal(ts, A, B, C_list, quant1, V)
		    % Compute winning set of
      		%  []A && <>[]B &&_i []<>C_i
    		% under (quant1, forall)-controllability
    		
    		if nargin<6
                V = [];
            end

            Vlist = {};
            Klist = {};

    		while true
    			Z = ts.pre(V, 1:ts.n_a, quant1, 'forall');
    			for i=1:length(ts.pg_U)
    				Z = union(Z, ...
    						  ts.pginv(ts.pg_U{i}, ts.pg_G{i}, V, A, quant1));
    			end
    			[Vt, Kt] = ts.win_intermediate(A, B, Z, C_list, quant1);

    			if length(Vt) == length(V)
    				break
    			end
                
                Vlist{end+1} = Vt;
                Klist{end+1} = Kt;

    			V = Vt;
    		end
    	end
	end
end