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
		fast_post = containers.Map();
		fast_pre = containers.Map();
		fast_enabled = false;

        % Debugging flag
        b_debug = false
	end

	methods
		function ts = TransSyst(n_s, n_a)
			% Create a TransSyst with n_s states and n_a actions
			ts.n_s = uint32(n_s);
			ts.n_a = uint32(n_a);
			ts.state1 = uint32([]);
			ts.state2 = uint32([]);
			ts.action = uint32([]);
		end

		function create_fast(ts)
			% Compute and store quick-access backward and forward
			% transition maps
            if ts.fast_enabled
                return
            end
			ts.fast_post = containers.Map(num2cell(1:ts.n_s * ts.n_a), ...
											cell(1, ts.n_s * ts.n_a));
			ts.fast_pre = containers.Map(num2cell(1:ts.n_s), ...
											cell(1, ts.n_s));

			for i=1:ts.num_trans()
				s1 = ts.state1(i);
				s2 = ts.state2(i);
				a = ts.action(i);
  				ts.fast_post((a-1)*ts.n_s + s1) = ...
   					[ts.fast_post((a-1)*ts.n_s + s1) s2];

   				ts.fast_pre(s2) = [ts.fast_pre(s2) s1];
			end

			ts.fast_enabled = true;
		end


		function add_transition(ts, s1, s2, a)
			% Add transition from s1 to s2 under action a

            if ts.b_debug
                assert(1 <= s1 && s1 <= ts.n_s)
                assert(1 <= s2 && s2 <= ts.n_s)
                assert(1 <= a && a <= ts.n_a)
            end

			ts.state1(end+1) = s1;
			ts.state2(end+1) = s2;
			ts.action(end+1) = a;
			ts.fast_enabled = false;
		end

		function split_state(ts, s1)
			% Split state number n:
			% remove all transitions pertaining to s1
			for i=ts.num_trans():-1:1
				if ts.state1(i) == s1 || ts.state2(i) == s1
					ts.state1(i) = [];
					ts.state2(i) = [];
					ts.action(i) = [];
				end
			end

			% move last (outside) state forward
			for i=1:ts.num_trans()
				if ts.state2(i) == ts.n_s
					ts.state2(i) = ts.n_s+1;
				end
			end

			% update progress groups
			for i = 1:length(ts.pg_G)
				if ismember(s1, ts.pg_G{i})
					ts.pg_G{i} = union(ts.pg_G{i}, ts.n_s);
				end
			end

			% increase state counter
			ts.n_s = ts.n_s + 1;
			ts.fast_enabled = false;
		end

		function add_progress_group(ts, U, G)
			% Add progress group G under modes U

            if ts.b_debug
                assert(all(1 <= U) && all(U <= ts.n_a))
                assert(all(1 <= G) && all(G <= ts.n_s))
            end

			ts.pg_U{end+1} = U;
			ts.pg_G{end+1} = G;
		end

		function ret = num_trans(ts)
			% Number of transitions
			ret = length(ts.state1);
		end
	end
end