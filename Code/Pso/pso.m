% Particle Swarm Optimization adjusted for version 4 of the parallel package
% fun: a reference to a funciton taking `x` and a list of arguments
% args: the list of arguments to `fun`
% ...
function [solution, residuum] = pso(fun, args, lb, ub, opt)
	% add missing fields to opt struct
	if nargin == 3
		opt = get_opt(length(lb));
	else		
		opt = get_opt(length(lb), opt);
	end
	
	swarm = get_swarm(lb, ub, opt);
  
	res = zeros(swarm.n_iter, 1);
	sol = zeros(swarm.n_iter, length(lb));
	iter = 1;

	while iter <= swarm.n_iter  
		% Randomization
	    r_loc = rand(size(swarm.x));
	    r_glob = rand(size(swarm.x));
    
	    % Velocity and new position using previously defined random variables
	    % v_new = a * v_alt + b_loc * r_loc * (p_loc - x_old) + b_glob * r_glob * (p_glob - x_old)
		swarm.v = ...
			swarm.a      * swarm.v + ...
			swarm.b_loc  * r_loc  .* (swarm.p_loc  - swarm.x) + ...
			swarm.b_glob * r_glob .* (swarm.p_glob - swarm.x);
		swarm.x = swarm.x + swarm.v;
		
		% Boundary handling
		switch (opt.boundary_handling)
			case 'nearest'
				% Set new position on the search space boundary if it is on the outside
				for i = 1:size(swarm.x, 2)
					swarm.x(:,i) = max(swarm.x(:, i), lb(i));
					swarm.x(:,i) = min(swarm.x(:, i), ub(i));
				end
			case 'nearest_and_velocity_correction'
				for i = 1:size(swarm.x, 2)
					% Set particles that are outside by onto the boundary and set their new velocity correspondingly
					% x_old = x_new - v
					% x_corrected = x_old + v_corrected
					%
					% v_corrected = x_corrected - x_old
					%             = x_corrected - x_new + v
					x_corrected = min(max(swarm.x(:, i), lb(i)), ub(i));
					swarm.v(:, i) = x_corrected - swarm.x(:, i) + swarm.v(:, i);
					swarm.x(:, i) = x_corrected;
				end
			case 'reinitialize'
				% Reinitialize particle randomly if it is outside of the search space
				for i = 1:size(swarm.x, 2)
					idx_outside = swarm.x(:, i) < lb(i);
					swarm.x(idx_outside, :) = rand(sum(idx_outside), length(lb)) * (ub(i) - lb(i)) + lb(i);

					idx_outside = swarm.x(:, i) > ub(i);
					swarm.x(idx_outside, :) = rand(sum(idx_outside), length(lb)) * (ub(i) - lb(i)) + lb(i);
				end
			otherwise
				error('Invalid boundary handling strategy');
		end

	    % Update local and global attractors by evaluating function fun at current position
		t1 = tic;
		swarm = update_attractors(swarm, fun, args, opt);
		if and(opt. parallel, toc(t1) < 1)
			opt.parallel = false;
			fprintf("Switched to serial execution since problem is small\n");
		end

		% Store current best result (e.g. best global attractor)
		res(iter) = min(swarm.y_p_glob);
		sol(iter,:) = swarm.p_glob(find(swarm.y_p_glob == res(iter), 1), :);

    folder = ["../../Results/", args{2}{1}.folderName];
    fid = fopen([folder, "/protokoll_pso.txt"], "a");
    fprintf(fid, "Res: %f ", res(iter));
    for j=1:length(sol(iter,:))
      fprintf(fid, "%s: %f | ", opt.parameter_names{j}, sol(iter,j));
    end
    fprintf(fid, "\n");
    fclose(fid);  
    
		% Logging, Visualization
		% Progress
		if (swarm.visu.progress.n != 0)
			if (mod(iter, floor(swarm.n_iter / swarm.visu.progress.n)) == 0)
				fprintf("Iteration %i/%i\nCurrent solution: ", iter, swarm.n_iter);
				fprintf("%f ", sol(iter, :));
				fprintf("\n");
				pause(0); % TODO: show(), drawnow()?
			end
		end

		% 2D plot
		if (swarm.visu.plot_2d.n != 0)
			if (mod(iter, floor(swarm.n_iter / swarm.visu.plot_2d.n)) == 0)
				plot_swarm_2d(swarm, lb, ub, swarm.visu.plot_2d.dim1, swarm.visu.plot_2d.dim2, opt.parameter_names);
				pause(0);
			end
		end

		% Course of global attractors
		if (swarm.visu.plot_y_p_glob.n != 0)
			if (mod(iter, floor(swarm.n_iter / swarm.visu.plot_y_p_glob.n)) == 0)
				figure(swarm.visu.plot_y_p_glob.figure_id);
				clf;				
				idx_nothing_happened = find(abs((res(1:iter) - res(iter)) / res(iter)) < 0.02, 1);
				if numel(idx_nothing_happened)
					hold on;
					semilogy([idx_nothing_happened, iter],...
					res([idx_nothing_happened, iter]), "linewidth", 4);
				end
				semilogy(res(1:iter), "k");
				xlim([0,length(res)]);
				pause(0);
			end
		end

		% Terminate if res does not change over the last 10% in relation to all Iterations
		idx_nothing_happened = find(abs((res(1:iter) - res(iter)) / res(iter)) < 0.02, 1);
		if numel(idx_nothing_happened)
			if iter - idx_nothing_happened > swarm.n_iter / 5
				fprintf("PSO terminiated due to small changes of res\n");
				res = res(1:iter);
				sol = sol(1:iter, :);
				break;
			end
		end
		iter = iter + 1;
	end

	if opt.visu
		fprintf("Final residuum: %f\n", res(end));
	end

	if nargout >= 1
		solution = sol(end, :);
	end
	if nargout >= 2
		residuum = res(end);
	end	
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function swarm = update_attractors(swarm, fun, args, opt)
	% Compute function value y at current position
	##  y = zeros(swarm.n, 1);
	##  for i = 1:swarm.n
	##    y(i) = fun(swarm.x(i, :));
	##  end
	% create cell array of argument cell arrays
	args_cells = cell(length(swarm.x),1);
	args_cells(:) = {args};
	data_cells = mat2cell(swarm.x, ones(1, swarm.n), size(swarm.x, 2));
	% run fun parallel or seriel
	if opt.parallel
		% number of cores to use - maximum of 30 to not overload the server
		nparallel = min(nproc()-3, 30);
		t1 = tic;
		% run in parallel
		% y = cellfun(fun, data_cells, args_cells);
		y = parcellfun(nparallel, fun, data_cells, args_cells, "VerboseLevel", 1);
		##      [a,b] = parcellfun (nparallel, runOptimOnce , x0Array, "UniformOutput", false);
	else
		if (opt.vectorized_function_evaluation == false)
			y = cellfun(fun, data_cells, args_cells);
		else
			y = fun(swarm.x, args);
		end
	end
	
	% Update current best local position p_loc and its value y_p_loc
	##	  for i = 1:swarm.n
	##	    if (y(i) < swarm.y_p_loc(i))
	##	      swarm.p_loc(i, :) = swarm.x(i, :);
	##	      swarm.y_p_loc(i) = y(i);
	##	    end
	##	  end
	##	tmp = min([swarm.y_p_loc, y], [], 2);
	idx_better = y < swarm.y_p_loc;
	swarm.p_loc(idx_better, :) = swarm.x(idx_better, :);
	swarm.y_p_loc(idx_better) = y(idx_better);
	
	% Update current best global position p_glob and its value y_p_glob
	for i = 1:swarm.n
		neighbors_of_i = swarm.neighbors.i(swarm.neighbors.j == i);
		[~, best_neighbor_of_i] = min(swarm.y_p_loc(neighbors_of_i));
		swarm.p_glob(i, :) = swarm.p_loc(neighbors_of_i(best_neighbor_of_i), :);
		swarm.y_p_glob(i) = swarm.y_p_loc(neighbors_of_i(best_neighbor_of_i));
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_swarm_2d(swarm, lb, ub, dim1, dim2, parameter_names)
	figure(swarm.visu.plot_2d.figure_id);
	clf;

	n_plots = length(dim1);
	n_plots_rows = round(sqrt(n_plots));
	n_plots_cols = ceil(n_plots / n_plots_rows);

	[~, idx_best_particle] = min(swarm.y_p_loc); % Best particle so far
	
	% Iterate dim/2 subplots
	for i = 1:n_plots
		subplot(n_plots_rows, n_plots_cols, i);
		hold on;

		% Bounds
		plot([lb(dim1(i)) ub(dim1(i)) ub(dim1(i)) lb(dim1(i)) lb(dim1(i))],...
			 [lb(dim2(i)) lb(dim2(i)) ub(dim2(i)) ub(dim2(i)) lb(dim2(i))], 'k');

		% Particles
		plot(swarm.x(:, dim1(i)), swarm.x(:, dim2(i)), '.k');
		plot(swarm.p_loc(idx_best_particle, dim1(i)), swarm.p_loc(idx_best_particle, dim2(i)), 'or');

		% Visible area
		factor = 1.1;
		dim1_limits = [lb(dim1(i)) ub(dim1(i))];
		dim1_limits = (dim1_limits - mean(dim1_limits)) * factor + mean(dim1_limits);
		dim2_limits = [lb(dim2(i)) ub(dim2(i))];
		dim2_limits = (dim2_limits - mean(dim2_limits)) * factor + mean(dim2_limits);
		xlim(dim1_limits);
		ylim(dim2_limits);

		% Labels
		xlabel(parameter_names(dim1(i)));
		ylabel(parameter_names(dim2(i)));

		% Other
		grid on;
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function opt = get_opt(n_dim, opt_in)
	opt.visu = true;
	opt.n_particles = 10;
	opt.n_iter = 100;
	opt.parallel = false;		
	opt.topology = "ring";
	opt.coupling = 3;
	
	opt.a = 0.72984;
	opt.b_loc = 1.496172;
	opt.b_glob = 1.496172;

	% Vectorized function evaluation is faster but not possible in certain cases, thus disabled by default
	opt.vectorized_function_evaluation = false;

	% Boundary handling
	%opt.boundary_handling = 'nearest';
	opt.boundary_handling = 'nearest_and_velocity_correction';
	%opt.boundary_handling = 'reinitialize';

	% Parameter names
	opt.parameter_names = cell(1, n_dim);
	for i = 1:n_dim
		opt.parameter_names(i) = ["Parameter ", num2str(i)];
	end
	
	% Copy fields that have been specified
	fields = fieldnames(opt_in);		
	for i=1:length(fields)
		opt = setfield(opt, fields{i}, getfield(opt_in, fields{i}));
	end
	
	if opt.parallel
		try
			pkg load parallel;
		catch
			opt.parallel = false;
		end_try_catch	
	end
end
