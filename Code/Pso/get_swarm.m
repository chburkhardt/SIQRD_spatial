function swarm = get_swarm(lb, ub, opt)
	swarm = [];
	dim = length(lb);
	
	% Number of iterations
	swarm.n_iter = opt.n_iter;
	
	% Number of particles
	swarm.n = opt.n_particles;
	
	% Swarm parameters:
	% v_old: Previous speed
	% v_new: Next speed
	% b_loc: Spring stiffness to local attractor (scalar)
	% b_glob: Spring stiffness to global attractor (scalar)
	% r_loc: Randomization of local attractor (vector)
	% r_glob: Randomization of global attractor (vector)
	% p_loc: best solution of the particle itself (vector)
	% p_glob: best solution found by the particle and its neighbors (vector)
	%
	% Equation:
	% v_new = a * v_old + b_loc * r_loc * (p_loc - x_old) + b_glob * r_glob * (p_glob - x_old)
	% x_new = x_old + v_new
	swarm.a = opt.a;
	swarm.b_loc = opt.b_loc;
	swarm.b_glob = opt.b_glob;
	
	% Position and velocity
	swarm.x = rand(swarm.n, dim) .* repmat(ub - lb, swarm.n, 1) + repmat(lb, swarm.n, 1);
	swarm.v = zeros(swarm.n, dim);
	
	% Local and global attractors p and their values y
	swarm.p_loc = swarm.x;
	swarm.p_glob = swarm.x;
	swarm.y_p_loc = inf(swarm.n, 1);
	swarm.y_p_glob = inf(swarm.n, 1);
	
	% Neighborhood topologies:
	switch opt.topology		
		case "full" % Complete graph
			[swarm.neighbors.i, swarm.neighbors.j, ~] = find(ones(swarm.n));
		case "ring" % Ring
			coupling = opt.coupling; % itself and its neighbors
			if !mod(coupling,2)
				error("Only symmetric coupling allowed");
			end			
			swarm.neighbors.i = zeros(1, swarm.n * coupling);
			swarm.neighbors.j = zeros(1, swarm.n * coupling);			
			for k = 1:swarm.n
				swarm.neighbors.i(1+(k-1)*coupling:k*coupling) = k;
				swarm.neighbors.j(1+(k-1)*coupling:k*coupling) =...
				mod((k-(coupling-1)/2:k+(coupling-1)/2)-1, swarm.n)+1;
			end			
	end

	if opt.visu
		% Visualization
		% Print progress n times during the computation.
		swarm.visu.progress.n = min(20, swarm.n_iter);
		
		% Plot particle position in dim1 and dim2 n times during the computation.
		% dim/2 subplots will be displayed in order to show every dimension.
		swarm.visu.plot_2d.dim1 = [1:2:dim];
		swarm.visu.plot_2d.dim2 = [2:2:dim];
		if (mod(dim, 2) == 1)
			% Append dummy dimension if dimension is odd
			swarm.visu.plot_2d.dim2(end + 1) = 1;
		end
		swarm.visu.plot_2d.n = min(20, swarm.n_iter);
		swarm.visu.plot_2d.figure_id = 1;
		
		% Plot course of best residual n times during the computation.
		swarm.visu.plot_y_p_glob.n = min(100, swarm.n_iter);
		swarm.visu.plot_y_p_glob.figure_id = 3;
	else
		% Visualization
		% Print progress n times during the computation.
		swarm.visu.progress.n = 0;
		
		% Plot particle position in dim1 and dim2 n times during the computation.
		swarm.visu.plot_2d.n = 0;
		swarm.visu.plot_2d.figure_id = 1;
		
		% Plot course of best residual n times during the computation.
		swarm.visu.plot_y_p_glob.n = 0;
		swarm.visu.plot_y_p_glob.figure_id = 3;		
	end
end
