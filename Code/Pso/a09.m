function a09
	tic
	
	%problem_name = "sphere";
	%problem_name = "rosenbrock";
	%problem_name = "rastrigin";
	problem_name = "schwefel";
	dim = 16;
	
	problem = get_problem(problem_name, dim);
	
	opt.parallel = true;
	opt.visu = true;
	opt.n_particles = 2000;
	opt.n_iter = 2000;
	opt.a = 0.72984;
	opt.b_loc = 1.496172;
	opt.b_glob = 1.496172;
	opt.coupling = 3;

	opt.vectorized_function_evaluation = true;

	lb = problem.bounds.lower * ones(1, dim);
	ub = problem.bounds.upper * ones(1, dim);
	
	[solution, residuum] = pso(problem.f, lb, ub, opt);
	fprintf("Calculation took %is\nfunction min was: %1.12f\nSolution\n", toc, residuum);
	fprintf("%f; ", solution);
	fprintf("\n");
end
