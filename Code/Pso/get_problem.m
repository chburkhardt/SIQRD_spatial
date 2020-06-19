function problem = get_problem(type, dim)
  problem = [];
  switch(type)
	case "sphere"
		problem.bounds.lower = -100;
		problem.bounds.upper =  100;
		problem.dim = dim;
		problem.f = @sphere_fun;
		problem.optimum.value = 0;
		problem.optimum.position = zeros(1, dim);
	case "rosenbrock"
		problem.bounds.lower = -30;
		problem.bounds.upper =  30;
		problem.dim = dim;
		problem.f = @rosenbrock;
		problem.optimum.value = 0;
		problem.optimum.position = ones(1, dim);
	case "rastrigin"
		problem.bounds.lower = -5.12;
		problem.bounds.upper =  5.12;
		problem.dim = dim;
		problem.f = @rastrigin;
		problem.optimum.value = 0;
		problem.optimum.position = zeros(1, dim);
	case "schwefel"
		problem.bounds.lower = -500;
		problem.bounds.upper =  500;
		problem.dim = dim;
		problem.f = @schwefel;
		problem.optimum.value = -dim * 418.9829;
		problem.optimum.position = 420.9687 * ones(1, dim);
  end
end
