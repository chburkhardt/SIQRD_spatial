function y = rastrigin(x)
	if iscell(x)
		x = x{1};
	end
	y = 10 * size(x, 2) + sum(x.^2 - 10*cos(2*pi*x), 2);
end
