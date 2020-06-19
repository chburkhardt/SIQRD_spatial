function y = sphere_fun(x)
	if iscell(x)
		x = x{1};
	end
	y = sum(x.^2, 2);
end
