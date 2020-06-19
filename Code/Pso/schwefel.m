function y = schwefel(x)
	if iscell(x)
		x = x{1};
	end
	y = sum(-x .* sin(sqrt(abs(x))), 2) + (size(x, 2) * 418.982887272434);
end
