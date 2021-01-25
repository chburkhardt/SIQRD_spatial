% a helper function providing an interface for optCovid
% used for parcellfun in e.g. pso.m
function output = evalSirLocalOneArgument(x, args)
  [paraNames, parameterArray, RKIdata] = args{:};
  addpath 'Optimize';
  output = norm(optCovid(x, paraNames, parameterArray, RKIdata, false));
endfunction
