% a helper function providing an interface for optCovid
% used for parcellfun in e.g. pso.m
function output = evalSirLocalOneArgumentLocal(x, args)
  [paraNames, parameterArray, RKIdata] = args{:};
  %addpath '../Spatial/Optimize';
  output = norm(optCovidLocal(x, paraNames, parameterArray, RKIdata, false));
endfunction
