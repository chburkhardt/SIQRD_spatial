## Copyright (C) 2020 andre
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} testBeta (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: andre <andre@DESKTOP-OAO29OE>
## Created: 2020-06-23

function retval = testBeta (input1, input2)
parameterArray = read_parameter("../../Results/parameter.txt"); 
parameter = parameterArray{1};
darkFigure = parameter.darkFigure;

courses = sir_eqn_spatial("totalCourseOfDisease", parameter);
lags = sir_eqn_spatial("lags");
plot(lags(1:end-1), courses.infective);
nu = 0.345;

betaFun = @(tau) exp(-nu*tau).*interp1(lags(1:end-1), courses.infective, tau, "linear");

intsum = 1/integral(betaFun, 1, 28);
fprintf("1/int = %f\n", intsum);

fprintf("integral of gamma_i = %f\n", sum(filter([1,1]/2,1, courses.infective)));

endfunction
