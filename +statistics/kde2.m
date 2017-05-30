function [bw, d, x, y] = kde2(data, varargin)
%KDE2 Bivariate kernel density estimator with diagonal bandwidth matrix.
%
%   (From the original author) Fast and accurate state-of-the-art bivariate
%   kernel density estimator with diagonal bandwidth matrix. The kernel is
%   assumed to be Gaussian. The two bandwidth parameters are chosen
%   optimally without ever assuming a parametric model for the data or any
%   rule-of-thumb. Unlike many other procedures, this one is immune to
%   accuracy failures in the estimation of multi-modal densities with
%   widely separated modes.
%
%   BW = KDE2(DATA) calculates the kernel density estimation with the two
%   optimal bandwidth for a bivariate Gaussian kernel of format [bandwith
%   X, bandwith Y].
%   [BW, D, X, Y] = KDE2(DATA) calculates the bandwidth and the density D
%   over a grid of the same size as DATA. X, Y are the meshgrid over which
%   the variable D is computed.
%   [...] = KDE2(..., PARMS) contains additional configurations.
%
%   Parameters
%   ----------
%   'Size'  Control the output size, default to the size of DATA.
%   'XLim'  X bounding box over which the density is computed, default to
%           the range of DATA.
%   'YLim'  Y bounding box over which the density is computed, default to
%           the range of DATA.
%   
%   Reference
%   ---------
%   Kernel density estimation via diffusion
%       Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010)
%       Annals of Statistics, Volume 38, Number 5, pages 2916-2957.

[nrows, ncols] = size(data);

%% parse the input
p = inputParser;
addOptional(p, 'Size', [nrows, ncols]);
addParameter(p, 'XLim', [0, ncols]);
addParameter(p, 'YLim', [0, nrows]);
parse(p, varargin{:});

%% pad the data to a square with power of 2 
n = 2^ceil(log2(max(nrows, ncols)));
data = image.impad(data, [n, n], 'center');

end 
