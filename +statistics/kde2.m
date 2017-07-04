function [bw, d, vx, vy] = kde2(data, varargin)
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
%   over a grid of the same size as DATA, where DATA is the binning result.
%   X, Y are the meshgrid over which the variable D is computed.
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
nsmpl = sum(data(:));

%% parse the input
p = inputParser;
addOptional(p, 'Size', [nrows, ncols]);
addParameter(p, 'XLim', [1, ncols]);
addParameter(p, 'YLim', [1, nrows]);
parse(p, varargin{:});

%TODO: use Size arguments
xlim = p.Results.XLim;
ylim = p.Results.YLim;

%% pad the data to a square with power of 2 
n = 2^ceil(log2(max(nrows, ncols)));
[data, anchor] = image.impad(data, [n, n], 'center');

%% DCT of the initial data
a = dct2d(data);

%% square of the optimal bandwidth
scaling = [xlim(2)-xlim(1), ylim(2)-ylim(1)];

I = (0:n-1).^2;
a2 = a.^2;

t_star = root(@(t)(t-evolve(t)), nsmpl);
    function t = root(f, N)
    %ROOT Find the smallest root when there are more than one exist.

        N = 50*(N<=50)+1050*(N>=1050)+N*((N<1050)&(N>50));
        tol = 10^-12+0.01*(N-50)/1000;
        flag = 0;
        while flag == 0
            try
                t = fzero(f, [0, tol]);
                flag = 1;
            catch
                % double the search interval
                tol = min(tol*2, .1); 
            end
            
            % if all the trial fails...
            if tol <= .1
                t = fminbnd(@(x)(abs(f(x))), 0, .1); 
                flag = 1;
            end
        end
    end

    function [out, t] = evolve(t0)
        sum_func = func([0, 2], t0) + func([2, 0], t0) + 2*func([1, 1], t0);
        t = (2*pi*nsmpl*sum_func)^(-1/3);
        out = (t0-t)/t;
    end

p_02 = func([0, 2], t_star);
p_20 = func([2, 0], t_star); 
p_11 = func([1, 1], t_star);
t_y = ( p_02^(3/4) / (4*pi*nsmpl*p_20^(3/4)*(p_11+sqrt(p_20*p_02))) )^(1/3);
t_x = ( p_20^(3/4) / (4*pi*nsmpl*p_02^(3/4)*(p_11+sqrt(p_20*p_02))) )^(1/3);
    function out = func(s, t_func)
    %FUNC The iterative function for the diffusion equation.
    
        if sum(s) <= 4
            sum_func = func([s(1)+1, s(2)], t_func) + func([s(1), s(2)+1], t_func); 
            const = ( 1 + 1/2^(sum(s)+1) ) / 3;
            t = (-2*const*K(s(1))*K(s(2))/nsmpl/sum_func)^(1/(2+sum(s)));

            out = psi(s, t);
        else
            out = psi(s, t_func);
        end
        
            %
            % helper functions
            %
            function y = K(x)
                y = (-1)^x * prod((1:2:2*x-1)) / sqrt(2*pi);
            end

            function out = psi(s, t)
                w = exp(-I*pi^2*t) .* [1, (1/2)*ones(1, n-1)];
                % x/y weighting
                wx = w.*(I.^s(1));
                wy = w.*(I.^s(2));
                
                out = (-1)^sum(s) * (wy*a2*wx') * pi^(2*sum(s));
            end
    end

% smooth the DCT of initial data using t_star
a_t = exp(- I.' * pi^2 * (t_x+t_y) / 2) .* a;

% the bandwidth
bw = sqrt([t_x, t_y]) .* scaling;

%% calculate the density
if nargout > 1
    d = idct2d(a_t) * (numel(a_t)/prod(scaling));
    
    % crop the result
    d = d(anchor(1):anchor(1)+nrows-1, anchor(2):anchor(2)+ncols-1);
    
    % remote negative values
    d(d < 0) = eps;
    
    % generate the grid
    xrange = xlim(1):scaling(1)/(ncols-1):xlim(2);
    yrange = ylim(1):scaling(2)/(nrows-1):ylim(2);
    [vx, vy] = meshgrid(xrange, yrange);
end

end 

function data = dct2d(data)
%DCT2D Computes the 2-D discrete cosine transform.

[nrows,ncols] = size(data);
if nrows ~= ncols
    error('statistics:kde2:dct2d', 'Data is not a square array.');
end

% compute weights for DFT coefficients
w = [1; 2*(exp(-1i*(1:nrows-1)*pi/(2*nrows))).'];
weight = w(:, ones(1, ncols));

data = dct1d(dct1d(data)')';
    function y = dct1d(x)
        % re-order the elements of the input data
        x = [x(1:2:end, :); x(end:-2:2, :)];

        % multiply FFT by weights
        y = real(weight .* fft(x));
    end
end

function data = idct2d(data)
%IDCT2D Computes the 2-D inverse discrete cosine transform.

[nrows, ncols] = size(data);
if nrows ~= ncols
    error('statistics:kde2:idct2d', 'Data is not a square array.');
end

% compute weights for DFT coefficients
w = exp(1i*(0:nrows-1)*pi/(2*nrows)).';
weights = w(:,ones(1,ncols));

data = idct1d(idct1d(data)');
    function y=idct1d(x)
        xp = real(ifft(weights .* x));
        
        % create filler
        y = zeros(nrows, ncols);
        
        % re-order the columns of the output results
        y(1:2:nrows, :) = xp(1:nrows/2, :);
        y(2:2:nrows, :) = xp(nrows:-1:nrows/2+1, :);
    end
end
