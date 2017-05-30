close all;
clearvars;

I = imread('Averaged shifted histograms_20x.tif');

% sampling size
sz = 32;

%% create template
T = zeros(sz);

ht = figure('Name', 'Averaged NPC', 'NumberTitle', 'off');
ht.Visible = 'off'; 
ht.ToolBar = 'none';

% show the template
imagesc(T);
axis image;
colormap(gray);
% move to corner 
ht.Position = [0, 0, 256, 256];
movegui(ht, 'southeast');

ht.UserData.List = [];

ht.Visible = 'on';

%% show the selector 
hs = figure('Name', 'Select samples...', 'NumberTitle', 'off');
hs.Visible = 'off';

% show the image
imagesc(log(I));
axis image;
colormap(gray);
colorbar;

% set to the largest possible size
scsz = 0.7 * util.screensize;
imsz = size(I);
ratio = min(scsz / imsz);
hs.Position = [0, 0, imsz([2, 1])*ratio];
movegui(hs, 'center');

% settled, show the figure
hs.Visible = 'on';

% attach the data
hs.UserData.Data = I;
hs.UserData.Size = sz;
hs.UserData.Result = [];
hs.UserData.Viewer = ht;

%% select the samples
% register the listener
hs.WindowButtonDownFcn = @mouseDown;

%% support functions
function mouseDown(source, event)

hs = source;

% get mouse action
action = hs.SelectionType;
% identify right click
if strcmp(action, 'alt')
    p = get(gca, 'CurrentPoint');
    x = p(1, 1); 
    y = p(1, 2);
    fprintf('(x, y) = (%.3f, %.3f)\n', x, y);
    
    hold on;
    plot(x, y, 'yx');
    hold off;
    
    % get the region
    cI = cropregion(hs.UserData.Data, [x, y], hs.UserData.Size);
    
    % save the position
    ls = hs.UserData.Viewer.UserData.List;
    ls = [ls; x, y];
    hs.UserData.Viewer.UserData.List = ls;
    
    % find the center
    sz = size(cI);
    c = findcenter(cI);
    % find shift
    s = c - sz/2;
    % align the image
    cI = imtranslate(cI, -s, 'OutputView', 'same');
    
    % sum
    if isempty(hs.UserData.Result)
        hs.UserData.Result = cI;
    else
        I = hs.UserData.Result;
        %I = (I+cI) / 2;
        
        % align the image
        I = I + cI;
        hs.UserData.Result = I;
    end
    
    % show the result
    figure(hs.UserData.Viewer);
    imagesc(hs.UserData.Result);
    
    fprintf('%d samples averaged\n', size(hs.UserData.Viewer.UserData.List, 1));
    
    % refocus
    figure(hs);
end

end

function cI = cropregion(data, pos, sz)

pos = floor(pos-sz/2);
cI = data(pos(2):pos(2)+sz-1, pos(1):pos(1)+sz-1);

end

function c = findcenter(I)

%% convert to particle list
% number of pixels
np = numel(I);
% number of potential dots
nd = sum(I(:));

% create the list
list = zeros([1, nd]);

id = 1;
% iterate through the pixels
for ip = 1:np
    % repeatance
    n = I(ip);
    
    if n > 0
        % record the position
        list(id:id+n-1) = ip;
        % increase the index
        id = id+n;
    end
end
% records should match the pixel sum
assert(nd == (id-1));

% convert the assignment to coordinates
[y, x] = ind2sub(size(I), list);

% pack as a single variable
data = [x.', y.'];

%% calculate the kde
sz = size(I);
[~, density, ~, ~] = kde2d( ...
    data, ...
    max(sz), ...    % sampling grid size (squared)
    [0, 0], sz ...  % sampling range
); 

%% k-mean
[vx, vy] = meshgrid(1:sz(1), 1:sz(2));
X = [vx(:), vy(:), density(:)];
[~, c] = kmeans(X, 1);

c = c(1:2);

end