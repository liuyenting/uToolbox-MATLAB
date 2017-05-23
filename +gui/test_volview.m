close all;
clearvars;

%% load demo data
n = 8;
magic_data = zeros([2*n, 2*n, n]);
for i = 1:n
    magic_data(:, :, i) = magic(2*n);
end

%% show the data
h = gui.volview.VolView('Title', 'Hello world!');
h.show(magic_data);
