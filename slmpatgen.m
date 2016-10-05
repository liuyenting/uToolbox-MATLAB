function slmpatgen
%SLMPATGEN Generates a pattern for the SLM module.
% Fill in all the values or load them from a .mat file, click 'Apply' to
% generate a new pattern. Click 'Save' to save the image in a .bmp file.

close all;

% Create the figure window.
f = figure('Position', [0, 0, 1070, 640]);
    set(f, 'Name', 'SLM Pattern Generator', 'NumberTitle', 'off');
    set(f, 'Visible', 'off');
    set(f, 'MenuBar', 'none', 'ToolBar', 'none');
    
% Viewer.
hViewer = axes('Units', 'pixels', 'Box', 'on', 'Position', [20, 20, 800, 600]);
    %set(hViewer, 'XColor', get(f, 'Color'));
    set(hViewer, 'XTick', []);
    %set(hViewer, 'YColor', get(f, 'Color'));
    set(hViewer, 'YTick', []);

% SLM parameters.
pSlm = uipanel('Units', 'pixels', 'Position', [830, 546, 220, 80]);
    set(pSlm, 'Title', 'SLM', 'FontSize', 14);
    
    tSlmSize = uicontrol('Style', 'text', 'Parent', pSlm, 'Position', [10, 40, 90, 20]);
        set(tSlmSize, 'String', 'Resolution', 'FontSize', 14, 'HorizontalAlignment', 'right');
    hSlmPxSz = uicontrol('Style', 'popupmenu', 'Parent', pSlm, 'Position', [110, 40, 100, 20]);
        set(hSlmPxSz, 'String', {'QXGA', 'Full HD'}, 'FontSize', 14);
    
    tSlmPxSz = uicontrol('Style', 'text', 'Parent', pSlm, 'Position', [10, 10, 90, 20]);
        set(tSlmPxSz, 'String', 'Pixel Size', 'FontSize', 14, 'HorizontalAlignment', 'right');
    hSlmPxSz = uicontrol('Style', 'edit', 'Parent', pSlm, 'Position', [110, 10, 70, 20]);
        set(hSlmPxSz, 'FontSize', 14);
    tSlmPxSzUnit = uicontrol('Style', 'text', 'Parent', pSlm, 'Position', [180, 10, 40, 20]);
        set(tSlmPxSzUnit, 'String', 'um', 'FontSize', 14, 'HorizontalAlignment', 'left');

% Optics parameters.
pOptics = uipanel('Units', 'pixels', 'Position', [830, 426, 220, 110]);
    set(pOptics, 'Title', 'Optics', 'FontSize', 14);
    
    tOpWavLen = uicontrol('Style', 'text', 'Parent', pOptics, 'Position', [10, 100, 90, 20]);
        set(tOpWavLen, 'String', 'Wavelength', 'FontSize', 14, 'HorizontalAlignment', 'right');
    hOpWavLen = uicontrol('Style', 'edit', 'Parent', pOptics, 'Position', [110, 100, 70, 20]);
        set(hOpWavLen, 'FontSize', 14);
    tOpWavLenUnit = uicontrol('Style', 'text', 'Parent', pOptics, 'Position', [180, 100, 40, 20]);
        set(tOpWavLenUnit, 'String', 'nm', 'FontSize', 14, 'HorizontalAlignment', 'left');
    
    tOpMag = uicontrol('Style', 'text', 'Parent', pOptics, 'Position', [10, 70, 90, 20]);
        set(tOpMag, 'String', 'Magnitude', 'FontSize', 14, 'HorizontalAlignment', 'right');
    hOpMag = uicontrol('Style', 'edit', 'Parent', pOptics, 'Position', [110, 70, 70, 20]);
        set(hOpMag, 'FontSize', 14);
    
    tOpNA = uicontrol('Style', 'text', 'Parent', pOptics, 'Position', [10, 25, 60, 20]);
        set(tOpNA, 'String', 'Aperture', 'FontSize', 14, 'HorizontalAlignment', 'right');
        
        tOpApOd = uicontrol('Style', 'text', 'Parent', pOptics, 'Position', [70, 40, 30, 20]);
            set(tOpApOd, 'String', 'OD', 'FontSize', 14, 'HorizontalAlignment', 'right');
        hOpApOd = uicontrol('Style', 'edit', 'Parent', pOptics, 'Position', [110, 40, 70, 20]);
            set(hOpApOd, 'FontSize', 14);
        tOpApOdUnit = uicontrol('Style', 'text', 'Parent', pOptics, 'Position', [180, 40, 40, 20]);
            set(tOpApOdUnit, 'String', 'mm', 'FontSize', 14, 'HorizontalAlignment', 'left');
        
        tOpApId = uicontrol('Style', 'text', 'Parent', pOptics, 'Position', [70, 10, 30, 20]);
            set(tOpApId, 'String', 'ID', 'FontSize', 14, 'HorizontalAlignment', 'right');
        hOpApId = uicontrol('Style', 'edit', 'Parent', pOptics, 'Position', [110, 10, 70, 20]);
            set(hOpApId, 'FontSize', 14);
        tOpApIdUnit = uicontrol('Style', 'text', 'Parent', pOptics, 'Position', [180, 10, 40, 20]);
            set(tOpApIdUnit, 'String', 'mm', 'FontSize', 14, 'HorizontalAlignment', 'left');

% Pattern parameters.
pPat = uipanel('Units', 'pixels', 'Position', [830, 246, 220, 170]);
    set(pPat, 'Title', 'Pattern', 'FontSize', 14);
    
    tPatN = uicontrol('Style', 'text', 'Parent', pPat, 'Position', [10, 130, 90, 20]);
        set(tPatN, 'String', 'N Beams', 'FontSize', 14, 'HorizontalAlignment', 'right');
    hPatN = uicontrol('Style', 'edit', 'Parent', pPat, 'Position', [110, 130, 70, 20]);
        set(hPatN, 'FontSize', 14);
        
    tPatSh = uicontrol('Style', 'text', 'Parent', pPat, 'Position', [10, 85, 70, 20]);
        set(tPatSh, 'String', 'Shift', 'FontSize', 14, 'HorizontalAlignment', 'right');
        
        tPatShX = uicontrol('Style', 'text', 'Parent', pPat, 'Position', [80, 100, 20, 20]);
            set(tPatShX, 'String', 'X', 'FontSize', 14, 'HorizontalAlignment', 'right');
        hPatShX = uicontrol('Style', 'edit', 'Parent', pPat, 'Position', [110, 100, 70, 20]);
            set(hPatShX, 'FontSize', 14);
        tPatShXUnit = uicontrol('Style', 'text', 'Parent', pPat, 'Position', [180, 100, 40, 20]);
            set(tPatShXUnit, 'String', 'um', 'FontSize', 14, 'HorizontalAlignment', 'left');
        
        tPatShY = uicontrol('Style', 'text', 'Parent', pPat, 'Position', [80, 70, 20, 20]);
            set(tPatShY, 'String', 'Y', 'FontSize', 14, 'HorizontalAlignment', 'right');
        hPatShY = uicontrol('Style', 'edit', 'Parent', pPat, 'Position', [110, 70, 70, 20]);
            set(hPatShY, 'FontSize', 14);
        tPatShYUnit = uicontrol('Style', 'text', 'Parent', pPat, 'Position', [180, 70, 40, 20]);
            set(tPatShYUnit, 'String', 'um', 'FontSize', 14, 'HorizontalAlignment', 'left');
    
    tPatTilt = uicontrol('Style', 'text', 'Parent', pPat, 'Position', [10, 40, 90, 20]);
        set(tPatTilt, 'String', 'Tilt', 'FontSize', 14, 'HorizontalAlignment', 'right');
    hPatTilt = uicontrol('Style', 'edit', 'Parent', pPat, 'Position', [110, 40, 70, 20]);
        set(hPatTilt, 'FontSize', 14);
    tPatTiltUnit = uicontrol('Style', 'text', 'Parent', pPat, 'Position', [180, 40, 40, 20]);
        set(tPatTiltUnit, 'String', 'deg', 'FontSize', 14, 'HorizontalAlignment', 'left');
    
    tPatTh = uicontrol('Style', 'text', 'Parent', pPat, 'Position', [10, 10, 90, 20]);
        set(tPatTh, 'String', 'Threshold', 'FontSize', 14, 'HorizontalAlignment', 'right');
    hPatTh = uicontrol('Style', 'slider', 'Parent', pPat, 'Position', [110, 10, 90, 20]);

hApply = uicontrol('Style', 'pushbutton', 'Position', [830, 80, 100, 20]);
    set(hApply, 'String', 'Apply', 'FontSize', 14);
    
hSave = uicontrol('Style', 'pushbutton', 'Position', [830, 50, 100, 20]);
    set(hSave, 'String', 'Save', 'FontSize', 14);

hExit = uicontrol('Style', 'pushbutton', 'Position', [830, 20, 100, 20]);
    set(hExit, 'String', 'Exit', 'FontSize', 14);
    
% Show the settled figure.
f.Visible = 'on';

end