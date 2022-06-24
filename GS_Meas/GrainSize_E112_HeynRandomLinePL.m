function [G_PL, MIC, intersection_count, nlines, total_line_length] = GrainSize_E112_HeynRandomLinePL(ebsd, varargin)

[~, G_PL, ~, MIC, ~, ~, ~, ~, ~, nlines, total_line_length] = grainsize_linint_random(ebsd, 50, varargin{:});

intersection_count = total_line_length / MIC;  

end