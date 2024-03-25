function str = double2latex(a, precision)
% Given a double, returns a string 
%
% Nicolas Boumal, June 5, 2016

    if ~exist('precision', 'var') || isempty(precision)
        precision = 2;
    end

    format = ['%' sprintf('.%de', precision)];
    
    parts = strsplit(sprintf(format, a), 'e');
    
    % Simplify the notation of the exponent
    parts{2} = sprintf('%d', str2double(parts{2}));
    
    str = [parts{1} ' \cdot 10^{' parts{2} '}'];

end
