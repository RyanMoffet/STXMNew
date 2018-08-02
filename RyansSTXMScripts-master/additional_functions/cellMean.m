function meanArray = cellMean(c, weight)
% meanArray = cellMean(c, [weight=1])
% mean over the elements of a cell c, keeping matrix structures of cell
% elements etc. Use weight if given.

% based on http://stackoverflow.com/q/5197692/321973, courtesy of gnovice
% (http://stackoverflow.com/users/52738/gnovice)
% extended to weighted averaging by Tobias Kienzler
% (see also http://stackoverflow.com/q/5231406/321973)

dim = ndims(c{1});          %# Get the number of dimensions for your arrays
if ~exist('weight', 'var') || isempty(weight); weight = 1; end;
eins = ones(size(c{1})); % that's german for "one", creative, I know...
if ~iscell(weight)
    % ignore length if all elements are equal, this is case 1
    if isequal(weight./max(weight(:)), ones(size(weight)))
        weight = repmat(eins, [size(eins)>0 length(c)]);
    elseif isequal(numel(weight), length(c)) % case 2: per cell-array weigth
        weight = repmat(shiftdim(weight, -3), [size(eins) 1]);
    else
        error(['Weird weight dimensions: ' num2str(size(weight))]);
    end
else % case 3, insert some dimension check here if you want
    weight = cat(dim+1,weight{:});
end;

M = cat(dim+1,c{:});        %# Convert to a (dim+1)-dimensional matrix
sumc = sum(M.*weight,dim+1);
sumw = sum(weight,dim+1);
meanArray = sumc./sumw;  %# Get the weighted mean across arrays