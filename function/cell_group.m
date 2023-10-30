function [current_cc] = cell_group(datadir, row, col, site, nuc_col, cyto_col, framerate, movie_leng)

shot        = [num2str(row),'_', num2str(col), '_', num2str(site)];
if ~exist([datadir,'cellinfo_',shot,'.mat'],'file')
    cells   = cellinfo_V2(datadir, row, col, site, nuc_col, cyto_col, framerate);
else
    load([datadir,'cellinfo_',shot,'.mat'],'cells');
end

current_cc  = cell(1000,movie_leng);

for i = 1:length(cells)

    if cells(i).mother>0   
        index           = cells(i).first_mito;
        validcell       = find(cellfun('isempty',current_cc(:,index))==1);
        if isempty(validcell)
            validcell   = 1;
        end
        current_cc{validcell(1),index}  = cells(i);
    end
end

