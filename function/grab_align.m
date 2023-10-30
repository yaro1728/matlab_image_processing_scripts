function [cdk2_trace, other_trace] = grab_align(varargin)

% grab cells within a certain time range, then align to point of interest(mitosis, drug addition, mitosis skip, R point)
% 
% calls cellinfo.m
% Mingwei Min 2017.5.16

%% initialize parser
p = inputParser;

% add the function name to the input parser
p.FunctionName = mfilename;

% default tag/value pairs
p.addParameter('datapath',       '\\spencerstorage.int.colorado.edu\LabShare\IXMicroImages-goodNames\Yao\20170802_175010_097\Data_MM\', @isstr);
p.addParameter('rowmat',         1,                        @isnumeric);
p.addParameter('colmat',         1,                        @isnumeric);
p.addParameter('sitemat',        1,                        @isnumeric);
p.addParameter('cyto_col',       10,                        @isnumeric);
p.addParameter('nuc_col',        11,                        @isnumeric);
p.addParameter('other_col',      7,                        @isnumeric);%timelapse 3 color if 6=median, if 7=mean
p.addParameter('framerate',      5,                        @isnumeric);
p.addParameter('drug_frame',     1,                        @isnumeric);
p.addParameter('movie_leng',     354,                      @isnumeric);
p.addParameter('grab_time',      [-20 50],               @isnumeric);
p.addParameter('grab_cdk2state', 'all',                    @isstr);
p.addParameter('POI_name',       'mitosis',                @isstr);
p.addParameter('grab_generation','current',                @isstr);
p.addParameter('POI2_name',      'drug_addition',          @isstr);

% parse the variable argument list
p.parse(varargin{:});

% set struct elements to short variable names 
datapath        = p.Results.datapath;
rowmat          = p.Results.rowmat;
colmat          = p.Results.colmat;
sitemat         = p.Results.sitemat;
cyto_col        = p.Results.cyto_col;
other_col       = p.Results.other_col;
nuc_col         = p.Results.nuc_col;
framerate       = p.Results.framerate;
movie_leng      = p.Results.movie_leng;
drug_frame      = p.Results.drug_frame;
grab_time       = p.Results.grab_time;
grab_cdk2state  = p.Results.grab_cdk2state;
POI_name        = p.Results.POI_name;
grab_generation = p.Results.grab_generation;
POI2_name       = p.Results.POI2_name;

%% initiate an empty matrix for storing CDK2 trace and other trace of interest
cdk2_trace  = ones(5000,2*movie_leng+1)*NaN;
other_trace = ones(5000,2*movie_leng+1)*NaN;

count = 1;

for row = rowmat
    for col = colmat
        for site = sitemat
            % get cell info
            shot = [num2str(row),'_', num2str(col), '_', num2str(site)];
            if ~exist([datapath,'cellinfo_',shot,'.mat'],'file')
                cells = cellinfo_V2(datapath, row, col, site, nuc_col, cyto_col, framerate);
            else
                load([datapath,'cellinfo_',shot,'.mat'],'cells');
            end
            
            % get cells within the category of interest
            cellOfInterest=[];
            switch grab_cdk2state
                case 'cdk2inc'
                    cdk_state       = cat(1, cells.cdk_state);
                    cellOfInterest  = find(cdk_state==1);
                    
                case 'cdk2low'
                    cdk_state       = cat(1, cells.cdk_state);
                    cellOfInterest  = find(cdk_state==0);
                    
                case 'cdk2emerge'
                    cdk_state       = cat(1, cells.cdk_state);
                    cellOfInterest  = find(cdk_state==2);
                    
                case 'mitosis_skip'
                    mitosis_skip    = cat(1, cells.mito_skip);
                    cellOfInterest  = find(mitosis_skip>0);
                    
                case 'all'
                    cellOfInterest  = 1:length(cells);
            end
                    
            % get point of interest, the point to align traces to
            switch POI_name
                case 'mitosis'
                    if strcmp(grab_generation, 'current')
                        POI = cat(1,cells.trace_start).*(cat(1,cells.mother)>0);
                    else
                        POI = cat(1,cells.trace_end);
                    end
                case 'mitosis_skip'
                    POI     = cat(1,cells.mito_skip);
                case 'R_point'
                    POI     = cat(1,cells.rpoint);
                case 'drug_addition'
                    POI     = drug_frame;
            end
            
            % get anthother point of interest, for determination of
            % grabbing interval
            switch POI2_name
                case 'drug_addition'
                    POI2     = drug_frame;                
                case 'mitosis'
                    if strcmp(grab_generation, 'current')
                        POI2 = cat(1,cells.trace_start).*(cat(1,cells.mother)>0);
                    else
                        POI2 = cat(1,cells.trace_end);
                    end
                case 'mitosis_skip'
                    POI2     = cat(1,cells.mito_skip);
                case 'R_point'
                    POI2     = cat(1,cells.rpoint);
            end
            
            % get cells within interested time range
            nomal_trace     = find(cat(1,cells.noisy)==0 & cat(1,cells.momQuality)>0);
            POI2_to_POI     = POI2 - POI;
            if length(POI2_to_POI)>1
                cellOfInterest2 = find(POI2_to_POI>=grab_time(1)*framerate & POI2_to_POI<=grab_time(2)*framerate & ~POI2==0 & ~POI==0);
                cellOfInterest  = intersect(intersect(cellOfInterest2,cellOfInterest),nomal_trace);
            end
            
            % get cells of generation of interest
            switch grab_generation
                case 'previous'
                    cells_mom       = cells(cellOfInterest);
                    cellOfInterest  = cat(1,cells_mom.mother);
                case 'current'
                    cellOfInterest  = cellOfInterest;
                case 'next'
                    cells_daughter  = cells(cellOfInterest);
                    cellOfInterest  = cat(1,cells_daughter.daughter);
            end
    
            %% align to point of interest
            if ~isempty(cellOfInterest)
                for i = 1:length(cellOfInterest)
                    cellid      = cellOfInterest(i); 
                    if length(POI2_to_POI)>1
                        raw         = lineage_linker(cellid, cells);
                        zero_point  = POI(cellid);
                    else
                        raw         = cells(cellid).raw;
                        zero_point  = POI;
                    end
                    % re-align traces
                    
                    new_start   = movie_leng+1-zero_point+1;
                    new_end     = movie_leng+1-zero_point+movie_leng;
                    cdk2_trace(count,new_start:new_end)     = raw(:,cyto_col)./raw(:,nuc_col);
                    other_trace(count,new_start:new_end)	= raw(:,other_col);
                    count = count +1;
                end
            end

        end
    end
end
cdk2_trace  = cdk2_trace(1:count-1,:);
other_trace = other_trace(1:count-1,:);
cdk2_trace(cdk2_trace==0)   = NaN;
other_trace(other_trace==0) = NaN;