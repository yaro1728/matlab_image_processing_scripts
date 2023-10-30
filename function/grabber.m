function [cellInfoOfInterest] = grabber(varargin)

% grab cells of interest within a certain time range (relative to mitosis, drug addition, mitosis skip, R point)
% the output of the function is struct with all feature from cellinfo.
%
% calls cellinfo.m
% Mingwei Min 2017.5.23

%% initialize parser
p = inputParser;

% add the function name to the input parser
p.FunctionName = mfilename;

% default tag/value pairs
p.addParameter('datapath',       '\\spencerstorage.int.colorado.edu\LabShare\IXMicroImages-goodNames\Yao\20170802_175010_097\Data_MM\', @isstr);
p.addParameter('rowmat',         1,                        @isnumeric);
p.addParameter('colmat',         1:12,                        @isnumeric);
p.addParameter('sitemat',        1,                        @isnumeric);
p.addParameter('cyto_col',       10,                        @isnumeric);
p.addParameter('nuc_col',        11,                        @isnumeric);
p.addParameter('other_col',      7,                        @isnumeric);
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

%% initiate an empty structure for information of cells of interest
cellInfoOfInterest  =  struct([]);

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
            cellOfInterest = [];
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
                        POI = cat(1,cells.first_mito);
                    else
                        POI = cat(1,cells.last_mito);
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
                        POI2 = cat(1,cells.first_mito);
                    else
                        POI2 = cat(1,cells.last_mito);
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

            cellInfoOfInterest = [cellInfoOfInterest,cells(cellOfInterest)];


        end
    end
end
