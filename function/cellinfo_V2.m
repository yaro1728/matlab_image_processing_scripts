function cells = cellinfo_V2(datadir, row, col, site, nuc_col, cyto_col, framerate)

% version 2, added in 
% 1. test for mitosis skip (potential senescence entry)
% 2. and identification of R point
% 3. cdk2 state of mother
% 4. test if the trace is noisy
% Mingwei Min, 2017.5.17

% addition in 2017.5.23: set NA for mother quality instead of empty

shot        = [num2str(row),'_', num2str(col), '_', num2str(site)];
load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');

cellnum     = length(genealogy);
cells       = struct([]);
callthresh  = 0.5;
callthresh2 = 0.7;

for i=1:cellnum
    cells(i).raw    = tracedata(i,:,:);
    cells(i).raw    = squeeze(cells(i).raw);
    raw             = cells(i).raw;
    trace_length    = sum(~isnan(raw(:,1)));
    trace_start     = find(~isnan(raw(:,1)),1,'first');
    trace_end       = find(~isnan(raw(:,1)),1,'last');
    cells(i).cdk2_trace     = raw(:,cyto_col)./raw(:,nuc_col);
    cdk2_trace      = cells(i).cdk2_trace;
    cdk2_trace(trace_start:trace_end)   = smooth(cdk2_trace(trace_start:trace_end));
    
    %% set default sensor quality
    cells(i).momQuality = 2; % 2 as the sensor quality cannot be decided    
    
    %% test if the CDK2 sensor is lost
    if (nanmean(raw(:,cyto_col))+nanmean(raw(:,nuc_col)))>20
        cells(i).sensor  = 1;
    else
        cells(i).sensor  = 0;
    end
    
    %% test if the cdk2 trace is noisy or beyond normal range
    if max(cdk2_trace)>3.5 | min(cdk2_trace)<0 | max(diff(cdk2_trace))>0.3
        cells(i).noisy  = 1;
    else
        cells(i).noisy  = 0;
    end
    
    %% if the cell is a mother
    if ismember(i,genealogy)
        cells(i).last_mito  = trace_end+1;
        cells(i).daughter   = find(genealogy==i);
        
        % if it's a mother, test if the sensor is working fine
        mom_quality_check   = min(trace_start, max(1,trace_end - framerate*3));
        cells(i).momQuality = max(cdk2_trace(mom_quality_check:end))>1;
        
        % if it's a mother, it does not skip mitosis
        cells(i).mito_skip  = 0;
    else
        cells(i).last_mito  = 0;
        cells(i).daughter   = 0;
        % if it's not a monther, check if it goes from
        % G2-like to G1-like CDK2 state without a mitosis
        if cells(i).sensor
            cells(i).mito_skip = mitosis_skip_check(cdk2_trace, raw(:,nuc_col), raw(:,3),framerate,~isnan(genealogy(i)));
        else
            cells(i).mito_skip = 0;
        end
        
    end

    %% if the cell is a daughter
    if ~isnan(genealogy(i))
        cells(i).first_mito  = trace_start;
        cells(i).mother      = genealogy(i);
        cells(i).mom_cdk2    = cells(genealogy(i)).cdk_state;
    % catergorize CDK2 group  
        if cells(genealogy(i)).momQuality==0 % if mother sensor fail
           categorization       = 3; % do not catergorize
           cells(i).momQuality  = 0; % its own sensor fail too
        else
            calltime    = trace_start + 3*framerate;
            if calltime>trace_end
                categorization  = 3;
            else            
                categorization  = (cdk2_trace(calltime))>callthresh;
                if categorization == 1 
                    test_point_min = min(calltime+4*framerate,trace_end);
                    if cdk2_trace(test_point_min)<callthresh2
                        categorization = 3;
                    end
                else
                    testpoint = calltime + framerate;
                    while testpoint<trace_end
                        if cdk2_trace(testpoint) > 0.7
                            categorization = 2;
                        end
                        testpoint = testpoint + framerate;
                    end
                end
            end
        end
    else
        cells(i).first_mito = 0;
        cells(i).mother     = 0;
        categorization      = 3;
    end

    %% test for R point
    if categorization>0
        cells(i).rpoint = rpoint_check_V2(cdk2_trace);
    else
        cells(i).rpoint = 0;
    end
    
    cells(i).trace_lenght   = trace_length;
    cells(i).trace_start    = trace_start;
    cells(i).trace_end      = trace_end;
    cells(i).cdk_state      = categorization;
end
save([datadir,'cellinfo_',shot,'.mat'],'cells');
