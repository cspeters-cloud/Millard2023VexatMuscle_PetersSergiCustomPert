%%
% SPDX-FileCopyrightText: 2023 Matthew Millard <millard.matthew@gmail.com>
%
% SPDX-License-Identifier: MIT
%
%
%%

function figPub = plotRatFibrilTRSS2017(...
                    figPub,...
                    ratFibrilModelsDefault,...
                    ratFibrilModelsFitted,...
                    fittingInfo,...
                    benchRecordFitted,...
                    expTRSS2017,...
                    simConfig,...
                    fittingConfig,...
                    plotConfig,...
                    pubPlotOptions)




%
% Plot the force-length relation
%


expfl.lce = [];
expfl.fN   = []; 

mdlfl.lce = [];
mdlfl.fN   = [];

assert(ratFibrilModelsDefault(1).curves.useCalibratedCurves==1,...
       'Error: the calibrated curves should be used');

curveFalN = ratFibrilModelsDefault(1).curves.activeForceLengthCurve; 
lceOptMdl = ratFibrilModelsDefault(1).musculotendon.optimalFiberLength;    

for idxTrial = simConfig.trials
    lce  = expTRSS2017.activeLengtheningData(idxTrial).x(1,1);
    fN   = expTRSS2017.activeLengtheningData(idxTrial).y(1,1);

    expfl.lce  = [expfl.lce;lce];
    expfl.fN    = [expfl.fN;fN];   

    fNMdl = calcBezierYFcnXDerivative(lce./lceOptMdl,curveFalN,0);
    mdlfl.lce  = [mdlfl.lce;lce];
    mdlfl.fN    = [mdlfl.fN;fNMdl];
            
end    

%
% fal-fitting
%
subplot('Position',reshape(plotConfig.subPlotPanel(1,1,:),1,4));

    curveFitted = ratFibrilModelsFitted(1).curves.activeForceLengthCurve;    
    curveOrig = ratFibrilModelsDefault(1).curves.activeForceLengthCurve;

    sampleFlNOrig=calcBezierYFcnXCurveSampleVector(...
                    curveOrig,100,curveOrig.xEnd);
    
    sampleFlNFitted=calcBezierYFcnXCurveSampleVector(...
                    curveFitted,100,curveFitted.xEnd);

    plot(sampleFlNOrig.x,sampleFlNOrig.y,'--',...
         'Color',[1,1,1].*0.75,'DisplayName','$$f^L(\ell^M)$$');
    hold on;
    plot(sampleFlNFitted.x,sampleFlNFitted.y,'-',...
         'Color',[1,1,1].*0.75,'DisplayName','$$f^L(\ell^M) (fitted)$$');
    hold on;

    for idx=simConfig.trials        
        plot(expfl.lce(idx)./lceOptMdl,...
             expfl.fN(idx),...
             '+','Color',[0,0,0],...
             'MarkerFaceColor',plotConfig.lineColors.exp(idx,:),...
             'DisplayName',...
             expTRSS2017.activeLengtheningData(idx).seriesName);
        hold on;
    end
    xlabel('Norm. Length ($$\ell/\ell_o$$)');
    ylabel(expTRSS2017.activeLengtheningData(3).yName);
    title('Force-Length Relation Fitting');  
    xlim([pubPlotOptions.lceNMin,pubPlotOptions.lceNMax]);
    box off;
    hold on;   

%
% fv-fitting
%        
subplot('Position',reshape(plotConfig.subPlotPanel(1,2,:),1,4));

    plot(sampleFlNFitted.x,sampleFlNFitted.y,'-',...
         'Color',[1,1,1].*0.75,'DisplayName','$$f^L(\ell^M)$$');
    hold on;        

    curveFvN = ratFibrilModelsDefault(1).curves.fiberForceVelocityCurve;
    curveFvNFitted = ratFibrilModelsFitted(1).curves.fiberForceVelocityCurve;

    fvN = calcBezierYFcnXDerivative(0.11,curveFvN,0);
    fvNFitted = calcBezierYFcnXDerivative(0.11,curveFvNFitted,0);
    
    plot(sampleFlNFitted.x,sampleFlNFitted.y.*fvN,'--',...
         'Color',[1,1,1].*0.75,'DisplayName','$$f^L(\ell^M)f^V(v^M)$$');
    hold on;  
    plot(sampleFlNFitted.x,sampleFlNFitted.y.*fvNFitted,'-',...
         'Color',[1,1,1].*0.75,'DisplayName','$$f^L(\ell^M)f^V(v^M)$$ (fitted)');
    hold on;  

    for idx=simConfig.trials
        hdlVis='off';
        if(idx==1)
            hdlVis = 'on';
        end

        plot(expTRSS2017.activeLengtheningData(idx).x./lceOptMdl,...
             expTRSS2017.activeLengtheningData(idx).y,...
             '-','Color',plotConfig.lineColors.exp(idx,:),...
             'DisplayName','$$f^{EXP}$$ TRSS2017',...
             'HandleVisibility',hdlVis);
        hold on; 
        idxKey = expTRSS2017.activeLengtheningData(idx).keyIndices(...
                    fittingConfig.idxFvKey);
        plot(expTRSS2017.activeLengtheningData(idx).x(idxKey)./lceOptMdl,...
             expTRSS2017.activeLengtheningData(idx).y(idxKey),...
             'x','Color',[0,0,0],...
             'HandleVisibility','off');
        hold on; 
    end

    xlabel('Norm. Length ($$\ell/\ell_o$$)');
    ylabel(expTRSS2017.activeLengtheningData(3).yName);
    title('Force-Velocity Curve Fitting');  
    xlim([pubPlotOptions.lceNMin,pubPlotOptions.lceNMax]);
    box off;
    hold on;           

%
% Experimental analysis plot
%

subplot('Position',reshape(plotConfig.subPlotPanel(2,1,:),1,4));
    fill([min(sampleFlNFitted.x);...
          max(sampleFlNFitted.x);...
          fliplr(sampleFlNFitted.x')'],...
          [0;...
           0;...
          fliplr(sampleFlNFitted.y')'],...
          [1,1,1].*0.85,'EdgeAlpha',0,...
          'HandleVisibility','off');
    hold on;  

for idx=simConfig.trials
    hVis='off';
    if(idx==1)
        hVis='on';
    end
    plot(expTRSS2017.activeLengtheningData(idx).x./lceOptMdl,...
         expTRSS2017.activeLengtheningData(idx).y,...
         '-','Color',plotConfig.lineColors.exp(idx,:),...
         'DisplayName','TRSS2017: $$f^{EXP}_i$$',...
         'HandleVisibility',hVis);
    hold on;  
    text(expTRSS2017.activeLengtheningData(idx).x(end)./lceOptMdl,...
         expTRSS2017.activeLengtheningData(idx).y(end),...
         sprintf('%s%i%s','$$f^{EXP}_{',idx,'}$$'),...
         'HorizontalAlignment','left',...
         'VerticalAlignment','baseline',...
         'FontSize',8);
    hold on;
end


for idx=simConfig.trials     
    subplot('Position',reshape(plotConfig.subPlotPanel(2,1,:),1,4)); 

        idxKey=expTRSS2017.activeLengtheningData(idx).keyIndices(...
                fittingConfig.idxFvKey);
        lceNKey = expTRSS2017.activeLengtheningData(idx).x(idxKey)/lceOptMdl;
        fceNKey = expTRSS2017.activeLengtheningData(idx).y(idxKey); 

        fvCurve = ratFibrilModelsFitted(idx).curves.fiberForceVelocityCurve;
        falCurve = ratFibrilModelsFitted(idx).curves.activeForceLengthCurve;
        lceNV = benchRecordFitted.normFiberLength(:,idx);
        vceNV = benchRecordFitted.normFiberVelocity(:,idx);



        fvNV    = zeros(size(lceNV));
        falNV   = zeros(size(lceNV));
        for i=1:1:length(lceNV)
            falNV(i,1)=calcBezierYFcnXDerivative(lceNV(i,1),falCurve,0);
            fvNV(i,1)=calcBezierYFcnXDerivative(vceNV(i,1),fvCurve,0);
        end
        fxeNV = fvNV.*falNV;

        txtName = expTRSS2017.activeLengtheningData(idx).seriesName;
        i0=strfind(txtName,'Exp.');
        txtName(1,i0:4)='Sim.';


        idxValid = find(lceNV >= lceNKey);

        hdlVis='off';
        if(idx==1)
            hdlVis='on';                
        end
        

        plot(lceNV(idxValid),...
            fxeNV(idxValid),...
             '--','Color',plotConfig.lineColors.simXE(idx,:),...
             'DisplayName','Crossbridges: $$f^{XE}_i=f^L(\ell)f^V(v)$$',...
             'HandleVisibility',hdlVis);
        hold on;   
        text(lceNV(idxValid(1,1)),...
            fxeNV(idxValid(1,1)),...
            sprintf('%s%i%s','$$f^{XE}_{',idx,'}$$'),...
            'FontSize',8,...
            'HorizontalAlignment','left',...
            'VerticalAlignment','top');
        hold on

        idxKeyF0 = expTRSS2017.activeLengtheningData(idx).keyIndices(...
                        fittingConfig.idxFlKey);
        idxKeyF2 = expTRSS2017.activeLengtheningData(idx).keyIndices(...
                         fittingConfig.idxFvKey);

        plot(expTRSS2017.activeLengtheningData(idx).x(idxKeyF0)./lceOptMdl,...
             expTRSS2017.activeLengtheningData(idx).y(idxKeyF0),...
             'o','Color',[0,0,0],'HandleVisibility','off');
        hold on;
        plot(expTRSS2017.activeLengtheningData(idx).x(idxKeyF2)./lceOptMdl,...
             expTRSS2017.activeLengtheningData(idx).y(idxKeyF2),...
             'x','Color',[0,0,0],'HandleVisibility','off');
        hold on;
        
        xlim([pubPlotOptions.lceNMin,pubPlotOptions.lceNMax]);

        box off;


    subplot('Position',reshape(plotConfig.subPlotPanel(2,1,:),1,4));                

        idxKey=expTRSS2017.activeLengtheningData(idx).keyIndices(...
                fittingConfig.idxFvKey);
        lceNKey = expTRSS2017.activeLengtheningData(idx).x(idxKey)/lceOptMdl;
        fceNKey = expTRSS2017.activeLengtheningData(idx).y(idxKey);           
              
        lceNVU= unique(lceNV);

        lceNExp = expTRSS2017.activeLengtheningData(idx).x ./ lceOptMdl;
        [lceNExpU,iu] = unique(lceNExp);
        fceNExp = expTRSS2017.activeLengtheningData(idx).y;
        fceNExpU = fceNExp(iu);

        titinAnalysis(idx).lceN = lceNVU;
        titinAnalysis(idx).fceN = interp1(lceNExpU,fceNExpU,lceNVU,'linear');

        

        [lceNU,iU]=unique(lceNV);
        fxeNU=fxeNV(iU);

        idxZeroTitin = find(titinAnalysis(idx).lceN <= lceNKey);  
        titinAnalysis(idx).fxeN = interp1(lceNU,fxeNU,titinAnalysis(idx).lceN);
        titinAnalysis(idx).f2N  = titinAnalysis(idx).fceN-titinAnalysis(idx).fxeN;
        titinAnalysis(idx).f2N(idxZeroTitin)=0;
        titinAnalysis(idx).f2N(titinAnalysis(idx).f2N<=0)=0;

        hVis='off';
        if(idx==1)
            hVis='on';
        end
        plot(titinAnalysis(idx).lceN, ...
             titinAnalysis(idx).f2N,...
             '-','Color',plotConfig.lineColors.calcTitinF(idx,:),...
             'DisplayName','Titin: $$f^T_i=f^{EXP}_i-f^{XE}_i$$',...
             'HandleVisibility',hVis);

        text(titinAnalysis(idx).lceN(end),...
             titinAnalysis(idx).f2N(end),...
             sprintf('%s%i%s','$$f^{T}_{',idx,'}$$'),...
             'HorizontalAlignment','left',...
             'VerticalAlignment','baseline',...
             'FontSize',8);
              
        hold on;
        box off;
        xlim([pubPlotOptions.lceNMin,pubPlotOptions.lceNMax]);
        ylim([0,pubPlotOptions.fceNMax]);

        xlabel('Norm. Length ($$\ell/\ell_o$$)');
        ylabel(expTRSS2017.activeLengtheningData(3).yName);
        title({'Estimated crossbridge and titin forces','during active-lengthening'});

        if(idx==3)
            legend('Location','NorthWest');
        end

    subplot('Position',reshape(plotConfig.subPlotPanel(2,2,:),1,4));  

        idxKValid = find(titinAnalysis(idx).f2N >...
                    pubPlotOptions.stiffnessLowerForceBound);

        xSpan = max(titinAnalysis(idx).lceN)-min(titinAnalysis(idx).lceN);
        [sp,f2NS,rho] = spaps(...
            titinAnalysis(idx).lceN(idxKValid),...
            titinAnalysis(idx).f2N(idxKValid),1e-6);

        f2NSFit = fnval(sp,titinAnalysis(idx).lceN(idxKValid));
        titinAnalysis(idx).f2NS=titinAnalysis(idx).f2N;
        titinAnalysis(idx).f2NS(idxKValid)=f2NSFit;


        plot(titinAnalysis(idx).lceN,...
             titinAnalysis(idx).f2N,...
             '-','Color',plotConfig.lineColors.calcTitinF(idx,:));
        hold on;
        plot(titinAnalysis(idx).lceN,...
             titinAnalysis(idx).f2NS,...
             '-','Color',[1,1,1].*0.75);
        hold on;
        box off;
        xlabel('Norm. Length ($$\ell/\ell_o$$)');
        ylabel('Norm. Force ($$f/f_o$$)');
        title('Calc. Titin forces vs. smoothed version');
            

    subplot('Position',reshape(plotConfig.subPlotPanel(3,1,:),1,4));  

        titinAnalysis(idx).k2N = ...
            calcCentralDifferenceDataSeries(...
                titinAnalysis(idx).lceN,...
                titinAnalysis(idx).f2N);

        titinAnalysis(idx).k2NS = ...
            calcCentralDifferenceDataSeries(...
                titinAnalysis(idx).lceN,...
                titinAnalysis(idx).f2NS);
        lceNSpan = max(titinAnalysis(idx).lceN)...
                  -min(titinAnalysis(idx).lceN);

        hVis='off';
        if(idx==1)
            hVis='on';
        end

        plot(titinAnalysis(idx).lceN-titinAnalysis(idx).lceN(1,1), ...
             titinAnalysis(idx).f2N,...
             '-','Color',plotConfig.lineColors.calcTitinF(idx,:),...
             'DisplayName',...
             sprintf('%s%i%s','$$f^{T}_{',idx,'}$$'));
        hold on;

        text(lceNSpan,titinAnalysis(idx).f2N(end),...
             sprintf('%1.1f%s',titinAnalysis(idx).f2N(end),'$$f_o$$'),...
             'HorizontalAlignment','right',...
             'VerticalAlignment','bottom',...
             'FontSize',8);
        hold on;

        box off;

        xlabel('Norm. Length Change $$(\ell/\ell_o)-\ell_i$$');
        ylabel('Norm. Force ($$f/f_o$$)');
        title('Estimated titin force-length relation');  
        
        xlim([0,lceNSpan]);
        ylim([0,pubPlotOptions.fceNMax]);        
        if(idx==3)
            legend('Location','NorthWest');
        end
    subplot('Position',reshape(plotConfig.subPlotPanel(3,2,:),1,4));
    
        idxKValid = find(titinAnalysis(idx).f2N >...
                    pubPlotOptions.stiffnessLowerForceBound);

        if(pubPlotOptions.useSmoothedStiffnessData==1)
            plot(titinAnalysis(idx).lceN(idxKValid(2:end))-titinAnalysis(idx).lceN(1,1), ...
                 titinAnalysis(idx).k2NS(idxKValid(2:end)),...
                 '-','Color',plotConfig.lineColors.calcTitinK(idx,:),...
                 'DisplayName',...
                 sprintf('%s%i%s','$$k^{T}_{',idx,'}$$'));
            hold on;
            text(lceNSpan,titinAnalysis(idx).k2NS(end),...
                 sprintf('%1.1f',titinAnalysis(idx).k2NS(end)),...
                 'HorizontalAlignment','right',...
                 'VerticalAlignment','bottom',...
                 'FontSize',8);
            hold on;                
        end
        if(pubPlotOptions.plotRawStiffnessData==1 ...
                || pubPlotOptions.useSmoothedStiffnessData==0)

            seriesColor = plotConfig.lineColors.calcTitinK(idx,:);
            if(pubPlotOptions.plotRawStiffnessData==1)
                seriesColor = [1,1,1].*0.5;
            end

            plot(titinAnalysis(idx).lceN(idxKValid)-titinAnalysis(idx).lceN(1,1), ...
                 titinAnalysis(idx).k2N(idxKValid),...
                 '-','Color',seriesColor,...
                 'DisplayName',...
                 sprintf('%s: %1.1f %s','$$\ell_i$$',...
                        titinAnalysis(idx).lceN(1,1),'$$\ell_o$$'));
            hold on;
            text(lceNSpan,titinAnalysis(idx).k2N(end),...
                 sprintf('%1.1f',titinAnalysis(idx).k2N(end)),...
                 'HorizontalAlignment','right',...
                 'VerticalAlignment','bottom',...
                 'FontSize',8);
            hold on;  
         
        end
        
        box off;
        xlim([0,lceNSpan]);
        ylim([0,pubPlotOptions.kceNMax]);

        xlabel('Norm. Length Change $$(\ell/\ell_o)-\ell_i$$');
        ylabel('Norm. Slope $$(f/f_o)/(\ell/\ell_o)$$');
        title('Estimated titin stiffness-length relation');  
        if(idx==3)
            legend('Location','SouthEast');
        end

end

%
% Generate plots that compare the simulation results to the experiments
%

% load(fullfile(projectFolders.output_structs_TRSS2017,...
%         ['benchRecordVexat_TRSS2017_fitted',fittingConfig.trialStr,'.mat']));
% tmp = load(fullfile(projectFolders.output_structs_FittedModels,...
%         ['ratTRSS2017EDLFibrilActiveTitinFitted',fittingConfig.trialStr,'.mat']));
% ratFibrilModelsFitted=tmp.ratFibrilModelsFitted;



curveFitted = ratFibrilModelsFitted(1).curves.activeForceLengthCurve;    
sampleFlNFitted=calcBezierYFcnXCurveSampleVector(...
                curveFitted,100,curveFitted.xEnd);

subplot('Position',reshape(plotConfig.subPlotPanel(4,1,:),1,4));     
    fill([min(sampleFlNFitted.x);...
          max(sampleFlNFitted.x);...
          fliplr(sampleFlNFitted.x')'],...
          [0;...
           0;...
          fliplr(sampleFlNFitted.y')'],...
          [1,1,1].*0.85,'EdgeAlpha',0,...
          'HandleVisibility','off');
    hold on; 

for idx=simConfig.trials     
    hdlVis='off';
    if(idx==1)
        hdlVis = 'on';
    end

    plot(expTRSS2017.activeLengtheningData(idx).x./lceOptMdl,...
         expTRSS2017.activeLengtheningData(idx).y,...
         '-','Color',plotConfig.lineColors.exp(idx,:),...
         'DisplayName','$$f^{EXP}$$ TRSS2017',...
         'HandleVisibility',hdlVis);
    hold on; 

    hVis='off';
    if(idx==1)
        hVis='on';
    end
    plot(titinAnalysis(idx).lceN, ...
         titinAnalysis(idx).f2N,...
         '-','Color',plotConfig.lineColors.calcTitinF(idx,:),...
         'DisplayName','Titin: $$f^T_i=f^{EXP}_i-f^{XE}_i$$',...
         'HandleVisibility',hVis);

    text(titinAnalysis(idx).lceN(end),...
         titinAnalysis(idx).f2N(end),...
         sprintf('%s%i%s','$$f^{T}_{',idx,'}$$'),...
         'HorizontalAlignment','left',...
         'VerticalAlignment','baseline',...
         'FontSize',8);        
    
   plot( benchRecordFitted.normFiberLength(:,idx), ...
         benchRecordFitted.normFiberForce(:,idx),...
         '-','Color',plotConfig.lineColors.simF(idx,:),...
         'DisplayName','Sim: $$f^{*}_i$$',...
         'HandleVisibility',hVis);

   hold on;

   plot( benchRecordFitted.normFiberLength(:,idx), ...
         benchRecordFitted.normDistalTitinForce(:,idx),...
         '-','Color',plotConfig.lineColors.simTitinF(idx,:),...
         'DisplayName','Sim: $$f^{T*}_i$$',...
         'HandleVisibility',hVis);

   hold on;

end


box off;
xlim([pubPlotOptions.lceNMin,pubPlotOptions.lceNMax]);
ylim([0,pubPlotOptions.fceNMax]); 
xlabel('Norm. Length ($$\ell/\ell_o$$)');
ylabel(expTRSS2017.activeLengtheningData(3).yName);
title({'Simulated and estimated crossbridge','and titin forces during active-lengthening'});


subplot('Position',reshape(plotConfig.subPlotPanel(5,1,:),1,4));     


    for idx=simConfig.trials

        plot(titinAnalysis(idx).lceN-titinAnalysis(idx).lceN(1,1), ...
             titinAnalysis(idx).f2N,...
             '-','Color',plotConfig.lineColors.calcTitinF(idx,:),...
             'DisplayName',...
             sprintf('%s%i%s','$$f^{T}_{',idx,'}$$'));
        hold on;

        text(lceNSpan,titinAnalysis(idx).f2N(end),...
             sprintf('%1.1f%s',titinAnalysis(idx).f2N(end),'$$f_o$$'),...
             'HorizontalAlignment','right',...
             'VerticalAlignment','bottom',...
             'FontSize',8);
        hold on;

    end

    for idx=simConfig.trials

       

       plot( benchRecordFitted.normFiberLength(:,idx)...
              -benchRecordFitted.normFiberLength(1,idx), ...
             benchRecordFitted.normDistalTitinForce(:,idx),...
             '-','Color',plotConfig.lineColors.simTitinF(idx,:),...
             'DisplayName',sprintf('%s%i%s','Sim: $$f^{T*}_',idx,'$$'),...
             'HandleVisibility','on');

       hold on;            

       lceNSpan = benchRecordFitted.normFiberLength(end,idx)...
                 -benchRecordFitted.normFiberLength(1,idx);
       text(lceNSpan,...
            benchRecordFitted.normDistalTitinForce(end,idx),...
             sprintf('%1.1f%s',benchRecordFitted.normDistalTitinForce(end,idx),'$$f_o$$'),...
             'HorizontalAlignment','left',...
             'VerticalAlignment','top',...
             'FontSize',8);
       hold on;


    end    

box off;
xlim([0,lceNSpan]);
ylim([0,pubPlotOptions.fceNMax]);

xlabel('Norm. Length Change $$(\ell/\ell_o)-\ell_i$$');
ylabel('Norm. Force ($$f/f_o$$)');
title('Estimated and simulated titin force-length relation');  

legend('Location','NorthWest');


subplot('Position',reshape(plotConfig.subPlotPanel(5,2,:),1,4));     


for idx=simConfig.trials

    idxKValid = find(titinAnalysis(idx).f2N >...
                pubPlotOptions.stiffnessLowerForceBound);

    if(pubPlotOptions.useSmoothedStiffnessData==1)
        plot(titinAnalysis(idx).lceN(idxKValid(2:end))-titinAnalysis(idx).lceN(1,1), ...
             titinAnalysis(idx).k2NS(idxKValid(2:end)),...
             '-','Color',plotConfig.lineColors.calcTitinK(idx,:),...
             'DisplayName',...
             sprintf('%s%i%s','$$k^{T}_{',idx,'}$$'));
        hold on;
        text(lceNSpan,titinAnalysis(idx).k2NS(end),...
             sprintf('%1.1f',titinAnalysis(idx).k2NS(end)),...
             'HorizontalAlignment','right',...
             'VerticalAlignment','bottom',...
             'FontSize',8);
        hold on;                
    end
    if(pubPlotOptions.useSmoothedStiffnessData==0)

        seriesColor = plotConfig.lineColors.calcTitinK(idx,:);

        plot(titinAnalysis(idx).lceN(idxKValid)-titinAnalysis(idx).lceN(1,1), ...
             titinAnalysis(idx).k2N(idxKValid),...
             '-','Color',seriesColor,...
             'DisplayName',...
             sprintf('%s: %1.1f %s','$$\ell_i$$',...
                    titinAnalysis(idx).lceN(1,1),'$$\ell_o$$'));
        hold on;
        text(lceNSpan,titinAnalysis(idx).k2N(end),...
             sprintf('%1.1f%s',titinAnalysis(idx).k2N(end),'$$f_o$$'),...
             'HorizontalAlignment','right',...
             'VerticalAlignment','bottom',...
             'FontSize',8);
        hold on;  
     
    end
end

for idx=simConfig.trials

   k2NSim = calcCentralDifferenceDataSeries(...
                benchRecordFitted.normFiberLength(:,idx),...
                benchRecordFitted.normDistalTitinForce(:,idx));

   k2NSim(isnan(k2NSim)) = 0;
   k2NSim(isinf(k2NSim)) = 0;
   

   plot( benchRecordFitted.normFiberLength(:,idx)...
          -benchRecordFitted.normFiberLength(1,idx), ...
         k2NSim,...
         '-','Color',plotConfig.lineColors.simTitinK(idx,:),...
         'DisplayName',sprintf('%s%i%s','Sim: $$k^{T*}_',idx,'$$'),...
         'HandleVisibility','on');

   hold on;            

   lceNSpan = benchRecordFitted.normFiberLength(end,idx)...
             -benchRecordFitted.normFiberLength(1,idx);

   text(lceNSpan,...
        k2NSim(end),...
         sprintf('%1.1f%s',k2NSim(end),'$$k_o$$'),...
         'HorizontalAlignment','left',...
         'VerticalAlignment','top',...
         'FontSize',8);
   hold on;

end    

box off;
xlim([0,lceNSpan]);
ylim([0,pubPlotOptions.kceNMax]);

xlabel('Norm. Length Change $$(\ell/\ell_o)-\ell_i$$');
ylabel('Norm. Slope $$(f/f_o)/(\ell/\ell_o)$$');
title('Estimated titin stiffness-length relation');  
if(idx==3)
    legend('Location','NorthWest');
end


