function make_report_figures(selection)
%MAKE_REPORT_FIGURES Recreate native figures from saved MATLAB calculations.
% No optimization or acoustic recalculation is performed here.
% Pass 'Figure_6_1' etc. to regenerate only one figure.
if nargin<1, selection="all"; end
root=fileparts(mfilename('fullpath'));
saved=load(fullfile(root,'data','figure_data.mat'),'figureData');
figureData=saved.figureData;
figIDs=["Figure_6_1","Figure_6_2","Figure_6_3","Figure_6_4","Figure_6_5", ...
    "Figure_S1","Figure_S2","Figure_S3"];
palette=[0 .447 .698; .835 .369 0; .0 .62 .45; .8 .475 .655; .337 .706 .914; .902 .624 0];
baselineColour=palette(1,:); finalColour=palette(2,:);
styles={'-','--','-.',':','-','--'};
for name=figIDs
    if selection~="all" && selection~=name, continue; end
    destination=fullfile(root,'figures'); if ~isfolder(destination), mkdir(destination); end
    switch name
        case "Figure_6_1"
            data=figureData.parametric;
            fig=new_figure(17.8,20); layout=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
            symbols=["h","r","e"]; g=data.baseline.geometry_mm;
            for v=1:3
                ax=nexttile(layout); hold(ax,'on'); band(ax,[0 1]); scan=data.scans(v); handles=[]; labels=strings(0);
                for j=scan.representative
                    isBase=j==scan.baseline_index;
                    if isBase, colour=baselineColour; width=1.9; style='-'; else
                        index=find(scan.representative==j);
                        colour=palette(mod(index-1,5)+2,:); width=1.4; style=styles{mod(index-1,6)+1};
                    end
                    handles(end+1)=plot(ax,data.frequency_Hz,scan.curves(:,j),'Color',colour,'LineStyle',style,'LineWidth',width); %#ok<AGROW>
                    label=sprintf('%s = %.4g mm',symbols(v),scan.values_mm(j));
                    if isBase, label=[label ' (baseline)']; end
                    labels(end+1)=string(label); %#ok<AGROW>
                end
                other=setdiff(1:3,v); title(ax,sprintf('(%c) %s = %.5g mm; %s = %.5g mm','a'+v-1, ...
                    symbols(other(1)),g(other(1)),symbols(other(2)),g(other(2))),'FontWeight','normal');
                format_axis(ax,'Frequency (Hz)','Absorption coefficient',[200 2500],[0 1]);
                legend(ax,handles,labels,'Location','eastoutside','Box','off','FontSize',9);
            end
        case "Figure_6_2"
            map=figureData.depth_extension_map;
            data.out=figureData.baseline;
            h=unique(map.cavity_depth_mm); e=unique(map.neck_extension_mm);
            he=edges(h); ee=edges(e);
            fig=new_figure(17.8,11.5); ax=axes(fig); hold(ax,'on');
            for j=1:numel(map.J)
                iH=find(h==map.cavity_depth_mm(j)); iE=find(e==map.neck_extension_mm(j));
                x=[he(iH) he(iH+1) he(iH+1) he(iH)]; y=[ee(iE) ee(iE) ee(iE+1) ee(iE+1)];
                if map.is_feasible(j)
                    patch(ax,x,y,map.J(j),'EdgeColor','none');
                else
                    patch(ax,x,y,[.78 .78 .78],'EdgeColor',[.25 .25 .25]);
                    plot(ax,x([1 3]),y([1 3]),'k-','LineWidth',.7);
                end
            end
            colormap(ax,parula(128)); cb=colorbar(ax); cb.Label.String='Mean absorption, 800-1200 Hz';
            g=data.out.geometry_mm;
            marker=plot(ax,g(1),g(3),'ko','MarkerFaceColor','w','MarkerSize',7,'LineWidth',1.3,'Clipping','off');
            legend(ax,marker,'Single-type baseline','Location','southoutside','Box','off');
            title(ax,sprintf('Fixed neck radius r = %.5g mm',g(2)),'FontWeight','normal');
            format_axis(ax,'Cavity depth, h (mm)','Neck extension, e (mm)',[h(1) h(end)],[e(1) e(end)]);
        case "Figure_6_3"
            a.out=figureData.baseline; b.out=figureData.final;
            fig=new_figure(17.8,15); layout=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
            for j=1:2
                ax=nexttile(layout); hold(ax,'on'); band(ax,[0 1]);
                p(1)=plot(ax,a.out.frequency_Hz,a.out.response.absorption,'Color',baselineColour,'LineWidth',1.6);
                p(2)=plot(ax,b.out.frequency_Hz,b.out.response.absorption,'Color',finalColour,'LineWidth',1.7);
                if j==1, xr=[200 2500]; else, xr=[800 1200]; end
                format_axis(ax,'Frequency (Hz)','Absorption coefficient',xr,[0 1]);
                title(ax,sprintf('(%c)','a'+j-1),'FontWeight','normal');
                legend(ax,p,{sprintf('Baseline, J = %.6f',a.out.J),sprintf('Final mixed design, J = %.6f',b.out.J)}, ...
                    'Location','southoutside','Orientation','horizontal','Box','off');
            end
        case "Figure_6_4"
            a.out=figureData.baseline; b.out=figureData.final;
            fig=new_figure(17.8,15); layout=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
            for j=1:2
                ax=nexttile(layout); hold(ax,'on');
                if j==1, transform=@real; target=1; yl='Re(Y)'; else, transform=@imag; target=0; yl='Im(Y)'; end
                p(1)=plot(ax,a.out.frequency_Hz,transform(a.out.response.normalized_admittance),'Color',baselineColour,'LineWidth',1.6);
                p(2)=plot(ax,b.out.frequency_Hz,transform(b.out.response.normalized_admittance),'Color',finalColour,'LineWidth',1.6);
                yline(ax,target,':','Color',[.3 .3 .3],'LineWidth',1,'HandleVisibility','off');
                format_axis(ax,'Frequency (Hz)',['Normalized admittance, ' yl],[700 1300],[]);
                title(ax,sprintf('(%c)','a'+j-1),'FontWeight','normal');
                legend(ax,p,{'Baseline','Final mixed design'},'Location','southoutside','Orientation','horizontal','Box','off');
            end
        case "Figure_6_5"
            a.out=figureData.final; out=a.out; K=size(out.geometry_mm,1);
            fig=new_figure(17.8,20); layout=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
            for j=1:3
                ax=nexttile(layout); hold(ax,'on'); p=[]; labels=strings(0);
                for k=1:K
                    if j==1, value=real(out.Y_type(:,k)); yl='Type contribution, Re(Y_k)'; elseif j==2
                        value=imag(out.Y_type(:,k)); yl='Type contribution, Im(Y_k)'; else
                        value=out.eta(:,k); yl='Dissipated incident-energy fraction';
                    end
                    p(end+1)=plot(ax,out.frequency_Hz,value,'Color',palette(k+2,:),'LineStyle',styles{k},'LineWidth',1.6); %#ok<AGROW>
                    labels(end+1)=sprintf('Type %d, n = %d',k,out.counts(k)); %#ok<AGROW>
                end
                if j==3
                    p(end+1)=plot(ax,out.frequency_Hz,out.response.absorption,'Color',finalColour,'LineStyle','--','LineWidth',1.3);
                    labels(end+1)="Final absorption";
                end
                yline(ax,0,':','Color',[.5 .5 .5],'HandleVisibility','off');
                format_axis(ax,'Frequency (Hz)',yl,[800 1200],[]);
                title(ax,sprintf('(%c)','a'+j-1),'FontWeight','normal');
                legend(ax,p,labels,'Location','eastoutside','Box','off');
            end
        case "Figure_S1"
            t=figureData.one_factor_summary; t.variable=string(t.variable);
            fig=new_figure(17.8,22); layout=tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');
            variables=["cavity_depth_mm","neck_radius_mm","neck_extension_mm"];
            xlabels={'Cavity depth, h (mm)','Neck radius, r (mm)','Extension, e (mm)'};
            fields={'J','peak_frequency_Hz','peak_absorption'};
            yl={'Mean absorption','Peak frequency (Hz)','Peak absorption'};
            for metric=1:3
                for v=1:3
                    ax=nexttile(layout); rows=t.variable==variables(v) & t.is_feasible;
                    plot(ax,t.value_mm(rows),t.(fields{metric})(rows),'-o','Color',baselineColour,'MarkerSize',3,'LineWidth',1.3);
                    format_axis(ax,xlabels{v},yl{metric},[],[]);
                    title(ax,sprintf('(%c)','a'+3*(metric-1)+v-1),'FontWeight','normal');
                    if metric~=2, ylim(ax,[0 1]); end
                end
            end
        case "Figure_S2"
            t=figureData.search_runs;
            fig=new_figure(17.8,15); layout=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
            ax=nexttile(layout); hold(ax,'on'); p=[];
            searchColours=palette([1 4 2],:);
            for K=1:3
                values=figureData.search_best{K};
                p(end+1)=plot(ax,1:numel(values),values,'Color',searchColours(K,:),'LineWidth',1.5,'LineStyle',styles{K}); %#ok<AGROW>
            end
            format_axis(ax,'Geometry evaluations (local searches)','Best-so-far mean absorption',[],[]);
            legend(ax,p,{'K = 1','K = 2','K = 3'},'Location','southoutside','Orientation','horizontal','Box','off');
            title(ax,'(a) Search history','FontWeight','normal');
            ax=nexttile(layout); hold(ax,'on');
            for K=1:3
                ix=find(t.K_requested==K);
                scatter(ax,ix,t.final_J(ix),28,searchColours(K,:),'filled');
            end
            format_axis(ax,'Completed local-search run','Final mean absorption',[],[]);
            title(ax,'(b) Completed start results','FontWeight','normal');
        case "Figure_S3"
            a.out=figureData.final; out=a.out;
            fig=new_figure(17.8,10.5); ax=axes(fig); hold(ax,'on'); band(ax,[0 1]); p=[]; labels=strings(0);
            for k=1:size(out.geometry_mm,1)
                p(end+1)=plot(ax,out.frequency_Hz,figureData.counterfactual_absorption(:,k),'Color',palette(k+2,:),'LineStyle',styles{k},'LineWidth',1.5); %#ok<AGROW>
                labels(end+1)=sprintf('16 identical cells of Type %d',k); %#ok<AGROW>
            end
            p(end+1)=plot(ax,out.frequency_Hz,out.response.absorption,'Color',finalColour,'LineWidth',1.8);
            labels(end+1)="Final mixed panel";
            format_axis(ax,'Frequency (Hz)','Absorption coefficient',[200 2500],[0 1]);
            legend(ax,p,labels,'Location','southoutside','Box','off','NumColumns',2);
    end
    drawnow;
    exportgraphics(fig,fullfile(destination,name+'.png'),'Resolution',600);
    close(fig);
    fprintf('Exported %s\n',name);
end
end

function fig=new_figure(width,height)
fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 width height], ...
    'DefaultAxesFontName','Arial','DefaultAxesFontSize',9,'DefaultTextFontName','Arial','DefaultTextFontSize',9);
end

function format_axis(ax,xlabelText,ylabelText,xrange,yrange)
set(ax,'FontName','Arial','FontSize',9,'LineWidth',.8,'Box','off','TickDir','out','Layer','top');
xlabel(ax,xlabelText); ylabel(ax,ylabelText);
if ~isempty(xrange), xlim(ax,xrange); end
if ~isempty(yrange), ylim(ax,yrange); end
end

function band(ax,ylimits)
patch(ax,[800 1200 1200 800],[ylimits(1) ylimits(1) ylimits(2) ylimits(2)], ...
    [.92 .94 .96],'EdgeColor','none','HandleVisibility','off');
end

function output=edges(values)
values=values(:)'; output=[values(1) (values(1:end-1)+values(2:end))/2 values(end)];
end
