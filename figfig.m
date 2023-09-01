function  figfig (name)
% appends last figure in a PDF file, with the spesified filename.
% name='tt' 'ii.png'
if size(name,2)>3
    if strcmpi(name(end-3:end),'.png')
        a_name=[name];
    else
        if strcmpi(name(end-3:end),'.pdf')
            a_name=[name];
        else
            a_name=[name '.pdf'];
        end
    end
else
    a_name=[name '.pdf'];
end

a=dir;
fil=strcmpi({a.name},a_name);
i=1;
if isempty(find(fil,1))
    fname =   a_name;
else
    while find(fil)>0
        if strcmpi(a_name(end-3:end),'.pdf')
            fname = sprintf([a_name(1:end-4) '_%d.pdf'], i);
        elseif strcmpi(a_name(end-3:end),'.png')
            fname = sprintf([a_name(1:end-4) '_%d.png'], i);
        end
        fil=strcmpi({a.name},fname);
        i=i+1;
    end
end

% styling
fig=gcf; fig2=gca;
try
    fig2.XColor=[0 0 0 1]; fig2.YColor=[0 0 0 1]; fig2.ZColor=[0 0 0 1];
    %fig2.MinorGridColor=[0 0 0]; %  fig2.GridColor=[0 0 0];fig2.GridAlpha=1; 
    % jj = findall(fig2, 'Type', 'text');
    % jj(1).String=strjoin(jj(1).String)
end
% try
%     for ff=1:size(fig2.Children,1)
%         if contains(class(fig2.Children(ff,:)),'text','IgnoreCase',true)
%             fig2.Children(ff).Color=[0 0 0];
%         end
%     end
% end
try
    fig2.Color='none'; %[1 1 1 1];
end
try
    fig.GraphicsSmoothing='on'; fig.InvertHardcopy='off';
end
try
fig2.FontColor=[0 0 0 1];
end

fig2.FontName='Helvetica' ; % fig2.FontSize=11;
fig.PaperUnits='centimeters'; fig.Units='centimeters';
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% change normalized units to data units when possible
% normalized units get messed up when upsampling - taken from: 
% https://github.com/djoshea/matlab-save-figure/blob/master/saveFigure.m
axh1 = findall(fig2, 'Type', 'Axes');
objNormalizedUnits1 = findall(axh1, 'Units', 'normalized', '-not', 'Type', 'Axes');
set(objNormalizedUnits1, 'Units', 'data');
axh = findall(fig, 'Type', 'Axes');
objNormalizedUnits = findall(axh, 'Units', 'normalized', '-not', 'Type', 'Axes');
set(objNormalizedUnits, 'Units', 'data');


%saving the Figure
if strcmpi(fname(end-3:end),'.pdf')
    fig.Color='none';
    %print(fig,'-dpdf','-painters','-fillpage','-r600',fname); %'-bestfit',
    exportgraphics(fig,fname,'Append',true,'BackgroundColor','none','ContentType','vector')
    %,'Append',true,'BackgroundColor','none','ContentType','vector'
    % print(fig,'-dpsc','-painters','-r600','-bestfit','-append','yoyo.ps');
elseif strcmpi(fname(end-3:end),'.png')
    try
    fig.Color=[1 1 1 1];
    end
    %print(fig,'-dpng','-r600',fname);
    exportgraphics(fig,fname,'BackgroundColor','none','ContentType','image')
end
%print(gcf,'-depsc2','-painters','-r600','12.eps');

% if spesified filename exist the append. 
% needs append_pdf function from export_fig on github https://github.com/altmany/export_fig
% export_fig require to install https://www.ghostscript.com/download/gsdnld.html 
if i>1 && strcmpi(fname(end-3:end),'.pdf')
    try
        append_pdfs(a_name, fname)
        delete(fname)
    catch
        warning(['PDF file appending is impossible - PDF file may be opened elsewhere'  newline 'file was saved with index name'])
   end
end
try fig.Color=[1 1 1 1]; end 
try fig2.Color=[1 1 1 1]; end


end

%% JUNK
% % trick Matlab into rendering everything at higher resolution
% jc = findobjinternal(fig, '-isa', 'matlab.graphics.primitive.canvas.JavaCanvas', '-depth', 1);

% if isa(jc, 'matlab.graphics.GraphicsPlaceholder')
%     %warning('Could not determine screen DPI: JavaCanvas not found');
%     renderDPI = 600;
% else
%     origDPI = jc.ScreenPixelsPerInch;
%     if origDPI < 120
%         renderDPI = 120;  % origDPI * p.Results.upsample;
%         if ~isempty(jc)
%             %             jc.OpenGL = 'off';
%             jc.ScreenPixelsPerInch = 100;
%             if exist('AutoAxis', 'class')
%                 AutoAxis.updateFigure();
%             end
%         end
%     else
%         renderDPI = 600; %origDPI;
%     end
% end




%append_pdfs(outputFilename, inputFilename1, inputFilename2, ...)


%
% set(fig,'Units','inches');
% screenposition = get(fig,'Position');
% set(fig,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print -dpdf -painters -append name
%
% fig.PaperPositionMode='auto';
% print(name,'-dpdf','-painters','-r600','-bestfit',strcat(plot_path,plot_name))
% print(name,'-dpdf','-fillpage','-append')

