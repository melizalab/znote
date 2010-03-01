function [] = zedit(signalfile, varargin)
%
% ZEDIT Front-end GUI for editing the feature composition
% of starling motifs.
%
% ZEDIT(signal) Starts the GUI with the input <signal>
% using default parameters
%
% ZEDIT(signal, params) Starts the GUI with an optional param
% structure. See ZEDIT_PARAMS for the available parameters and their
% default values.
%
% Copyright C Daniel Meliza, Z Chi 2010.  Licensed for use under Creative
% Commons Attribution-Noncommercial-Share Alike 3.0 United States
% License (http://creativecommons.org/licenses/by-nc-sa/3.0/us/).

if nargin < 1
  error(['Need an input file']);
end

if exist(signalfile,'file') ~= 2
     error(['Input file ' signalfile ' does not exist']);
end

varargin{:}
[p pcmroot pcmext] = fileparts(signalfile);
params = zedit_params;

% open the GUI
figtag = [mfilename '_' pcmroot];
fig = findobj('tag',figtag);
if isempty(fig)
     fig = open([mfilename '.fig']);
     set(fig,'tag',figtag,'name',signalfile);%,'renderer','opengl');
     handles = guihandles(fig);
     h = findobj(fig,'style','pushbutton');
     set(h,'Callback',@cb_btn);
     h = findobj(fig,'style','listbox');
     set(h,'Callback',@cb_list);
     set(handles.thresh,'Callback',@update_thresh);
else
     figure(fig);
     handles = guidata(fig);
end

% store the params and handles
params.pcmroot = pcmroot;
params.pcmext  = pcmext;
handles.params = params;
guidata(fig, handles);

% update callbacks

% populate fields
update_fields(fig);
% calculate the starting PSD
update_psd(fig);



%%%%%%%%%%%%%%%% SUBFUNCTIONS

function [] = create_label(fig)
% computes connected components based on the parameter structure
params  = guivalue(fig,'params');
[L,PSD] = zlabel([params.pcmroot params.pcmext], params);
update_psd(fig, PSD);
if ~isempty(L)
  % store the label matrix in the guidata
  labels  = guivalue(fig,'labels');
  %setname = sprintf('lbl%2.1f_%d_%d_%d', params.thresh, params.df, ...
  %                  params.dt, params.min_size);
  setname = sprintf('lbl%d_%d_%d', params.df, params.dt, params.min_size);
  labels.(setname) = L;
  guivalue(fig,'labels',labels);
  % update the label list
  update_labelset(fig, setname);
end



function [] = update_alpha(fig)
% based on the selected labelset, and selected features, sets the
% alpha property of the spectrogram to values that result in only
% the selected features being visible against a light
% background. If there are no no features, or no spectrogram, we
% don't do anything.
BACKGND_ALPHA = 200;

Z       = 1;
h       = guivalue(fig,'list_labelset');
labelset= getselected(h);
if ~isempty(labelset)
     h  = guivalue(fig,'list_features');
     featureset = str2num(char(getselected(h)));
     if ~isempty(featureset)
          % load the labels
          lblmat = getfield(guivalue(fig,'labels'), labelset);
          Z      = zeros(size(lblmat));
          Z(ismember(lblmat, featureset)) = 1;
          Z(lblmat==0) = BACKGND_ALPHA;
     end
end
h_spec   = findobj(gca,'tag','spec');
set(h_spec(1),'AlphaData',Z);

          


function [] = update_labelset(fig, selected)
% updates the labelset list with the current contents of the labels
% guidata
h    = guivalue(fig,'list_labelset');
lbl  = guivalue(fig,'labels');
if isempty(lbl) || isempty(fieldnames(lbl))
     set(h,'String',{},'Value',0);
     selected = [];
else
     lset = fieldnames(lbl);
     if nargin < 2
          selected = getselected(h);
     end
     set(h,'String',lset);
     ind = strmatch(selected,lset);
     if ~isempty(ind)
          set(h,'Value',ind);
     end
end
% and load the labelset into the feature list
update_featurelist(fig, selected);




function [] = update_featurelist(fig, setname, selected_features)
% updates the feature list with the contents of the currently
% selected labelset
if nargin < 2
     h      = guivalue(fig,'list_labelset');
     setname= getselected(h);
end

h        = guivalue(fig,'list_features');
if isempty(setname)
     set(h,'String',{},'Value',0);
else
     % need to check that the labelfile is valid
     labels   = getfield(guivalue(fig,'labels'),setname);
     h_spec   = findobj(fig,'tag','spec');
     if ~isempty(h_spec) && all(size(get(h_spec,'CData'))==size(labels))
          features = strtrim(cellstr(num2str(unique(labels))));
          features = features(2:end);
          set(h,'String',features,'Value',1:length(features));
          if nargin >= 3 
               setselected(h,selected_features);
          end
     else
          set(h,'String','Wrong dimensions','Value',1);
     end
end
% update the alpha of the spectrogram based on selected features
update_alpha(fig);




function [] = update_psd(fig, spectro)
% computes the time-frequency PSD. Can pass it a already-computed
% spectro, which is useful when computing labels.
params  = guivalue(fig,'params');
if nargin < 2
     [spectro,T,F] = zpsd([params.pcmroot params.pcmext], params);
else
  h_spec = findobj(fig,'tag','spec');
  T = get(h_spec,'XData');
  F = get(h_spec,'YData');
end

% plot it to the axes
ax = guivalue(fig,'axes');
delete(kids(ax));

imagesc(F,T,spectro,'tag','spec','parent',ax);
xlabel('Time (ms)');ylabel('F (Hz)');
axis xy, axis tight, box on;
% update subsidiary plots
update_thresh(fig);



function [] = update_fields(fig)
% cycles through the param struct and sets any objects it can find
% to the correct values
params = guivalue(fig,'params');
for i = fieldnames(params)'
     fn = i{1};
     h = findobj(fig,'tag',fn);
     if ~isempty(h)
          if isnumeric(params.(fn))
               set(h,'String',num2str(params.(fn)));
          else
               set(h,'String',params.(fn));
          end
     end
end



function [] = update_params(fig)
% has the opposite effect of update_fields, sets values in the
% params from the fields
params = guivalue(fig,'params');
for i = fieldnames(params)'
     fn = i{1};
     h = findobj(fig,'tag',fn);
     if ~isempty(h)
          val = get(h,'String');
          if isnumeric(params.(fn))
               params.(fn) = str2num(val);
          else
               params.(fn) = val;
          end
     end     
end
guivalue(fig,'params',params);



function [] = update_thresh(fig, event)
% updates the plots contingent on the threshhold value (the
% colorbar, its marker, and the contour plot)
% note, this can be used as a callback
h_ax   = guivalue(fig,'axes');
h_old  = findobj(h_ax, 'tag', 'contour') ;
delete(h_old);

h_spec = findobj(h_ax, 'tag', 'spec');

T = get(h_spec,'xdata');
F = get(h_spec,'ydata');
val = str2num(get(guivalue(fig,'thresh'),'String'));
hold on,contour(T,F,get(h_spec,'cdata'),[val val], 'k', 'tag', 'contour');

h_cb   = colorbar;
line(get(h_cb,'xlim'),[val val],'color','black','parent',h_cb);
set(h_cb,'buttondownfcn',@cb_thresh_slider);




%%%%%%% CALLBACKS

function [] = cb_thresh_slider(obj, event)
% handles click events on the colorbar by moving the threshhold to
% the new location
click_val = get(obj, 'currentpoint');
val       = click_val(1,2);
% store the threshold in the edit control
h         = guivalue(obj, 'thresh');
set(h,'String',num2str(val));
% and update the contour map
update_thresh(gcbf);

function [] = cb_btn(obj, event)
% button switchyard
tag = lower(get(obj, 'tag'));

% a lot of the buttons need these values
h        = guivalue(obj,'list_labelset');
labelset = getselected(h);     
h        = guivalue(obj,'list_features');
selected = str2num(getselected(h));

if strcmp(tag,'btn_calc_psd')
     % update the params from the fields
     update_params(gcbf);
     % update the psd
     update_psd(gcbf);

elseif strcmp(tag,'btn_get_labels')
     % update the params
     update_params(gcbf);
     % calculate the label
     create_label(gcbf);

elseif strcmp(tag,'btn_load_labelset')
     % reads a label file from disk
     pos = get(gcbf,'position');
     [fn pn] = uigetfile('*.bin','Location',pos(1:2));
     if (fn ~= 0)
          lblmat = bimatrix(fullfile(pn,fn),'int') + 1;
          [p fn e] = fileparts(fn);
          labels = guivalue(obj,'labels');
          labels.(fn) = lblmat;
          guivalue(obj,'labels',labels);
          update_labelset(gcbf,fn);
     end

elseif strcmp(tag,'btn_save_labelset')
     if isempty(labelset)
          errordlg('No label file to save!');
          return
     end
     pos = get(gcbf,'position');
     [fn pn] = uiputfile('*.bin','Location',pos(1:2));
     if (fn ~= 0)
          lblmat = getfield(guivalue(obj,'labels'),labelset) - 1;
          bomatrix(lblmat,fullfile(pn,fn),'int');
          fprintf('Wrote labelset %s to file %s\n', labelset, fn);
     end
     
elseif strcmp(tag,'btn_drop_labelset')
     if ~isempty(labelset)
          labels = guivalue(obj,'labels');
          labels = rmfield(labels, labelset);
          guivalue(obj,'labels',labels);
          update_labelset(gcbf);
     end
     
%elseif strcmp(tag, 'btn_merge_labelset')
     

elseif strcmp(tag,'btn_lasso_feat')
     % allow the user to lasso an area of the spectrogram. The
     % intersection of the lassoed area and the existing features
     % (under selection) will be created as a new feature
     if isempty(labelset) || isempty(selected)
          return
     end
     
     h_ax   = guivalue(obj,'axes');
     h_spec   = findobj(h_ax,'tag','spec');
     
     set(obj,'enable','off');
     [IN, h_line] = lasso(h_spec);
     set(obj,'enable','on');
     
     labels = guivalue(obj,'labels');
     isect  = ismember(labels.(labelset),selected) & IN;
     new_ind = max(unique(labels.(labelset)))+1;
     labels.(labelset)(isect) = new_ind;
     guivalue(obj,'labels',labels);
     fprintf('Generated new feature %d\n', new_ind);
     delete(h_line);
     
     update_featurelist(gcbf, labelset, new_ind);
     
elseif strcmp(tag,'btn_merge_feat')
     % pass the label file to zmerge
     if isempty(labelset) || isempty(selected)
          return
     end     
     
     labels   = guivalue(obj,'labels');
     labels.(labelset)   = zmerge(labels.(labelset), selected);
     guivalue(obj,'labels',labels);
     fprintf('Merged features [%s] into %d\n', num2str(selected'), selected(1));
     
     update_featurelist(gcbf, labelset, selected(1));
     
elseif strcmp(tag, 'btn_drop_feat')
     % pass the label file to zmerge
     if isempty(labelset) || isempty(selected)
          return
     end
     
     labels   = guivalue(obj,'labels');
     labels.(labelset)   = zmerge(labels.(labelset), selected,'delete');
     guivalue(obj,'labels',labels);
     fprintf('Deleted features [%s]\n', num2str(selected'));
     
     update_featurelist(gcbf, labelset, num2str(selected(1)+1));
     
elseif strcmp(tag,'btn_sort_feat')
     labels   = guivalue(obj,'labels');
     if isempty(labelset) || isempty(labels)
          return
     end
     labels.(labelset)   = zmerge(labels.(labelset),[], ...
                                      'renorm');
     guivalue(obj,'labels',labels);
     update_featurelist(gcbf, labelset);
     fprintf('Resorted features\n');
     
else
     fprintf('Callback not implemented for object %s\n', tag);
end


function [] = cb_list(obj, event)
% switchyard for list callbacks
tag = lower(get(obj, 'tag'));

if strcmp(tag,'list_labelset')
     update_featurelist(gcbf);
elseif strcmp(tag,'list_features')
     update_alpha(gcbf);
else
     fprintf('Callback not implemented for %s\n', tag);
end


