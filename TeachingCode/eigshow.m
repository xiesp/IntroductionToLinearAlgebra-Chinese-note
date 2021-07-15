function eigshow(arg)
%EIGSHOW Graphical demonstration of eigenvalues and singular values.
%
%   EIGSHOW presents a graphical experiment showing the effect on the
%   the unit circle of the mapping induced by various 2-by-2 matrices.
%   A pushbutton allows the choice of "eig" mode or "svd" mode.
%
%   In eig mode, the mouse can be used to move the vector x around the
%   unit circle.  The resulting trajectory of A*x is plotted.  The object
%   is to find vectors x so that A*x is parallel to x.  Each such x is an
%   eigenvector of A.  The length of A*x is the corresponding eigenvalue.
%
%   In svd mode, the mouse moves two perpendicular unit vectors, x and y.
%   The resulting A*x and A*y are plotted.  When A*x is perpendicular to
%   A*y, then x and y are right singular vectors, A*x and A*y are
%   multiples of left singular vectors, and the lengths of A*x and A*y
%   are the corresponding singular values.
%
%   The figure title is a menu of selected matrices, including some
%   with fewer than two real eigenvectors.  EIGSHOW(A) inserts A,
%   which must be 2-by-2, in the menu.
%
%   Here are some questions to consider:
%      Which matrices are singular?
%      Which matrices have complex eigenvalues?
%      Which matrices have double eigenvalues?
%      Which matrices have eigenvalues equal to singular values?
%      Which matrices have nondiagonal Jordan canonical forms?

%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 1997/11/21 23:25:37 $

if nargin == 0;
   initialize
elseif arg == 0
   action
elseif arg < 0
   setmode(arg)
else
   initialize(arg);
end

%------------------

function initialize(arg)

if nargin == 0
   arg = 6;
end

if isequal(get(gcf,'tag'),'eigshow');
   h = get(gcf,'userdata');
   mats = h.mats;
else
   set(gcf,'numbertitle','off','menubar','none')
   h.svd = 0;
   mats = {
      '[5/4 0; 0 3/4]'
      '[5/4 0; 0 -3/4]'
      '[1 0; 0 1]'
      '[0 1; 1 0]'
      '[0 1; -1 0]'
      '[1 -1;1 -1]'
      '[1 1;0 1]'
      '[1 3; 4 2]/4'
      '[1 3; 2 4]/4'
      '[3 1; 4 2]/4'
      '[3 1; -2 4]/4'
      '[2 4; 2 4]/4'
      '[2 4; -1 -2]/4'
      '[6 4; -1 2]/4'
      'randn(2,2)'
      };
end

if all(size(arg)==1)
   if (arg < length(mats))
      mindex = arg;
      A = eval(mats{mindex});
   else
      A = randn(2,2);
      S = ['[' sprintf('%4.2f %4.2f; %4.2f %4.2f',A) ']'];
      mindex = length(mats);
      mats = [mats(1:mindex-1); {S}; mats(mindex)];
   end
else
   A = arg;
   if isstr(A)
      S = A;
      A = eval(A);
   else
      S = ['[' sprintf('%4.2f %4.2f; %4.2f %4.2f',A) ']'];
   end
   if any(size(A) ~= 2)
      error('Matrix must be 2-by-2')
   end
   mats = [{S};  mats];
   mindex = 1;
end

clf
if h.svd, t = 'svd / (eig)';
    else, t = 'eig / (svd)';
end
uicontrol(...
   'style','pushbutton', ...
   'units','normalized', ...
   'position',[.86 .60 .12 .06], ...
   'string',t, ...
   'value',h.svd, ...
   'callback','eigshow(-1)');
uicontrol(...
   'style','pushbutton', ...
   'units','normalized', ...
   'position',[.86 .50 .12 .06], ...
   'string','help', ...
   'callback','helpwin eigshow')
uicontrol(...
   'style','pushbutton', ...
   'units','normalized', ...
   'position',[.86 .40 .12 .06], ...
   'string','close', ...
   'callback','close(gcf)')
uicontrol(...
   'style','popup', ...
   'units','normalized', ...
   'position',[.28 .92 .48 .08], ...
   'string',mats, ...
   'tag','mats', ...
   'fontname','courier', ...
   'fontweight','bold', ...
   'fontsize',14, ...
   'value',mindex, ...
   'callback','eigshow(get(gco,''value''))');

s = 1.1*max(1,norm(A));
axis([-s s -s s])
axis square
xcolor = [0 .6 0];
Axcolor = [0 0 .8];
h.A = A;
h.mats = mats;
h.x = initv([1 0]','x',xcolor);
h.Ax = initv(A(:,1),'Ax',Axcolor);
if h.svd
   h.y = initv([0 1]','y',xcolor);
   h.Ay = initv(A(:,2),'Ay',Axcolor);
   xlabel('Make A*x perpendicular to A*y','fontweight','bold')
   set(gcf,'name','svdshow')
else
   xlabel('Make A*x parallel to x','fontweight','bold')
   set(gcf,'name','eigshow')
end
set(gcf,'tag','eigshow', ...
   'userdata',h, ...
   'windowbuttondownfcn', ...
      'eigshow(0); set(gcf,''windowbuttonmotionfcn'',''eigshow(0)'')', ...
   'windowbuttonupfcn', ...
      'set(gcf,''windowbuttonmotionfcn'','''')')

%------------------

function h = initv(v,t,color)
h.mark = line(v(1),v(2),'marker','.','erase','none','color',color);
h.line = line([0 v(1)],[0 v(2)],'erase','xor','color',color);
h.text = text(v(1)/2,v(2)/2,t,'fontsize',12,'erase','xor','color',color);

%------------------

function action
h = get(gcf,'userdata');
pt = get(gca,'currentpoint');
x = pt(1,1:2)';
x = x/norm(x);
movev(h.x,x);
A = h.A;
movev(h.Ax,A*x);
if h.svd
   y = [-x(2); x(1)];
   movev(h.y,y);
   movev(h.Ay,A*y);
end

%------------------

function movev(h,v)
set(h.mark,'xdata',v(1),'ydata',v(2));
set(h.line,'xdata',[0 v(1)],'ydata',[0 v(2)]);
set(h.text,'pos',v/2);

%------------------

function setmode(arg)
h = get(gcf,'userdata');
h.svd = ~h.svd;
set(gcf,'userdata',h)
initialize(get(findobj(gcf,'tag','mats'),'value'))


