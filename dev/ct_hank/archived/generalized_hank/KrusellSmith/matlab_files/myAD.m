function x = myAD(x, varargin)
% Modified Package by SeHyoun Ahn Copyright (C) 2015 (SeHyoun Ahn) based on
% /*
%  * myAD and myA2D - Automatic Differentiation of 1st
%  * Copyright (C) 2006 Martin Fink. (martinfink "at" gmx.at)
%  *
%  * This library is free software; you can redistribute it and/or
%  * modify it under the terms of the GNU General Public License
%  * as published by the Free Software Foundation; either version 2.1
%  * of the License, or (at your option) any later version.
%  *
%  * This library is distributed in the hope that it will be useful,
%  * but WITHOUT ANY WARRANTY; without even the implied warranty of
%  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
%  * General Public License for more details.
%  *
%  * You should have received a copy of the GNU General Public License
%  * along with this library; if not, visit
%  * http://www.gnu.org/licenses/gpl.html or write to the Free Software
%  * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
%  */
% License. This is released with GPL 2.1.

if (isa(x, 'myAD'))
    return;
end

struct_x.values = x;
if (nargin < 2)
    struct_x.derivatives = speye(numel(x));
else
    struct_x.derivatives = varargin{1};
end

x = class(struct_x, 'myAD');