function nrb2iges (nurbs, filename)
% NRB2IGES : Write a NURBS curve or surface to an IGES file.
%
% Calling Sequence:
% 
%   nrb2iges (nurbs, filename);
% 
% INPUT:
% 
%   nurbs    : NURBS curve or surface, see nrbmak.
%   filename : name of the output file.
%  
% Description:
% 
%   The data of the nurbs structure is written in a file following the IGES
%   format. For a more in-depth explanation see, for example:
%   <http://engineeronadisk.com/V2/notes_design/engineeronadisk-76.html>.
%
%    Copyright (C) 2014 Jacopo Corno
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%    This file is based on nurbs2iges.m, (C) 2006 Fu Qiang, originally 
%    released under the MIT license.

  dt = datestr (now, 'yyyy.mm.dd');
  
  dim = numel (nurbs(1).order);
  
  % START SECTION
  S{1} = '';
  S{2} = 'IGES obtained from Nurbs toolbox.';
  S{3} = 'See <http://octave.sourceforge.net/nurbs/>.';
  S{4} = '';
    
  % GLOBAL SECTION
  G{1} = '1H,';                           % Parameter Deliminator Character
  G{2} = '1H;';                           % Record Delimiter Character
  G{3} = HString ('Nurbs toolbox');       % Product ID from Sender
  G{4} = HString (filename);              % File Name
  G{5} = HString ('Octave Nurbs');        % System ID
  G{6} = HString ('nrb2iges');            % Pre-processor Version
  G{7} = '32';                            % Number of Bits for Integers (No. of bits present in the integer representation of the sending system)
  G{8} = '75';                            % Single Precision Magnitude (Maximum power of 10 which may be represented as a single precision floating point number from the sending system)
  G{9} = '6';                             % Single Precision Significance (No. of significant digits of a single precision floating point number on the sending system)
  G{10}= '75';                            % Double Precision Magnitude (Maximum power of 10 which may be represented as a double precision floating point number from the sending system)
  G{11}= '15';                            % Double Precision Significance (No. of significant digits of a double precision floating point number on the sending system)
  G{12}= HString('Nurbs from Octave');    % Product ID for Receiver
  G{13}= '1.0';                           % Model Space Scale
  G{14}= '6';                             % Unit Flag (6 = metres)
  G{15}= HString('M');                    % Units  (metres = "M")
  G{16}= '1000';                          % Maximum Number of Line Weights
  G{17}= '1.0';                           % Size of Maximum Line Width
  G{18}= HString(dt);                     % Date and Time of file generation
  G{19}= '0.000001';                      % Minimum User-intended Resolution
  G{20}= '10000.0';                       % Approximate Maximum Coordinate
  G{21}= HString('Jacopo Corno');         % Name of Author
  G{22}= HString('GSCE - TU Darmstadt');  % Author's Organization
  G{23}= '3';                             % IGES Version Number (3 = IGES version 2.0)
  G{24}= '0';                             % Drafting Standard Code (0 = no standard)
  
  % Convert section array to lines (maximum lenght 72)
  SectionG = make_section (G, 72);
  
  % DIRECTORY ENTRY SECTION
  % Each directory entry consists of two, 80 character, fixed formatted lines
  D = [];
  for ii = 1:length (nurbs)
    switch (dim)
      case 1 % NURBS curve
        D(ii).type = 126;
      case 2 % NURBS surface
        D(ii).type = 128;
      otherwise
        error ('Only curves and surfaces can be saved in IGES format.')
    end
    D(ii).id = 2*ii - 1; % odd counter (see Parameter data section)
    D(ii).p_start = 0;
    D(ii).p_count = 0;
  end
  
  % PARAMETER DATA SECTION
  % The structure is a free formatted data entry from columns 1 to 64.
  % Each line of free formatted data consists of the entity type number
  % followed by the parameter data describing the entity.
  % Columns 65 to 72 are reserved for a parameter data index which is an
  % odd number counter, right justified in the field, which begins at the
  % number 1 and progresses in odd increments for each entity entered.
  % Column 73 is reserved for the letter ‘P’ to indicate the data element
  % belongs to the parameter data section. 
  % Columns 74 to 80 are reserved for the sequence number. Each line of 
  % data corresponds to the entity type as specified in the global section.
  SectionP = {};
  for ii = 1:length (nurbs)
    P = make_section_array (nurbs(ii)); % finish one entity
    % Convert section array to lines
    SP = make_section (P, 64);
    D(ii).p_count = length (SP);
    if (ii == 1)
        D(ii).p_start = 1;
    else
        D(ii).p_start = D(ii-1).p_start + D(ii-1).p_count;
    end
    SectionP{ii} = SP;
  end

  % SAVE
  fid = fopen (filename, 'w');
  
  % Save Start Section
  for ii = 1:length (S)
    fprintf (fid, '%-72sS%7d\n', S{ii}, ii);
  end
  
  % Save Global Section
  for ii = 1:length (SectionG)
    fprintf (fid, '%-72sG%7d\n', SectionG{ii}, ii);
  end
  
  % Save Directory Entry Section
  for i = 1:length (D)
    fprintf (fid, '%8d%8d%8d%8d%8d%8d%8d%8d%8dD%7d\n', ...
             D(i).type, D(i).p_start, 0, 0 ,0, 0, 0, 0, 0, i*2-1);
    fprintf (fid, '%8d%8d%8d%8d%8d%8s%8s%8s%8dD%7d\n', ...
             D(i).type, 0, 0, D(i).p_count, 0, ' ', ' ', ' ', 0, i*2);
  end
  
  % Save Parameter Data Section
  lines_p = 0;
  for jj = 1:length (D)
    sec = SectionP{jj};
    for ii = 1:length (sec)
      lines_p = lines_p + 1;
      fprintf (fid, '%-64s %7dP%7d\n', sec{ii}, D(jj).id, lines_p);
    end
  end
  
  % Save Terminate Section
  sec_t = sprintf ('%7dS%7dG%7dD%7dP%7d', length (S), length(SectionG), 2*length(D), lines_p);
  fprintf (fid, '%-72sT%7d\n', sec_t, 1);
  
  fclose(fid);
  
end

function P = make_section_array (nurbs)
  dim = numel (nurbs.order);
  
  % in IGES the control points are stored in the format [x, y, z, w]
  % instead of [w*x, w*y, w*z, w]
  for idim = 1:3
    nurbs.coefs(idim,:) = nurbs.coefs(idim,:) ./ nurbs.coefs(4,:);
  end
  
  P = {};
  switch dim
    case 1
      % Rational B-Spline Curve Entity
      cp = nurbs.coefs;
      deg = nurbs.order - 1;
      knots = nurbs.knots;
      uspan = [0 1];
      isplanar = ~any(cp(3,:));
      P{1} = '126';                       % NURBS curve
      P{2} = int2str (size (cp, 2) - 1);  % Number of control points
      P{3} = int2str (deg);               % Degree
      P{4} = int2str (isplanar);          % Curve on xy plane
      P{5} = '0';
      P{6} = '0';
      P{7} = '0';
      index = 8;
      for ii = 1:length (knots)
        P{index} = sprintf ('%f', knots(ii));
        index = index + 1;
      end
      for ii = 1:size (cp, 2)
        P{index} = sprintf ('%f', cp(4,ii));
        index = index + 1;
      end
      for ii = 1:size (cp, 2)
        P{index} = sprintf ('%f', cp(1,ii));
        index = index + 1;
        P{index} = sprintf ('%f', cp(2,ii));
        index = index + 1;
        P{index} = sprintf ('%f', cp(3,ii));
        index = index + 1;
      end
      P{index} = sprintf ('%f', uspan(1));
      index = index +1;
      P{index} = sprintf ('%f', uspan(2));
      index = index +1;
      P{index} = '0.0';
      index = index +1;
      P{index} = '0.0';
      index = index +1;
      if isplanar
        P{index} = '1.0';
      else
        P{index} = '0.0';
      end
      index = index + 1;
      P{index} = '0';
      index = index + 1;
      P{index} = '0';
    case 2
      % Rational B-Spline Surface Entity
      cp = nurbs.coefs;
      degU = nurbs.order(1) - 1;
      degV = nurbs.order(2) - 1;
      knotsU = nurbs.knots{1};
      knotsV = nurbs.knots{2};
      uspan = [0 1];
      vspan = [0 1];
      P{1} = '128';                       % NURBS surface
      P{2} = int2str (size (cp, 2) - 1);  % Number of control points in U
      P{3} = int2str (size (cp, 3) - 1);  % Number of control points in V
      P{4} = int2str (degU);              % Degree in U
      P{5} = int2str (degV);              % Degree in V
      P{6} = '0';
      P{7} = '0';
      P{8} = '0';
      P{9} = '0';
      P{10} = '0';
      index = 11;
      for ii = 1:length (knotsU)
        P{index} = sprintf ('%f', knotsU(ii));
        index = index + 1;
      end
      for ii = 1:length (knotsV)
        P{index} = sprintf ('%f', knotsV(ii));
        index = index + 1;
      end
      for jj = 1:size (cp, 3)
        for ii = 1:size (cp, 2)
          P{index} = sprintf ('%f', cp(4,ii,jj));
          index = index + 1;
        end
      end
      for jj = 1:size (cp, 3)
        for ii = 1:size (cp, 2)
          P{index} = sprintf ('%f',cp(1,ii,jj));
          index = index + 1;
          P{index} = sprintf ('%f',cp(2,ii,jj));
          index = index + 1;
          P{index} = sprintf ('%f',cp(3,ii,jj));
          index = index + 1;
        end
      end
      P{index} = sprintf('%f',uspan(1));
      index = index +1;
      P{index} = sprintf('%f',uspan(2));
      index = index +1;
      P{index} = sprintf('%f',vspan(1));
      index = index +1;
      P{index} = sprintf('%f',vspan(2));
      index = index +1;
      P{index} = '0';
      index = index + 1;
      P{index} = '0';
    otherwise
      
  end
end

function hs = HString (str)
% HString : Convert the string STR to the Hollerith format.

  hs = sprintf ('%dH%s', length(str), str);
end

function sec = make_section (fields, linewidth)
  sec = {};
  index = 1;
  line = '';
  num = length (fields);
  for i = 1:num
    if (i < num)
      newitem = [fields{i} ','];
    else
      newitem = [fields{i} ';'];
    end
    len = length (line) + length (newitem);
    if ( len > linewidth )
      % new line
      sec{index} = line;
      index = index + 1;
      line = '';
    end
    line = [line newitem];
  end
  sec{index} = line;
end
