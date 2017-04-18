% CREATE TxRx FILE

% INPUTS:
%     type : Cell array where each entry {i} is type 'TRX_EY_DIPOLE'
%       Tx : Cell array where Tx{i} = [x z m]
%       Rx : Cell array where each Rx{i} is an N X 3 array with unique
%            [xj zj] observation locations.
%  t_LINES : Time vector for all observed times in seconds (TRX_LINES) or []
%   t_LOOP : Time vector for all observed times in seconds (TRX_LOOP) or []
% t_DIPOLE : Time vector for all observed times in seconds (TRX_DIPOLE) or []
% Filename : File name string for output file
%     Flag : String for IGNORE data (I like NaN)



%=================================================================
%       TRANSMITTER PROPERTIES
%=================================================================

type = {'TRX_EY_DIPOLE','TRX_EY_DIPOLE','TRX_EY_DIPOLE'} ;

% [X-location Z-location Moment(Ids)]
Tx{1}  = [-10 -1 1] ;
Tx{2}  = [-10 -2 1] ;
Tx{3}  = [-10 -3 1] ;
Tx{4}  = [-10 -4 1] ;
Tx{5}  = [-10 -5 1] ;
Tx{6}  = [-10 -6 1] ;
Tx{7}  = [-10 -7 1] ;
Tx{8}  = [-10 -8 1] ;
Tx{9}  = [-10 -9 1] ;
Tx{10} = [-10 -10 1] ;

%=================================================================
%       RECEIVER PROPERTIES
%=================================================================

[X,Z] = ndgrid(10,linspace(-1,-10,10)) ; 
Rx{1} = [X(:) Z(:)] ;
Rx{2} = [X(:) Z(:)] ;
Rx{3} = [X(:) Z(:)] ;
Rx{4} = [X(:) Z(:)] ;
Rx{5} = [X(:) Z(:)] ;
Rx{6} = [X(:) Z(:)] ;
Rx{7} = [X(:) Z(:)] ;
Rx{8} = [X(:) Z(:)] ;
Rx{9} = [X(:) Z(:)] ;
Rx{10} = [X(:) Z(:)] ;

%=================================================================
%       FREQUENCY PROPERTIES
%=================================================================

% FREQUENCY VECTOR
f = [1e8 2e8 4e8]' ;

% FILENAME
Filename = 'TxRxFileGPR2D.txt' ;

% FLAG VALUE
Flag = NaN ;












