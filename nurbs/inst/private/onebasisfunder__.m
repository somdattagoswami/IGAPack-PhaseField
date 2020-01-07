function [N, Nder] = onebasisfunder__ (u, p, U)

%  __ONEBASISFUNDER__: Undocumented internal function
%
%   Copyright (C) 2012 Rafael Vazquez
%   This software comes with ABSOLUTELY NO WARRANTY; see the file
%   COPYING for details.  This is free software, and you are welcome
%   to distribute it under the conditions laid out in COPYING.

  N = zeros (size (u));
  Nder = zeros (size (u));

  for ii = 1:numel (u)
    if (~ any (U <= u(ii))) || (~ any (U > u(ii)))
      continue;
    elseif (p == 0)
      N(ii) = 1;
      Nder(ii) = 0;
      continue;
    else
      ln = u(ii) - U(1);
      ld = U(end-1) - U(1);
      if (ld ~= 0)
        aux = onebasisfun__ (u(ii), p-1, U(1:end-1))/ ld;
        N(ii) = N(ii) + ln * aux;
        Nder(ii) = Nder(ii) + p * aux;
      end

      dn = U(end) - u(ii);
      dd = U(end) - U(2);
      if (dd ~= 0)
        aux = onebasisfun__ (u(ii), p-1, U(2:end))/ dd;
        N(ii) = N(ii) + dn * aux;
        Nder(ii) = Nder(ii) - p * aux;
      end
    end
  end
  
end
