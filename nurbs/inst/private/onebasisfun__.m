function Nip = onebasisfun__ (u, p, U)

%  __ONEBASISFUN__: Undocumented internal function
%
%   Adapted from Algorithm A2.4 from 'The NURBS BOOK' pg74.
%
%   Copyright (C) 2009 Carlo de Falco
%   Copyright (C) 2012 Rafael Vazquez
%   This software comes with ABSOLUTELY NO WARRANTY; see the file
%   COPYING for details.  This is free software, and you are welcome
%   to distribute it under the conditions laid out in COPYING.

  Nip = zeros (size (u));
  N = zeros (p+1, 1);

  for ii = 1:numel(u)
    if ((u(ii) == U(1)) && (U(1) == U(end-1)) || ...
        (u(ii) == U(end)) && (U(end) == U(2)))
      Nip(ii) = 1;
      continue
    end
    if (~ any (U <= u(ii))) || (~ any (U > u(ii)))
      continue;
    end
    for jj = 1:p+1 % Initialize zero-th degree functions
      if (u(ii) >= U(jj) && u(ii) < U(jj+1))
        N(jj) = 1;
      else
        N(jj) = 0;
      end
    end
    for k = 1:p
      if (N(1) == 0)
        saved = 0;
      else
        saved = (u(ii) - U(1))*N(1) / (U(k+1)-U(1));
      end

      for jj = 1:p-k+1
        Uleft = U(1+jj);
        Uright = U(1+jj+k);
        if (N(jj+1) == 0)
          N(jj) = saved;
          saved = 0;
        else
          temp = N(jj+1)/(Uright-Uleft);
          N(jj) = saved + (Uright - u(ii))*temp;
          saved = (u(ii) - Uleft)*temp;
        end
      end
    end
    Nip(ii) = N(1);
  end


end
