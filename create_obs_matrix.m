function OM = create_obs_matrix(N,RD,OVP,fit)
  % Build an observation matrix of regularly shifter reports of duration RD with overlap OVP
  % over a time series of length N
  % If fit is set to 0 (default), the observation matrix will have only report vectors
  % with RD ones.
  % Otherwise, it will keep shifting the report till the whole time series is
  % covered


    % Check number of inputs.
    if nargin > 4
        error('obs_matrix:TooManyInputs', ...
            'requires at most 3 optional inputs');
    end

    if RD < 2
        error('For overlapping reports report duration should be at least 2');
    end

    if RD <= OVP
        error('For overlapping reports report duration should be greater than overlap');
    end

    % Fill in unset optional value of fit.
    if nargin == 3
            fit = 0;
    end

    OM = zeros(N,N);   % Obseervation matrix
    OM(1,:)=[ones(RD,1); zeros(N-RD,1)];

    first = RD-OVP+1;
    last = first+RD-1;
    prevlast = last;

    i = 1;

    while last <= N
       Tmp = OM(i+1,:);
       Tmp(first:last) = 1;
       OM(i+1,:) = Tmp;

       prevlast = last;
       first = last-OVP+1;
       last = first+RD-1;
       i=i+1;
    end

   if fit == 1 && RD > 1 && last > N && prevlast < N && first <= N
        Tmp = OM(i+1,:);
        first = prevlast-OVP+1;
        Tmp(first:N) = 1;
        OM(i+1,:) = Tmp;
        i = i+1;
   end

   OM = OM(1:i,:);

end
