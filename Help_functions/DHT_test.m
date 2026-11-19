%% example: generate array of signals in q range, pad zeros to N
 N    = 2^8;               % zero padding
 q    = eps:0.075:10;      % data grid
 a    = linspace(1,20);     
 data = sin(q'*a)./(q'*a); % set of data vectors
 filt = sin(pi*[1:numel(q)]'./(numel(q))).^2;
 data = bsxfun(@times,data, filt); % filter edges with sin^2 window
 [F, Fgrid] = DHT(0.5, data, q, N);
Y = fft2(data,N,100);

 subplot(1,3,1)
 imagesc(a,q,data);
 xlim([a(1) a(end)]);
 ylim([q(1) q(end)]);
 title('data')

 subplot(1,3,2)
 imagesc(a,Fgrid,F);
 xlim([a(1) a(end)]);
 ylim([a(1) a(end)]);
 title('DHT(data)')

  subplot(1,3,3)
 imagesc(abs(Y));
 xlim([a(1) a(end)]);
 ylim([a(1) a(end)]);
 title('DHT(data)')
