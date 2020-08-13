function [S,h,p] = synctest( X,varargin )
% [S,h,p] = synctest( X,(alpha,plotFlag) )
%
% Test for synchronization of neural spikes in the matrix X.
% Synchronization is measured as:
%       
%       var( E[X] ) / E[var( X )]
%
% This equation produces a number between 0 and 1, with 1 being completely
% synchronized. This is due to Jensen's inequality:
%
%       f(E[x]) <= E[f(x)]
% 
% The probability of having observed a synchrony S is calculated by
% bootstraping random permutations of the spikes for each neuron (while
% preserving first and second-order statistics).
%
% Inputs:
%   X - neuron x time matrix with each element (i,j)
%       representing a spike (1) or no spike (0) at that time
%
%   (alpha) - 1 - significance level threshold. Default = 0.95
%
%   (plotFlag) - boolean. If true, a figure with the shuffled synchrony 
%                distribution will be plotted. 
%
% Outputs:
%   S - the test statistic (level of synchronization), between [0,1]
%
%   h - boolean with 1 indicating rejection of the null 
%       (neurons are synchronized)
%
%   p - p-val of observing S given a null of no synchronization
%
% Limitations / assumptions:
%   1) currently assumes data in X are binary (0 or 1)
%
% By Jordan Sorokin, 1/25/2018
    
% check if under-sampled
nsp = sum( X~=0,2 );
[n,t] = size( X );
v = var( X,[],2 );
if std( v )/mean( v ) > 1
    warning( 'Undersampled data...interpret results with caution!' );
end

% check optional inputs
if nargin > 1 && ~isempty( varargin{1} )
    alpha = varargin{1};
    assert( alpha < 1, 'alpha must be in the set (0,1)' );
else
    alpha = 0.95; % 1-sided test
end

if nargin > 2 && ~isempty( varargin{2} )
    plotFlag = varargin{2};
else
    plotFlag = false;
end

% pre-allocate 
S_shuffled = zeros( 1,1000 );
X_shuffled = false( n,t );

% generate probability distributions of each neuron's ISI
% so that we may sample from this distribution when shuffling
[i,j] = find( X );
ISI = cell( 1,n );
pISI = cell( 1,n );
for neuron = 1:n
    isi = diff( j(i==neuron) );
    if numel( isi ) < 2
        pISI{neuron} = [];
        continue
    end
    
    % compute kde for this neuron's ISI
    [~,p,ISI{neuron}] = kde( isi,2^10,1,max(isi) );
    pISI{neuron} = cumsum( p / sum( p ) );
end
clear isi p i j v

% compute synchrony of the spike matrix
S = compute_synchrony( X );

% now bootstrap the spike times 1000 times
for bootstrap = 1:1000
    X_shuffled(:) = false;
    
    % randomly permute spikes 
    for neuron = 1:n

        % create our random relative spike times
        if isempty( pISI{neuron} )
            tSpikes = randperm( t,nsp(neuron) );
        else
            dSpikes = ISI{neuron}(sum( bsxfun( @ge,rand( 1,nsp(neuron) ),pISI{neuron} ) ) + 1 );
            tSpikes = round( cumsum( dSpikes ) ); 
            reShuffle = (tSpikes > t);
            tSpikes(reShuffle) = randperm(t,nnz( reShuffle )); % clamp to the max # of points
        end
        
        X_shuffled(neuron,tSpikes) = true;
    end
    
    S_shuffled(bootstrap) = compute_synchrony( X_shuffled );
end

% compute the test statistic / p-value
p = 1 - nnz( S > S_shuffled ) / 1000;
h = p <= (1-alpha);

% plot if desired
if plotFlag 
    figure;
    
    subplot(2,2,1);
    imagesc( ~X );
    title( 'raw data' );
    colormap( 'gray' );
    
    subplot(2,2,2);
    imagesc( ~X_shuffled );
    title( 'shuffled data' );
    colormap( 'gray' );
    
    subplot(2,1,2); hold on
    histogram( S_shuffled,'normalization','probability' );
    ylabel( 'P(S_{shuffled})' );
    xlabel( 'S_{shuffled}' );
    yl = get( gca,'ylim' );
    
    % overlay +/- 1 std 
    sd = std( S_shuffled );
    mu = mean( S_shuffled );
    plot([mu-sd,mu+sd;mu-sd,mu+sd],[yl(1),yl(2)],'k--');
    
    % overlay actual value of synchrony
    plot( [S,S],[yl(1),yl(2)],'r--' );
    set( gca,'tickdir','out','box','off' );
    title( ['p: ',num2str( p )] );
end
%%
    function y = compute_synchrony( x )
        x_hat = mean( x );
        sigma_x = var( x_hat );
        sigma_xi = var( x,[],2 );
        
        y = sigma_x / mean( sigma_xi );
    end
        
end