function plot_iri(IRI,YC)
%PLOT_IRI IRI plot
%
%   function plot_iri(IRI,YC)
%
%   IRI ... output of iri function
%   YC ... input profile of the iri function

    enlarge_figure(1,2); 
    h = plotyy((IRI(:,1)+IRI(:,2))/2,IRI(:,3),YC(:,1),YC(:,2));
    legend({'IRI','Road profile'});
    ylabel(h(1),'IRI');
    ylabel(h(2),'Elevation [m]');
    xlabel({'Stationing [m]',' '});
    title(sprintf('IRI: %0.2f \\pm %0.2f',mean(IRI(:,3)),std(IRI(:,3))));
end