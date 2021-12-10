% This example shows the results of two methods for irregular sampling (Sayers cannot be used)
% no overlap
road_profile_irreg = read_profile_2c('profile_2.txt');

seglen = 20; % segment length 20m
overlap = false; 
start = 478.5;
box_filter = true;

disp('Method comparison for irregular sampling, no overlap.');
tic;IRI1 = iri(road_profile_irreg,seglen,start,overlap,box_filter,1);t0 = toc;
tic;IRI2 = iri(road_profile_irreg,seglen,start,overlap,box_filter,2);t2 = toc; 
disp('----------------------------------------');
disp(['Time of evaluation for ode45, Sroubek & Sorel: ']);
disp([num2str([t0 t2]) ' s']);
disp('First ten values for methods ode45, Sroubek & Sorel:');
disp([IRI1(1:10,3)';IRI2(1:10,3)']); 
disp(['Max. IRI difference Sroubek & Sorel - ode45: ' num2str(max(abs(IRI2(:,3)-IRI1(:,3))))]);

figure;plot_iri(IRI2,road_profile);title('Sroubek & Sorel, irregular sampling, no overlap');
