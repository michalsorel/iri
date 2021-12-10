% This example shows the results of all three methods for regular sampling
road_profile = read_profile_2c('profile_1.txt');

seglen = 20; % segment length 20m
overlap = 0.5; % overlap 0.5m 
start = 478.5;
box_filter = true;

disp('Method comparison for regular sampling, overlap 0.5m.');
tic;IRI0 = iri(road_profile,seglen,start,overlap,box_filter,0);t0 = toc;
tic;IRI1 = iri(road_profile,seglen,start,overlap,box_filter,1);t1 = toc;
tic;IRI2 = iri(road_profile,seglen,start,overlap,box_filter,2);t2 = toc; 
disp('----------------------------------------');
disp(['Time of evaluation for Sayers, ode45, Sroubek & Sorel: ']);
disp([num2str([t0 t1 t2]) ' s']);
disp('First ten values for methods Sayers, ode45, Sroubek & Sorel:');
disp([IRI0(1:10,3)';IRI1(1:10,3)';IRI2(1:10,3)']); 
disp(['Max. IRI difference Sayers - Sroubek & Sorel: ' num2str(max(abs(IRI0(:,3)-IRI2(:,3))))]);
disp(['Max. IRI difference Sayers - ode45: ' num2str(max(abs(IRI0(:,3)-IRI1(:,3))))]);
disp(['Max. IRI difference Sroubek & Sorel - ode45: ' num2str(max(abs(IRI2(:,3)-IRI1(:,3))))]);

figure;plot_iri(IRI2,road_profile);title('Sroubek & Sorel, regular sampling, overlap 0.5m');
figure;plot_iri(IRI0,road_profile);title('Sayers, regular sampling, overlap 0.5m');
