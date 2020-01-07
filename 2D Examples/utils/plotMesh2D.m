function plotMesh2D(PHTelem,controlPts,geometry,writeElem)

p = geometry.p;
q = geometry.q;

if nargin < 4
    writeElem = 0;
end
% We define colors for up to 6 patches
colorArray = {'blue', 'red', 'green', 'cyan', 'magenta', 'black'};

numPts = 11; %number of plot points to use on each edge
uref = linspace(-1,1,numPts);
vref = linspace(-1,1,numPts);

% 1D bernstein polynomials evaluated at the Gauss points on the master element
[B_u] = bernstein_basis(uref,p);
[B_v] = bernstein_basis(vref,q);

R = zeros(numPts, numPts, (p+1)*(q+1));

basisCounter = 0;
for j=1:q+1
    for i=1:p+1
        basisCounter = basisCounter + 1;
        R(:,:,basisCounter) = B_u(:,i)*B_v(:,j)';
    end
end
numPatches = length(PHTelem);

for patchIndex = 1:numPatches
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            colorIndex = rem((patchIndex-1),6)+1;
            % Determine the edges of the element in the physical space
            coord_south = zeros(numPts,2);
            coord_east = zeros(numPts,2);
            coord_north = zeros(numPts,2);
            coord_west = zeros(numPts,2);
            
            nument = min(size(PHTelem{patchIndex}(i).C,1),length(PHTelem{patchIndex}(i).nodes));
            nodes = PHTelem{patchIndex}(i).nodes(1:nument);
            cpts = controlPts{patchIndex}(nodes, 1:2);
            wgts = controlPts{patchIndex}(nodes, 3);
%             plot(cpts(:,1),cpts(:,2),'.r')
            for j=1:numPts
                %compute the points on the element south edge
                N = squeeze(R(j,1,:));
                NS = (PHTelem{patchIndex}(i).C(1:nument,:))*N;
                NS = NS.*wgts;
                w_sum = sum(NS);
                NS = NS/w_sum;
                
                % Compute the points on the element east edge
                N = squeeze(R(numPts,j,:));
                NE = (PHTelem{patchIndex}(i).C(1:nument,:))*N;
                NE = NE.*wgts;
                w_sum = sum(NE);
                NE = NE/w_sum;
                
                % Compute the points on the element north edge
                N = squeeze(R(numPts-j+1,numPts,:));
                NN = (PHTelem{patchIndex}(i).C(1:nument,:))*N;
                NN = NN.*wgts;
                w_sum = sum(NN);
                NN = NN/w_sum;
                
                % Compute the points on the element west edge
                N = squeeze(R(1,numPts-j+1,:));
                NW = (PHTelem{patchIndex}(i).C(1:nument,:))*N;
                NW = NW.*wgts;
                w_sum = sum(NW);
                NW = NW/w_sum;
                
                coord_south(j,:) = NS'*cpts;
                coord_east(j,:) = NE'*cpts;
                coord_north(j,:) = NN'*cpts;
                coord_west(j,:) = NW'*cpts;
            end
            
            % Plot the edges of the element
            coords = [coord_south; coord_east; coord_north; coord_west];
            line(coords(:,1), coords(:,2), 'Color', colorArray{colorIndex})
            
            if writeElem
                %write the element number in the middle
                coord_mid = mean(coords);
                text(coord_mid(1), coord_mid(2), num2str(i), 'Color', colorArray{colorIndex})
            end
        end
    end
end
axis tight
hold on
drawnow
end
