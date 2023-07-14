%cut out the white part of the jet colormap to show up NaN data better
JET=jet(128);
jetA=JET(1:43,:);
jetB=JET(64:end,:);
jetNEW = [jetA; jetB];
