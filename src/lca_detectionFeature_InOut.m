function [feat, fo, EMu, ECov] = lca_detectionFeature_InOut(E,iIn,iOut)

        %Estimate LCA inconsistency statistics from OUTSIDE of region
        EMu = mean(E(iOut,:),1);
        ECov = cov(E(iOut,:));
        
        N = numel(iIn);
        
        %How far do statistics INSIDE region deviate?
        fo = mean(E(iIn,:),1)- EMu;
        feat = N*fo*pinv(ECov)*fo';