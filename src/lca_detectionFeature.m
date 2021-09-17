function [feat, fo, EMu, ECov] = lca_detectionFeature(E,inds)

        EMu = mean(E,1);
        ECov = cov(E);
        N = numel(inds);
        
        fo = mean(E,1)- EMu;
        feat = N*fo*pinv(ECov)*fo';