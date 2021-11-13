function cluster_info_all = es_getSPMclusterInfo(SPM_filename,thresh)

load(SPM_filename);

SPM.u = thresh;
SPM.k = 0;
SPM.thresDesc = 'none';
SPM.n = [];
SPM.Im = [];
SPM.pm = []; 
SPM.Ex = 0;

clear cluster_info_all
for c=1:length(SPM.xCon)
    SPM.Ic = c;
    
    [~,xSPM] = spm_getSPM(SPM);
    A = spm_clusters(xSPM.XYZ);
    clusters_id = unique(A);
    clear clusters_size clusters_XYZ clusters_XYZmm
    if ~isempty(clusters_id)
        for i=clusters_id
            clusters_size(i) = length(find(A==i));
            clusters_XYZ{i} = xSPM.XYZ(:,find(A==i));
            clusters_XYZmm{i} = xSPM.XYZmm(:,find(A==i));
        end
        [cluster_info_all(c).size,clusters_id_ordered] = sort(clusters_size,'descend');
        for i=1:length(clusters_id)
            cluster_info_all(c).XYZ{i} = clusters_XYZ{clusters_id_ordered(i)};
            cluster_info_all(c).XYZmm{i} = clusters_XYZmm{clusters_id_ordered(i)};
            [~,cluster_ind] = intersect(xSPM.XYZmm',cluster_info_all(c).XYZmm{i}','rows');
            [~,maxi_ind] = max(xSPM.Z(cluster_ind));
            cluster_info_all(c).XYZmm_peak{i} = xSPM.XYZmm(:,cluster_ind(maxi_ind));
        end
    end
    
end