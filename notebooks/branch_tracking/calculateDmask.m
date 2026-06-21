function [Dmask,out_15] = calculateDmask_x(BW_skeleton, x_value_upper, x_value_lower, x_value_inlet)

   
    % Identifizierung der Verzweigungs- und Endpunkte
    B = bwmorph(BW_skeleton, 'branchpoints');
    E = bwmorph(BW_skeleton, 'endpoints');
    [y, x] = find(E);
    B_loc = find(B);
    
    % Erweiterung von Verzweigungs- und Endpunkten für bessere Visualisierung
    E = imdilate(E, strel('disk', 10));
    B = imdilate(B, strel('disk', 10));
    % figure;
    % imshow(BW_skeleton);
    % figure;
    % imshowpair(E, B);
    % figure;
    % imshow(labeloverlay(double(BW_skeleton), E, 'Colormap', 'autumn', 'Transparency', 0));



    % Initialisierung von Dmask
    Dmask = false(size(BW_skeleton));
    
    % Abschnitt 1: Maskierung basierend auf y-Bereich
    for k = 1:numel(x)
        if x(k) >= x_value_upper && y(k) <= x_value_lower
            D = bwdistgeodesic(BW_skeleton, x(k), y(k));
            distanceToBranchPt = min(D(B_loc));
            Dmask(D < distanceToBranchPt) = true;
        end
    end
    
   
    
    % Finalisiertes Skelett und Ausgabe
    skelD = BW_skeleton - Dmask;
    out_15 = skelD;

    % Visualisierung des Ergebnisses
    % figure;
    % imshowpaiAbmkzr(out_15, BW_skeleton);
end
