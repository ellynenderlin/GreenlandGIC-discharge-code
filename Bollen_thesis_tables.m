%CALCULATE SUMMARY STATISTIC BY REGION AND TOTAL
% close all; clear all;

%load term structure
% cd /Users/katebollen/Documents/MS/data/DEMs/centerlines/
cd /Users/ellynenderlin/Research/NASA_GreenlandPeriph-Mapping/
load Greenland_GIC_centerlines.mat;

%loop through each region by BoxID and grab steady state D, 21st cent D,
%SMB adjustment, rate of H change, and r-squared for smb,v,t
%key to regionFlag: 1=WE, 2=SE, 3=CE, 4=NE, 5=NO


%set-up dummy matrices for each region
%WEST
westDearly = []; westDlate = []; westSMBearly = []; westSMBlate = []; 
west_dHdt = []; west_dHdtR = []; west_dHdv = []; west_dHdvR = []; west_dHds = []; west_dHdsR = [];
%SOUTHEAST
southeastDearly = []; southeastDlate = []; southeastSMBearly = []; southeastSMBlate = []; 
southeast_dHdt = []; southeast_dHdtR = []; southeast_dHdv = []; southeast_dHdvR = []; southeast_dHds = []; southeast_dHdsR = [];
%CENTRAL EAST
centraleastDearly = []; centraleastDlate = []; centraleastSMBearly = []; centraleastSMBlate = []; 
centraleast_dHdt = []; centraleast_dHdtR = []; centraleast_dHdv = []; centraleast_dHdvR = []; centraleast_dHds = []; centraleast_dHdsR = [];
%NORTHEAST
northeastDearly = []; northeastDlate = []; northeastSMBearly = []; northeastSMBlate = []; 
northeast_dHdt = []; northeast_dHdtR = []; northeast_dHdv = []; northeast_dHdvR = []; northeast_dHds = []; northeast_dHdsR = [];
%NORTH
northDearly = []; northDlate = []; northSMBearly = []; northSMBlate = []; 
north_dHdt = []; north_dHdtR = []; north_dHdv = []; north_dHdvR = []; north_dHds = []; north_dHdsR = [];
for i=1:length(term)
    disp([num2str(i),' of ',num2str(length(term)),' complete'])
    if term(i).MankoffFlag == 0
        if term(i).regionFlag == 1
            %average annual D
            Dearly = term(i).fluxD(1:14); Dlate = term(i).fluxD(15:end);
            dDearly = term(i).dD_SMB(1:14); dDlate = term(i).dD_SMB(15:end);
            %replace 0s w/ NaNs
            dDearly(Dearly==0) = NaN; dDlate(Dlate==0) = NaN;
            Dearly(Dearly==0) = NaN; Dlate(Dlate==0) = NaN;
            %extract stats
            westDearly = [westDearly; nanmedian(Dearly)];
            westDlate = [westDlate; nanmedian(Dlate)];
            westDearlymad = [westDearlymad; mad(Dearly,1)]; westDlatemad = [westDlatemad; mad(Dlate,1)];
            %average annual dD_SMB
            westSMBearly = [westSMBearly; nanmedian(dDearly)];
            westSMBlate = [westSMBlate; nanmedian(dDlate)];
            westSMBearlymad = [westSMBearlymad; mad(dDearly,1)]; westSMBlatemad = [westSMBlatemad; mad(dDlate,1)];
            clear Dearly dDearly Dlate dDlate;
            
            %thickness change correlation info
            if length(term(i).inlandZmed) > 1 || length(term(i).seawardZmed) > 1
                %dHdt
                if ~isnan(term(i).dZ_tsr2)
                    west_dHdt = [west_dHdt; term(i).dZ_tsfit.p1];
                    if length(term(i).inlandZmed) > 1
                        r = corrcoef(term(i).inlandZyrs', term(i).inlandZmed');
                        west_dHdtR = [west_dHdtR; r(1,2)]; clear r;
                    elseif length(term(i).seawardZmed) > 1
                        r = corrcoef(term(i).seawardZyrs', term(i).seawardZmed');
                        west_dHdtR = [west_dHdtR; r(1,2)]; clear r;
                    end
                end
                
                %dHdV
                if ~isnan(term(i).dZ_vr2)
                    west_dHdv = [west_dHdv; term(i).dZ_vfit.p1];
                    if length(term(i).inlandZmed) > 1
                        for j = 1:length(term(i).inlandZyrs)
                            Vdiff = abs(term(i).inlandZyrs(j) - term(i).gateVdateavg);
                            Vref = find(Vdiff==min(Vdiff));
                            vel(j) = term(i).fluxVavg(Vref);
                            clear Vdiff Vref;
                        end
                        r = corrcoef(vel', term(i).inlandZmed');
                        west_dHdvR = [west_dHdvR; r(1,2)]; clear r;
                    elseif length(term(i).seawardZmed) > 1
                        for j = 1:length(term(i).seawardZyrs)
                            Vdiff = abs(term(i).seawardZyrs(j) - term(i).gateVdateavg);
                            Vref = find(Vdiff==min(Vdiff));
                            vel(j) = term(i).fluxVavg(Vref);
                            clear Vdiff Vref;
                        end
                        r = corrcoef(vel', term(i).seawardZmed');
                        west_dHdvR = [west_dHdvR; r(1,2)]; clear r;
                    end
                    clear vel;
                end
                
                %dHds (s=SMB)
                if ~isnan(term(i).dZ_smbr2)
                    west_dHds = [west_dHds; term(i).dZ_smbfit.p1];
                    if length(term(i).inlandZmed) > 1
                        for j = 1:length(term(i).inlandZyrs)
                            SMBref = find(term(i).SMByrs == floor(term(i).inlandZyrs(j)));
                            smb(j) = term(i).SMB(SMBref);
                            clear SMBref;
                        end
                        r = corrcoef(smb', term(i).inlandZmed');
                        west_dHdsR = [west_dHdsR; r(1,2)]; clear r;
                    elseif length(term(i).seawardZmed) > 1
                        for j = 1:length(term(i).seawardZyrs)
                            SMBref = find(term(i).SMByrs == floor(term(i).seawardZyrs(j)));
                            smb(j) = term(i).SMB(SMBref);
                            clear SMBref;
                        end
                        r = corrcoef(smb', term(i).seawardZmed');
                        west_dHdsR = [west_dHdsR; r(1,2)]; clear r;
                    end
                    clear smb;
                end
                
            else
                west_dHdt = [west_dHdt; NaN]; west_dHdtR = [west_dHdtR; NaN];
                west_dHdv = [west_dHdv; NaN]; west_dHdvR = [west_dHdvR; NaN];
                west_dHds = [west_dHds; NaN]; west_dHdsR = [west_dHdsR; NaN];
            end
            
        elseif term(i).regionFlag == 2
            %average annual D
            Dearly = term(i).fluxD(1:14); Dlate = term(i).fluxD(15:end);
            dDearly = term(i).dD_SMB(1:14); dDlate = term(i).dD_SMB(15:end);
            %replace 0s w/ NaNs
            dDearly(Dearly==0) = NaN; dDlate(Dlate==0) = NaN;
            Dearly(Dearly==0) = NaN; Dlate(Dlate==0) = NaN;
            %extract stats
            southeastDearly = [southeastDearly; nanmedian(Dearly)];
            southeastDlate = [southeastDlate; nanmedian(Dlate)];
            southeastDearlymad = [southeastDearlymad; mad(Dearly,1)]; southeastDlatemad = [southeastDlatemad; mad(Dlate,1)];
            %average annual dD_SMB
            southeastSMBearly = [southeastSMBearly; nanmedian(dDearly)];
            southeastSMBlate = [southeastSMBlate; nanmedian(dDlate)];
            southeastSMBearlymad = [southeastSMBearlymad; mad(dDearly,1)]; southeastSMBlatemad = [southeastSMBlatemad; mad(dDlate,1)];
            clear Dearly dDearly Dlate dDlate;
            
            %thickness change correlation info
            if length(term(i).inlandZmed) > 1 || length(term(i).seawardZmed) > 1
                %dHdt
                if ~isnan(term(i).dZ_tsr2)
                southeast_dHdt = [southeast_dHdt; term(i).dZ_tsfit.p1];
                if length(term(i).inlandZmed) > 1
                    r = corrcoef(term(i).inlandZyrs', term(i).inlandZmed');
                    southeast_dHdtR = [southeast_dHdtR; r(1,2)]; clear r;
                elseif length(term(i).seawardZmed) > 1
                    r = corrcoef(term(i).seawardZyrs', term(i).seawardZmed');
                    southeast_dHdtR = [southeast_dHdtR; r(1,2)]; clear r;
                end
                end
                
                %dHdV
                if ~isnan(term(i).dZ_vr2)
                    southeast_dHdv = [southeast_dHdv; term(i).dZ_vfit.p1];
                    if length(term(i).inlandZmed) > 1
                        for j = 1:length(term(i).inlandZyrs)
                            Vdiff = abs(term(i).inlandZyrs(j) - term(i).gateVdateavg);
                            Vref = find(Vdiff==min(Vdiff));
                            vel(j) = term(i).fluxVavg(Vref);
                            clear Vdiff Vref;
                        end
                        r = corrcoef(vel', term(i).inlandZmed');
                        southeast_dHdvR = [southeast_dHdvR; r(1,2)]; clear r;
                    elseif length(term(i).seawardZmed) > 1
                        for j = 1:length(term(i).seawardZyrs)
                            Vdiff = abs(term(i).seawardZyrs(j) - term(i).gateVdateavg);
                            Vref = find(Vdiff==min(Vdiff));
                            vel(j) = term(i).fluxVavg(Vref);
                            clear Vdiff Vref;
                        end
                        r = corrcoef(vel', term(i).seawardZmed');
                        southeast_dHdvR = [southeast_dHdvR; r(1,2)]; clear r;
                    end
                    clear vel;
                end
                
                %dHds (s=SMB)
                if ~isnan(term(i).dZ_smbr2)
                    southeast_dHds = [southeast_dHds; term(i).dZ_smbfit.p1];
                    if length(term(i).inlandZmed) > 1
                        for j = 1:length(term(i).inlandZyrs)
                            SMBref = find(term(i).SMByrs == floor(term(i).inlandZyrs(j)));
                            smb(j) = term(i).SMB(SMBref);
                            clear SMBref;
                        end
                        r = corrcoef(smb', term(i).inlandZmed');
                        southeast_dHdsR = [southeast_dHdsR; r(1,2)]; clear r;
                    elseif length(term(i).seawardZmed) > 1
                        for j = 1:length(term(i).seawardZyrs)
                            SMBref = find(term(i).SMByrs == floor(term(i).seawardZyrs(j)));
                            smb(j) = term(i).SMB(SMBref);
                            clear SMBref;
                        end
                        r = corrcoef(smb', term(i).seawardZmed');
                        southeast_dHdsR = [southeast_dHdsR; r(1,2)]; clear r;
                    end
                    clear smb;
                end
                
            else
                southeast_dHdt = [southeast_dHdt; NaN]; southeast_dHdtR = [southeast_dHdtR; NaN];
                southeast_dHdv = [southeast_dHdv; NaN]; southeast_dHdvR = [southeast_dHdvR; NaN];
                southeast_dHds = [southeast_dHds; NaN]; southeast_dHdsR = [southeast_dHdsR; NaN];
            end
            
        elseif term(i).regionFlag == 3
            %average annual D
            Dearly = term(i).fluxD(1:14); Dlate = term(i).fluxD(15:end);
            dDearly = term(i).dD_SMB(1:14); dDlate = term(i).dD_SMB(15:end);
            %replace 0s w/ NaNs
            dDearly(Dearly==0) = NaN; dDlate(Dlate==0) = NaN;
            Dearly(Dearly==0) = NaN; Dlate(Dlate==0) = NaN;
            %extract stats
            centraleastDearly = [centraleastDearly; nanmedian(Dearly)];
            centraleastDlate = [centraleastDlate; nanmedian(Dlate)];
            centraleastDearlymad = [centraleastDearlymad; mad(Dearly,1)]; centraleastDlatemad = [centraleastDlatemad; mad(Dlate,1)];
            %average annual dD_SMB
            centraleastSMBearly = [centraleastSMBearly; nanmedian(dDearly)];
            centraleastSMBlate = [centraleastSMBlate; nanmedian(dDlate)];
            centraleastSMBearlymad = [centraleastSMBearlymad; mad(dDearly,1)]; centraleastSMBlatemad = [centraleastSMBlatemad; mad(dDlate,1)];
            clear Dearly dDearly Dlate dDlate;
            
            %thickness change correlation info
            if length(term(i).inlandZmed) > 1 || length(term(i).seawardZmed) > 1
                %dHdt
                if ~isnan(term(i).dZ_tsr2)
                centraleast_dHdt = [centraleast_dHdt; term(i).dZ_tsfit.p1];
                if length(term(i).inlandZmed) > 1
                    r = corrcoef(term(i).inlandZyrs', term(i).inlandZmed');
                    centraleast_dHdtR = [centraleast_dHdtR; r(1,2)]; clear r;
                elseif length(term(i).seawardZmed) > 1
                    r = corrcoef(term(i).seawardZyrs', term(i).seawardZmed');
                    centraleast_dHdtR = [centraleast_dHdtR; r(1,2)]; clear r;
                end
                end
                
                %dHdV
                if ~isnan(term(i).dZ_vr2)
                    centraleast_dHdv = [centraleast_dHdv; term(i).dZ_vfit.p1];
                    if length(term(i).inlandZmed) > 1
                        for j = 1:length(term(i).inlandZyrs)
                            Vdiff = abs(term(i).inlandZyrs(j) - term(i).gateVdateavg);
                            Vref = find(Vdiff==min(Vdiff));
                            vel(j) = term(i).fluxVavg(Vref);
                            clear Vdiff Vref;
                        end
                        r = corrcoef(vel', term(i).inlandZmed');
                        centraleast_dHdvR = [centraleast_dHdvR; r(1,2)]; clear r;
                    elseif length(term(i).seawardZmed) > 1
                        for j = 1:length(term(i).seawardZyrs)
                            Vdiff = abs(term(i).seawardZyrs(j) - term(i).gateVdateavg);
                            Vref = find(Vdiff==min(Vdiff));
                            vel(j) = term(i).fluxVavg(Vref);
                            clear Vdiff Vref;
                        end
                        r = corrcoef(vel', term(i).seawardZmed');
                        centraleast_dHdvR = [centraleast_dHdvR; r(1,2)]; clear r;
                    end
                    clear vel;
                end
                
                %dHds (s=SMB)
                if ~isnan(term(i).dZ_smbr2)
                    centraleast_dHds = [centraleast_dHds; term(i).dZ_smbfit.p1];
                    if length(term(i).inlandZmed) > 1
                        for j = 1:length(term(i).inlandZyrs)
                            SMBref = find(term(i).SMByrs == floor(term(i).inlandZyrs(j)));
                            smb(j) = term(i).SMB(SMBref);
                            clear SMBref;
                        end
                        r = corrcoef(smb', term(i).inlandZmed');
                        centraleast_dHdsR = [centraleast_dHdsR; r(1,2)]; clear r;
                    elseif length(term(i).seawardZmed) > 1
                        for j = 1:length(term(i).seawardZyrs)
                            SMBref = find(term(i).SMByrs == floor(term(i).seawardZyrs(j)));
                            smb(j) = term(i).SMB(SMBref);
                            clear SMBref;
                        end
                        r = corrcoef(smb', term(i).seawardZmed');
                        centraleast_dHdsR = [centraleast_dHdsR; r(1,2)]; clear r;
                    end
                    clear smb;
                end
                
            else
                centraleast_dHdt = [centraleast_dHdt; NaN]; centraleast_dHdtR = [centraleast_dHdtR; NaN];
                centraleast_dHdv = [centraleast_dHdv; NaN]; centraleast_dHdvR = [centraleast_dHdvR; NaN];
                centraleast_dHds = [centraleast_dHds; NaN]; centraleast_dHdsR = [centraleast_dHdsR; NaN];
            end
            
        elseif term(i).regionFlag == 4
            %average annual D
            Dearly = term(i).fluxD(1:14); Dlate = term(i).fluxD(15:end);
            dDearly = term(i).dD_SMB(1:14); dDlate = term(i).dD_SMB(15:end);
            %replace 0s w/ NaNs
            dDearly(Dearly==0) = NaN; dDlate(Dlate==0) = NaN;
            Dearly(Dearly==0) = NaN; Dlate(Dlate==0) = NaN;
            %extract stats
            northeastDearly = [northeastDearly; nanmedian(Dearly)];
            northeastDlate = [northeastDlate; nanmedian(Dlate)];
            northeastDearlymad = [northeastDearlymad; mad(Dearly,1)]; northeastDlatemad = [northeastDlatemad; mad(Dlate,1)];
            %average annual dD_SMB
            northeastSMBearly = [northeastSMBearly; nanmedian(dDearly)];
            northeastSMBlate = [northeastSMBlate; nanmedian(dDlate)];
            northeastSMBearlymad = [northeastSMBearlymad; mad(dDearly,1)]; northeastSMBlatemad = [northeastSMBlatemad; mad(dDlate,1)];
            clear Dearly dDearly Dlate dDlate;
            
            %thickness change correlation info
            if length(term(i).inlandZmed) > 1 || length(term(i).seawardZmed) > 1
                %dHdt
                if ~isnan(term(i).dZ_tsr2)
                northeast_dHdt = [northeast_dHdt; term(i).dZ_tsfit.p1];
                if length(term(i).inlandZmed) > 1
                    r = corrcoef(term(i).inlandZyrs', term(i).inlandZmed');
                    northeast_dHdtR = [northeast_dHdtR; r(1,2)]; clear r;
                elseif length(term(i).seawardZmed) > 1
                    r = corrcoef(term(i).seawardZyrs', term(i).seawardZmed');
                    northeast_dHdtR = [northeast_dHdtR; r(1,2)]; clear r;
                end
                end
                
                %dHdV
                if ~isnan(term(i).dZ_vr2)
                    northeast_dHdv = [northeast_dHdv; term(i).dZ_vfit.p1];
                    if length(term(i).inlandZmed) > 1
                        for j = 1:length(term(i).inlandZyrs)
                            Vdiff = abs(term(i).inlandZyrs(j) - term(i).gateVdateavg);
                            Vref = find(Vdiff==min(Vdiff));
                            vel(j) = term(i).fluxVavg(Vref);
                            clear Vdiff Vref;
                        end
                        r = corrcoef(vel', term(i).inlandZmed');
                        northeast_dHdvR = [northeast_dHdvR; r(1,2)]; clear r;
                    elseif length(term(i).seawardZmed) > 1
                        for j = 1:length(term(i).seawardZyrs)
                            Vdiff = abs(term(i).seawardZyrs(j) - term(i).gateVdateavg);
                            Vref = find(Vdiff==min(Vdiff));
                            vel(j) = term(i).fluxVavg(Vref);
                            clear Vdiff Vref;
                        end
                        r = corrcoef(vel', term(i).seawardZmed');
                        northeast_dHdvR = [northeast_dHdvR; r(1,2)]; clear r;
                    end
                    clear vel;
                end
                
                %dHds (s=SMB)
                if ~isnan(term(i).dZ_smbr2)
                    northeast_dHds = [northeast_dHds; term(i).dZ_smbfit.p1];
                    if length(term(i).inlandZmed) > 1
                        for j = 1:length(term(i).inlandZyrs)
                            SMBref = find(term(i).SMByrs == floor(term(i).inlandZyrs(j)));
                            smb(j) = term(i).SMB(SMBref);
                            clear SMBref;
                        end
                        r = corrcoef(smb', term(i).inlandZmed');
                        northeast_dHdsR = [northeast_dHdsR; r(1,2)]; clear r;
                    elseif length(term(i).seawardZmed) > 1
                        for j = 1:length(term(i).seawardZyrs)
                            SMBref = find(term(i).SMByrs == floor(term(i).seawardZyrs(j)));
                            smb(j) = term(i).SMB(SMBref);
                            clear SMBref;
                        end
                        r = corrcoef(smb', term(i).seawardZmed');
                        northeast_dHdsR = [northeast_dHdsR; r(1,2)]; clear r;
                    end
                    clear smb;
                end
                
            else
                northeast_dHdt = [northeast_dHdt; NaN]; northeast_dHdtR = [northeast_dHdtR; NaN];
                northeast_dHdv = [northeast_dHdv; NaN]; northeast_dHdvR = [northeast_dHdvR; NaN];
                northeast_dHds = [northeast_dHds; NaN]; northeast_dHdsR = [northeast_dHdsR; NaN];
            end
            
        else
            %average annual D
            Dearly = term(i).fluxD(1:14); Dlate = term(i).fluxD(15:end);
            dDearly = term(i).dD_SMB(1:14); dDlate = term(i).dD_SMB(15:end);
            %replace 0s w/ NaNs
            dDearly(Dearly==0) = NaN; dDlate(Dlate==0) = NaN;
            Dearly(Dearly==0) = NaN; Dlate(Dlate==0) = NaN;
            %extract stats
            northDearly = [northDearly; nanmedian(Dearly)];
            northDlate = [northDlate; nanmedian(Dlate)];
            northDearlymad = [northDearlymad; mad(Dearly,1)]; northDlatemad = [northDlatemad; mad(Dlate,1)];
            %average annual dD_SMB
            northSMBearly = [northSMBearly; nanmedian(dDearly)];
            northSMBlate = [northSMBlate; nanmedian(dDlate)];
            northSMBearlymad = [northSMBearlymad; mad(dDearly,1)]; northSMBlatemad = [northSMBlatemad; mad(dDlate,1)];
            clear Dearly dDearly Dlate dDlate;
            
            %thickness change correlation info
            if length(term(i).inlandZmed) > 1 || length(term(i).seawardZmed) > 1
                %dHdt
                if ~isnan(term(i).dZ_tsr2)
                north_dHdt = [north_dHdt; term(i).dZ_tsfit.p1];
                if length(term(i).inlandZmed) > 1
                    r = corrcoef(term(i).inlandZyrs', term(i).inlandZmed');
                    north_dHdtR = [north_dHdtR; r(1,2)]; clear r;
                elseif length(term(i).seawardZmed) > 1
                    r = corrcoef(term(i).seawardZyrs', term(i).seawardZmed');
                    north_dHdtR = [north_dHdtR; r(1,2)]; clear r;
                end
                end
                
                %dHdV
                if ~isnan(term(i).dZ_vr2)
                    north_dHdv = [north_dHdv; term(i).dZ_vfit.p1];
                    if length(term(i).inlandZmed) > 1
                        for j = 1:length(term(i).inlandZyrs)
                            Vdiff = abs(term(i).inlandZyrs(j) - term(i).gateVdateavg);
                            Vref = find(Vdiff==min(Vdiff));
                            vel(j) = term(i).fluxVavg(Vref);
                            clear Vdiff Vref;
                        end
                        r = corrcoef(vel', term(i).inlandZmed');
                        north_dHdvR = [north_dHdvR; r(1,2)]; clear r;
                    elseif length(term(i).seawardZmed) > 1
                        for j = 1:length(term(i).seawardZyrs)
                            Vdiff = abs(term(i).seawardZyrs(j) - term(i).gateVdateavg);
                            Vref = find(Vdiff==min(Vdiff));
                            vel(j) = term(i).fluxVavg(Vref);
                            clear Vdiff Vref;
                        end
                        r = corrcoef(vel', term(i).seawardZmed');
                        north_dHdvR = [north_dHdvR; r(1,2)]; clear r;
                    end
                    clear vel;
                end
                
                %dHds (s=SMB)
                if ~isnan(term(i).dZ_smbr2)
                    north_dHds = [north_dHds; term(i).dZ_smbfit.p1];
                    if length(term(i).inlandZmed) > 1
                        for j = 1:length(term(i).inlandZyrs)
                            SMBref = find(term(i).SMByrs == floor(term(i).inlandZyrs(j)));
                            smb(j) = term(i).SMB(SMBref);
                            clear SMBref;
                        end
                        r = corrcoef(smb', term(i).inlandZmed');
                        north_dHdsR = [north_dHdsR; r(1,2)]; clear r;
                    elseif length(term(i).seawardZmed) > 1
                        for j = 1:length(term(i).seawardZyrs)
                            SMBref = find(term(i).SMByrs == floor(term(i).seawardZyrs(j)));
                            smb(j) = term(i).SMB(SMBref);
                            clear SMBref;
                        end
                        r = corrcoef(smb', term(i).seawardZmed');
                        north_dHdsR = [north_dHdsR; r(1,2)]; clear r;
                    end
                    clear smb;
                end
                
            else
                north_dHdt = [north_dHdt; NaN]; north_dHdtR = [north_dHdtR; NaN];
                north_dHdv = [north_dHdv; NaN]; north_dHdvR = [north_dHdvR; NaN];
                north_dHds = [north_dHds; NaN]; north_dHdsR = [north_dHdsR; NaN];
            end
            
        end
        
    end
end

%statistics
%WEST
%sum time-averaged regional D 1985-1998 (pre) and 1999-2018 (post)
west_preD = nansum(westDearly).*917./(1000*10^9); %Gt/yr
west_postD = nansum(westDlate).*917./(1000*10^9); %Gt/yr
disp(['W discharge median = ',num2str(round(west_preD,2)),' (steady-state) & ',num2str(round(west_postD,2)),' (21st century) Gt/yr']);
disp(['W discharge MAD = ',num2str(sqrt(nansum(westDearlymad.^2)).*917./(1000*10^9)),' (steady-state) & ',num2str(sqrt(nansum(westDlatemad.^2)).*917./(1000*10^9)),' (21st century) Gt/yr']);
%sum time-averaged SMB correction
west_preSMBdD = nansum(westSMBearly).*917./(1000*10^9); %Gt/yr
west_postSMBdD = nansum(westSMBlate).*917./(1000*10^9); %Gt/yr
disp(['W discharge SMB correction median = ',num2str(round(west_preSMBdD,2)),' (steady-state) & ',num2str(round(west_postSMBdD,2)),' (21st century) Gt/yr']);
disp(['W discharge SMB correction MAD = ',num2str(sqrt(nansum(westSMBearlymad.^2)).*917./(1000*10^9)),' (steady-state) & ',num2str(sqrt(nansum(westSMBlatemad.^2)).*917./(1000*10^9)),' (21st century) Gt/yr']);
disp(['W % glaciers thinning = ', num2str((length((west_dHdt(west_dHdt < 0)))/length(west_dHdt))*100)]);
%median dH/dt and correlation coefficient
west_dHdtmed = nanmedian(west_dHdt); west_dHdtMAD = mad(west_dHdt(~isnan(west_dHdt)),1);
west_dHdtRmed = nanmedian(west_dHdtR); west_dHdtRMAD = mad(west_dHdtR(~isnan(west_dHdt)),1);
disp(['W dH/dt = ',num2str(round(west_dHdtmed,2)),'+/- ',num2str(west_dHdtMAD),' m/yr (R =',num2str(round(west_dHdtRmed,2)),')']);
%median dH/dspeed and correlation coefficient
west_dHdvmed = nanmedian(west_dHdv); west_dHdvMAD = mad(west_dHdv(~isnan(west_dHdv)),1);
west_dHdvRmed = nanmedian(west_dHdvR); west_dHdvRMAD = mad(west_dHdvR(~isnan(west_dHdv)),1);
disp(['W dH/dspeed = ',num2str(round(west_dHdvmed,2)),' m/yr (R =',num2str(round(west_dHdvRmed,2)),')']);
%median dH/dSMB and correlation coefficient
west_dHdsmed = nanmedian(west_dHds); west_dHdsMAD = mad(west_dHds(~isnan(west_dHds)),1);
west_dHdsRmed = nanmedian(west_dHdsR); west_dHdsRMAD = mad(west_dHdsR(~isnan(west_dHds)),1);
disp(['W dH/dSMB = ',num2str(round(west_dHdsmed,2)),' m/yr (R =',num2str(round(west_dHdsRmed,2)),')']);


%SOUTHEAST
southeast_preD = nansum(southeastDearly).*917./(1000*10^9); %Gt/yr
southeast_postD = nansum(southeastDlate).*917./(1000*10^9); %Gt/yr
disp(['SE discharge median = ',num2str(round(southeast_preD,2)),' (steady-state) & ',num2str(round(southeast_postD,2)),' (21st century) Gt/yr']);
disp(['SE discharge MAD = ',num2str(sqrt(nansum(southeastDearlymad.^2)).*917./(1000*10^9)),' (steady-state) & ',num2str(sqrt(nansum(southeastDlatemad.^2)).*917./(1000*10^9)),' (21st century) Gt/yr']);
%sum time-averaged SMB correction
southeast_preSMBdD = nansum(southeastSMBearly).*917./(1000*10^9); %Gt/yr
southeast_postSMBdD = nansum(southeastSMBlate).*917./(1000*10^9); %Gt/yr
disp(['SE discharge SMB correction median = ',num2str(round(southeast_preSMBdD,2)),' (steady-state) & ',num2str(round(southeast_postSMBdD,2)),' (21st century) Gt/yr']);
disp(['SE discharge SMB correction MAD = ',num2str(sqrt(nansum(southeastSMBearlymad.^2)).*917./(1000*10^9)),' (steady-state) & ',num2str(sqrt(nansum(southeastSMBlatemad.^2)).*917./(1000*10^9)),' (21st century) Gt/yr']);
disp(['SE % glaciers thinning = ', num2str((length((southeast_dHdt(southeast_dHdt < 0)))/length(southeast_dHdt))*100)]);
%median dH/dt and r^2
southeast_dHdtmed = nanmedian(southeast_dHdt); southeast_dHdtMAD = mad(southeast_dHdt(~isnan(southeast_dHdt)),1);
southeast_dHdtRmed = nanmedian(southeast_dHdtR); southeast_dHdtRMAD = mad(southeast_dHdtR(~isnan(southeast_dHdt)),1);
disp(['SE dH/dt = ',num2str(round(southeast_dHdtmed,2)),'+/- ',num2str(southeast_dHdtMAD),' m/yr (R =',num2str(round(southeast_dHdtRmed,2)),')']);
%median dH/dspeed and correlation coefficient
southeast_dHdvmed = nanmedian(southeast_dHdv); southeast_dHdvMAD = mad(southeast_dHdv(~isnan(southeast_dHdv)),1);
southeast_dHdvRmed = nanmedian(southeast_dHdvR); southeast_dHdvRMAD = mad(southeast_dHdvR(~isnan(southeast_dHdv)),1);
disp(['SE dH/dspeed = ',num2str(round(southeast_dHdvmed,2)),' m/yr (R =',num2str(round(southeast_dHdvRmed,2)),')']);
%median dH/dSMB and correlation coefficient
southeast_dHdsmed = nanmedian(southeast_dHds); southeast_dHdsMAD = mad(southeast_dHds(~isnan(southeast_dHds)),1);
southeast_dHdsRmed = nanmedian(southeast_dHdsR); southeast_dHdsRMAD = mad(southeast_dHdsR(~isnan(southeast_dHds)),1);
disp(['SE dH/dSMB = ',num2str(round(southeast_dHdsmed,2)),' m/yr (R =',num2str(round(southeast_dHdsRmed,2)),')']);

%CENTRAL EAST
centraleast_preD = nansum(centraleastDearly).*917./(1000*10^9); %Gt/yr
centraleast_postD = nansum(centraleastDlate).*917./(1000*10^9); %Gt/yr
disp(['CE discharge median = ',num2str(round(centraleast_preD,2)),' (steady-state) & ',num2str(round(centraleast_postD,2)),' (21st century) Gt/yr']);
disp(['CE discharge MAD = ',num2str(sqrt(nansum(centraleastDearlymad.^2)).*917./(1000*10^9)),' (steady-state) & ',num2str(sqrt(nansum(centraleastDlatemad.^2)).*917./(1000*10^9)),' (21st century) Gt/yr']);
%sum time-averaged SMB correction
centraleast_preSMBdD = nansum(centraleastSMBearly).*917./(1000*10^9); %Gt/yr
centraleast_postSMBdD = nansum(centraleastSMBlate).*917./(1000*10^9); %Gt/yr
disp(['CE discharge SMB correction median = ',num2str(round(centraleast_preSMBdD,2)),' (steady-state) & ',num2str(round(centraleast_postSMBdD,2)),' (21st century) Gt/yr']);
disp(['CE discharge SMB correction MAD = ',num2str(sqrt(nansum(centraleastSMBearlymad.^2)).*917./(1000*10^9)),' (steady-state) & ',num2str(sqrt(nansum(centraleastSMBlatemad.^2)).*917./(1000*10^9)),' (21st century) Gt/yr']);
disp(['CE % glaciers thinning = ', num2str((length((centraleast_dHdt(centraleast_dHdt < 0)))/length(centraleast_dHdt))*100)]);
%median dH/dt and r^2
centraleast_dHdtmed = nanmedian(centraleast_dHdt); centraleast_dHdtMAD = mad(centraleast_dHdt(~isnan(centraleast_dHdt)),1);
centraleast_dHdtRmed = nanmedian(centraleast_dHdtR); centraleast_dHdtRMAD = mad(centraleast_dHdtR(~isnan(centraleast_dHdt)),1);
disp(['CE dH/dt = ',num2str(round(centraleast_dHdtmed,2)),'+/- ',num2str(centraleast_dHdtMAD),' m/yr (R =',num2str(round(centraleast_dHdtRmed,2)),')']);
%median dH/dspeed and correlation coefficient
centraleast_dHdvmed = nanmedian(centraleast_dHdv); centraleast_dHdvMAD = mad(centraleast_dHdv(~isnan(centraleast_dHdv)),1);
centraleast_dHdvRmed = nanmedian(centraleast_dHdvR); centraleast_dHdvRMAD = mad(centraleast_dHdvR(~isnan(centraleast_dHdv)),1);
disp(['CE dH/dspeed = ',num2str(round(centraleast_dHdvmed,2)),' m/yr (R =',num2str(round(centraleast_dHdvRmed,2)),')']);
%median dH/dSMB and correlation coefficient
centraleast_dHdsmed = nanmedian(centraleast_dHds); centraleast_dHdsMAD = mad(centraleast_dHds(~isnan(centraleast_dHds)),1);
centraleast_dHdsRmed = nanmedian(centraleast_dHdsR); centraleast_dHdsRMAD = mad(centraleast_dHdsR(~isnan(centraleast_dHds)),1);
disp(['CE dH/dSMB = ',num2str(round(centraleast_dHdsmed,2)),' m/yr (R =',num2str(round(centraleast_dHdsRmed,2)),')']);

%NORTHEAST
northeast_preD = nansum(northeastDearly).*917./(1000*10^9); %Gt/yr
northeast_postD = nansum(northeastDlate).*917./(1000*10^9); %Gt/yr
disp(['NE discharge median = ',num2str(round(northeast_preD,2)),' (steady-state) & ',num2str(round(northeast_postD,2)),' (21st century) Gt/yr']);
disp(['NE discharge MAD = ',num2str(sqrt(nansum(northeastDearlymad.^2)).*917./(1000*10^9)),' (steady-state) & ',num2str(sqrt(nansum(northeastDlatemad.^2)).*917./(1000*10^9)),' (21st century) Gt/yr']);
%sum time-averaged SMB correction
northeast_preSMBdD = nansum(northeastSMBearly).*917./(1000*10^9); %Gt/yr
northeast_postSMBdD = nansum(northeastSMBlate).*917./(1000*10^9); %Gt/yr
disp(['NE discharge SMB correction median = ',num2str(round(northeast_preSMBdD,2)),' (steady-state) & ',num2str(round(northeast_postSMBdD,2)),' (21st century) Gt/yr']);
disp(['NE discharge SMB correction MAD = ',num2str(sqrt(nansum(northeastSMBearlymad.^2)).*917./(1000*10^9)),' (steady-state) & ',num2str(sqrt(nansum(northeastSMBlatemad.^2)).*917./(1000*10^9)),' (21st century) Gt/yr']);
disp(['NE % glaciers thinning = ', num2str((length((northeast_dHdt(northeast_dHdt < 0)))/length(northeast_dHdt))*100)]);
%median dH/dt and r^2
northeast_dHdtmed = nanmedian(northeast_dHdt); northeast_dHdtMAD = mad(northeast_dHdt(~isnan(northeast_dHdt)),1);
northeast_dHdtRmed = nanmedian(northeast_dHdtR); northeast_dHdtRMAD = mad(northeast_dHdtR(~isnan(northeast_dHdt)),1);
disp(['NE dH/dt = ',num2str(round(northeast_dHdtmed,2)),'+/- ',num2str(northeast_dHdtMAD),' m/yr (R =',num2str(round(northeast_dHdtRmed,2)),')']);
%median dH/dspeed and correlation coefficient
northeast_dHdvmed = nanmedian(northeast_dHdv); northeast_dHdvMAD = mad(northeast_dHdv(~isnan(northeast_dHdv)),1);
northeast_dHdvRmed = nanmedian(northeast_dHdvR); northeast_dHdvRMAD = mad(northeast_dHdvR(~isnan(northeast_dHdv)),1);
disp(['NE dH/dspeed = ',num2str(round(northeast_dHdvmed,2)),' m/yr (R =',num2str(round(northeast_dHdvRmed,2)),')']);
%median dH/dSMB and correlation coefficient
northeast_dHdsmed = nanmedian(northeast_dHds); northeast_dHdsMAD = mad(northeast_dHds(~isnan(northeast_dHds)),1);
northeast_dHdsRmed = nanmedian(northeast_dHdsR); northeast_dHdsRMAD = mad(northeast_dHdsR(~isnan(northeast_dHds)),1);
disp(['NE dH/dSMB = ',num2str(round(northeast_dHdsmed,2)),' m/yr (R =',num2str(round(northeast_dHdsRmed,2)),')']);

%NORTH
north_preD = nansum(northDearly).*917./(1000*10^9); %Gt/yr
north_postD = nansum(northDlate).*917./(1000*10^9); %Gt/yr
disp(['N discharge median = ',num2str(round(north_preD,2)),' (steady-state) & ',num2str(round(north_postD,2)),' (21st century) Gt/yr']);
disp(['N discharge MAD = ',num2str(sqrt(nansum(northDearlymad.^2)).*917./(1000*10^9)),' (steady-state) & ',num2str(sqrt(nansum(northDlatemad.^2)).*917./(1000*10^9)),' (21st century) Gt/yr']);
%sum time-averaged SMB correction
north_preSMBdD = nansum(northSMBearly).*917./(1000*10^9); %Gt/yr
north_postSMBdD = nansum(northSMBlate).*917./(1000*10^9); %Gt/yr
disp(['N discharge SMB correction median = ',num2str(round(north_preSMBdD,2)),' (steady-state) & ',num2str(round(north_postSMBdD,2)),' (21st century) Gt/yr']);
disp(['N discharge SMB correction MAD = ',num2str(sqrt(nansum(northSMBearlymad.^2)).*917./(1000*10^9)),' (steady-state) & ',num2str(sqrt(nansum(northSMBlatemad.^2)).*917./(1000*10^9)),' (21st century) Gt/yr']);
disp(['N % glaciers thinning = ', num2str((length((north_dHdt(north_dHdt < 0)))/length(north_dHdt))*100)]);
%median dH/dt and r^2
north_dHdtmed = nanmedian(north_dHdt); north_dHdtMAD = mad(north_dHdt(~isnan(north_dHdt)),1);
north_dHdtRmed = nanmedian(north_dHdtR); north_dHdtRMAD = mad(north_dHdtR(~isnan(north_dHdt)),1);
disp(['N dH/dt = ',num2str(round(north_dHdtmed,2)),'+/- ',num2str(north_dHdtMAD),' m/yr (R =',num2str(round(north_dHdtRmed,2)),')']);
%median dH/dspeed and correlation coefficient
north_dHdvmed = nanmedian(north_dHdv); north_dHdvMAD = mad(north_dHdv(~isnan(north_dHdv)),1);
north_dHdvRmed = nanmedian(north_dHdvR); north_dHdvRMAD = mad(north_dHdvR(~isnan(north_dHdv)),1);
disp(['N dH/dspeed = ',num2str(round(north_dHdvmed,2)),' m/yr (R =',num2str(round(north_dHdvRmed,2)),')']);
%median dH/dSMB and correlation coefficient
north_dHdsmed = nanmedian(north_dHds); north_dHdsMAD = mad(north_dHds(~isnan(north_dHds)),1);
north_dHdsRmed = nanmedian(north_dHdsR); north_dHdsRMAD = mad(north_dHdsR(~isnan(north_dHds)),1);
disp(['N dH/dSMB = ',num2str(round(north_dHdsmed,2)),' m/yr (R =',num2str(round(north_dHdsRmed,2)),')']);


% ALL REGIONS
disp('Obtain total discharge from GreenlandGIC_discharge_timeseries.m');
GIC_preSMBdD = nansum([westSMBearly; southeastSMBearly; centraleastSMBearly; northeastSMBearly; northSMBearly]).*917./(1000*10^9); %Gt/yr
GIC_postSMBdD = nansum([westSMBlate; southeastSMBlate; centraleastSMBlate; northeastSMBlate; northSMBlate]).*917./(1000*10^9); %Gt/yr
disp(['W discharge SMB correction median = ',num2str(round(GIC_preSMBdD,2)),' (steady-state) & ',num2str(round(GIC_postSMBdD,2)),' (21st century) Gt/yr']);
disp(['W discharge SMB correction MAD = ',...
    num2str(sqrt(nansum([westSMBearlymad; southeastSMBearlymad; centraleastSMBearlymad; northeastSMBearlymad; northSMBearlymad].^2)).*917./(1000*10^9)),' (steady-state) & ',...
    num2str(sqrt(nansum([westSMBlatemad; southeastSMBlatemad; centraleastSMBlatemad; northeastSMBlatemad; northSMBlatemad].^2)).*917./(1000*10^9)),' (21st century) Gt/yr']);
%median dH/dt and r^2
GIC_dHdtmed = nanmedian([west_dHdt; southeast_dHdt; centraleast_dHdt; northeast_dHdt; north_dHdt]); 
GIC_dHdtMAD = mad([west_dHdt; southeast_dHdt; centraleast_dHdt; northeast_dHdt; north_dHdt],1);
GIC_dHdtRmed = nanmedian([west_dHdtR; southeast_dHdtR; centraleast_dHdtR; northeast_dHdtR; north_dHdtR]); 
GIC_dHdtRMAD = mad([west_dHdtR; southeast_dHdtR; centraleast_dHdtR; northeast_dHdtR; north_dHdtR],1);
disp(['GIC dH/dt = ',num2str(round(GIC_dHdtmed,2)),'+/- ',num2str(GIC_dHdtMAD),' m/yr (R =',num2str(round(GIC_dHdtRmed,2)),')']);
%median dH/dspeed and correlation coefficient
GIC_dHdvmed = nanmedian([west_dHdv; southeast_dHdv; centraleast_dHdv; northeast_dHdv; north_dHdv]); 
GIC_dHdvMAD = mad([west_dHdv; southeast_dHdv; centraleast_dHdv; northeast_dHdv; north_dHdv],1);
GIC_dHdvRmed = nanmedian([west_dHdvR; southeast_dHdvR; centraleast_dHdvR; northeast_dHdvR; north_dHdvR]); 
GIC_dHdvRMAD = mad([west_dHdvR; southeast_dHdvR; centraleast_dHdvR; northeast_dHdvR; north_dHdvR],1);
disp(['GIC dH/dspeed = ',num2str(round(GIC_dHdvmed,2)),' m/yr (R =',num2str(round(GIC_dHdvRmed,2)),')']);
%median dH/dSMB and correlation coefficient
GIC_dHdsmed = nanmedian([west_dHds; southeast_dHds; centraleast_dHds; northeast_dHds; north_dHds]); 
GIC_dHdsMAD = mad([west_dHds; southeast_dHds; centraleast_dHds; northeast_dHds; north_dHds],1);
GIC_dHdsRmed = nanmedian([west_dHdsR; southeast_dHdsR; centraleast_dHdsR; northeast_dHdsR; north_dHdsR]); 
GIC_dHdsRMAD = mad([west_dHdsR; southeast_dHdsR; centraleast_dHdsR; northeast_dHdsR; north_dHdsR],1);
disp(['GIC dH/dSMB = ',num2str(round(GIC_dHdsmed,2)),' m/yr (R =',num2str(round(GIC_dHdsRmed,2)),')']);

