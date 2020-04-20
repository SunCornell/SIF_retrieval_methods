%%%%%%%%%%%%%%%%%%%%%%%
% Main and helper functions for SIF retrieval methods used in the manuscript
% "Systematic assessment of retrieval methods for canopy far-red solar-induced chlorophyll 
% fluorescence (SIF) using automated high-frequency field spectroscopy"
%
% Authors/Contributors: 
% Chang, Christine Y. (For questions, contact: cyc54@cornell.edu)
% Frankenberg, Christian
% Guanter, Luis
% Kornfeld, Ari
% Gu, Lianhong
% Sun, Ying
%
%%%%%%%%%%%%%%%%%%%%%%%

% Setup 
function [run, models] = run_setup(dataset, sfmprange, doasrange, sfmp, fld3, doas, svda)
    % set run settings
    run = [];
    models = [];
    run.fracDOY = dataset.fracDOY; %decimal Julian day, contained in dataset struct
    run.wl = dataset.wl; %wavelength associated with each pixel of QE Pro. this should exclude the first and last 10 pixels which are optically inactive
    run.Ei_ret = dataset.rad_sky; %incoming irradiance (E), contained in dataset struct. this is size (nmeas, pixels) where pixels excludes optically inactive pixels
    run.Eo_ret = dataset.rad_veg; %outgoing irradiance (L), see above
    run.nmeas = size(run.Ei_ret,1); %number of measurements
    run.hf = dataset.hf; %prescribed shape of fluorescence emission for DOAS and SVD. this is size (pixels,1)

    %which models to run? 0 = off; 1 = on
    models.sfmp = sfmp;
    models.fld3 = fld3;
    models.doas = doas;
    models.svd = svda;

	%fitting window range for SFM and DOAS
    run.sfmp_wlrange = sfmprange;  %e.g. [759, 767.76]
    run.doas_wlrange = doasrange;  %e.g. [745, 759]
end

% Example SVD settings
% if rolling > 0, abs(rolling) = the number of spectra in a centered moving window used for training the SVD
% if rolling = 0, it will use all spectra from the day for training the SVD
% if rolling < 0, it will use that number of spectra but evenly dispersed throughout the day for training the SVD

svd_settings = [];
svd_settings.wl_range = [759.5 761.5]; %narrow O2A fitting window
svd_settings.rolling = 5; %moving window of 5 spectra 

%%%%%%%%%%%%%%%%%%%%%%%
% SIF retrieval wrapper
function [result] = retrieve_SIF(run, models, svd_settings)

    result = [];

    %3FLD
    if models.fld3 == 1
        disp("Starting FLD3")
        tic;
        [SIF, ~, ~, ~, result.fld3_refl] = FLD3(run);
        result.fld3_fs0 = SIF/1000/pi; %apply 1000/pi to bi-hemispherical SIF only
        disp(toc/60);
    end
    
    result = daily_svdfit(run, models, svd_settings, result);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIF retrieval methods
function result = daily_svdfit(run, models, svd_settings, result)
    disp("Running daily retrievals for: SFM, DOAS, SVD")
    doy_list = unique(floor(run.fracDOY));

    %SFM: assume polynomial refl and fluo, 2nd order poly
    if models.sfmp == 1
        tic;
        idxr = find(run.wl > min(run.sfmp_wlrange) & run.wl < max(run.sfmp_wlrange));
		idx760 = run.wl(floor(run.wl)==760);
		idx760 = find(idxr == idx760(1) ); 
        SIF = zeros(run.nmeas,1);
        disp("Starting SFM ")
        for z=1:length(doy_list)
            doy = doy_list(z);
            msg = "[" + string(run.sfmp_wlrange(1))+"-"+string(run.sfmp_wlrange(2))+"]: DOY "+string(doy);
            disp(msg)
            idx_doy = find(floor(run.fracDOY) == doy);
            [sif, SFM] = sfm_poly2(run, idx760, idx_doy);
            SIF(idx_doy) = sif/1000/pi; %apply 1000/pi to bi-hemispherical SIF only
        end
        result.sfmp_fs0 = SIF;
        result.sfm = SFM; %meta
        disp(toc/60);
    end
    

	%
    if models.doas == 1
        tic;
        disp("Starting DOAS ")
        SIF = zeros(run.nmeas,1);
        for z=1:length(doy_list)
            doy = doy_list(z);
            msg = "[" + string(run.doas_wlrange(1))+"-"+string(run.doas_wlrange(2))+"]: DOY "+string(doy);
            disp(msg)
            idx_doy = find(floor(run.fracDOY) == doy);
            [FR, ~] = doasFit(run, idx_doy);
            SIF(idx_doy) = FR.SIF/1000/pi; %apply 1000/pi to bi-hemispherical SIF only
            refl(idx_doy,:) = FR.refl;
        end
        result.doas_fs0 = SIF;
        result.doas = FR; %meta
        result.doas_refl = refl;
        disp(toc/60);
    end
    
    
    if models.svd == 1
        disp("Starting SVD ")
        tic;
        SIF = zeros(run.nmeas,1) * nan;
        idxr  = find(run.wl > min(svd_settings.wl_range) & run.wl < max(svd_settings.wl_range)); 
        refl = zeros(run.nmeas, length(run.wl));
        for z=1:length(doy_list)
            doy = doy_list(z);
            msg = "["+string(svd_settings.wl_range(1))+"-"+string(svd_settings.wl_range(2))+"]: DOY "+string(doy);
            disp(msg)
            if svd_settings.rolling == 0
                svdtype = "all";
            else
                svdtype = "local";
            end
            idx_doy = find(floor(run.fracDOY) == doy);
            if ((svd_settings.rolling > 0) && (length(idx_doy) < svd_settings.rolling))
				% not enough measurements to perform the SVD! NAN
                SVD = [];
                SVD.svd_fs0 = ones(length(idx_doy),1) * nan;
                SVD.refl = ones(length(idx_doy),1024) * nan;
            else
                [SVD] = rolling_svd_(run.Ei_ret(idx_doy,:), run.Eo_ret(idx_doy,:), run.wl, svd_settings, run.hf); 
			end
            SIF(idx_doy) = SVD.svd_fs0/1000/pi; %apply 1000/pi to bi-hemispherical SIF only
            refl(idx_doy,:) = SVD.refl; 
            
        end
        result.svd.svd_fs0 = SIF;
        result.svd.svd_refl = refl;
        disp(toc/60);
    end
end

function [fluo, SFM] = sfm_poly2(run, idx760, idx_doy)
% Spectral Fitting Method from Meroni et al. 2009
    wavelength = run.wl;
    idxr = find(wavelength > min(run.sfmp_wlrange) & wavelength < max(run.sfmp_wlrange));
    nmeas_ = length(idx_doy);
    sky_spec = run.Ei_ret(idx_doy,:);
    sample_spec = run.Eo_ret(idx_doy,:);
   
    SFM = [];
    SFM.refl = zeros(nmeas_,length(idxr));
    SFM.fluo = zeros(nmeas_,1);
    SFM.SIF = zeros(nmeas_,1);
    SFM.reduced_chisq = zeros(nmeas_,1);
    SFM.rmse = zeros(nmeas_,1);
    SFM.rrmse = zeros(nmeas_,1);
    SFM.spec_meas = zeros(nmeas_,length(idxr))*NaN;
    SFM.spec_mod = zeros(nmeas_,length(idxr))*NaN;
    for i = 1:nmeas_
        dlambda = wavelength(idxr) - wavelength(idxr(1));
        p1 = dlambda.^2 .* sky_spec(i, idxr)';
        p2 = dlambda .* sky_spec(i, idxr)';
        p3 = sky_spec(i, idxr)';
        p4 = dlambda.^2;
        p5 = dlambda;
        p6 = ones(size(idxr));
        K = [p1,p2,p3,p4,p5,p6];
        b = K\sample_spec(i,idxr)';
        y_hat = K*b;
        y = sample_spec(i,idxr)';
        rmse(i) = sqrt(mean((y_hat - y).^2));
        rrmse(i) = rmse(i) ./ mean(y);
        fluo_ = b(4).*dlambda.^2 + b(5).*dlambda + b(6);
        refl(i,:) = (sample_spec(i,idxr) - fluo_') ./ sky_spec(i,idxr);
        fluo(i) = fluo_(idx760);
        nparams = size(K,2);
        sv_end = b(end);
        [~, SFM.reduced_chisq(i)] = reduced_chisq_screen(y', y_hat', length(idxr), nparams, sv_end, 1);
        [SFM.rmse(i), SFM.rrmse(i), ~] = error_calc(y', y_hat');
    
        SFM.spec_meas(i,:) = y';
        SFM.spec_mod(i,:) = y_hat';
    end
    [fluo] = quality_filter(fluo, SFM.reduced_chisq, SFM.rrmse, 0);
    SFM.refl = refl;
    SFM.wavewindow = wavelength(idxr);
end  

function [fluo, refl, Eidxo, Eidxi] = sFLD(run) 
%Standard Fraunhofer Line Discrimination method (sFLD) from Plascyk and Gabriel 1975, Damm et al. 2011
    nmeas_ = run.nmeas; 
    sky_spec = run.Ei_ret;
    sample_spec = run.Eo_ret;
    sky_wl = run.wl;
    sample_wl = run.wl;
    outer = [758, 759];
    inner = [760, 761];
    Eidxo  = sky_wl > outer(1) & sky_wl < outer(2);
    Eidxi  = sky_wl > inner(1) & sky_wl < inner(2); 
    Lidxo  = sample_wl > outer(1) & sample_wl < outer(2);
    Lidxi  = sample_wl > inner(1) & sample_wl < inner(2); 
    
    fluo = zeros(nmeas_,1);
    refl = zeros(nmeas_,1);
    for i = 1:nmeas_
        sky_in = mean(sky_spec(i,Eidxi));
        sky_out = mean(sky_spec(i,Eidxo));
        sample_in = mean(sample_spec(i,Lidxi));
        sample_out = mean(sample_spec(i,Lidxo));
        fluo(i) = (sky_out.*sample_in - sky_in.*sample_out)./(sky_out - sky_in); 
        refl(i) = ((sample_out - sample_in)./(sky_out - sky_in));
    end
end

function [fluo, EidxL, Eidxi, EidxR, refl] = FLD3(run) 
%Modified Fraunhofer Line Discrimination method (3FLD) from Maier et al. 2003, Damm et al. 2011
    nmeas_ = run.nmeas; 
    sky_spec = run.Ei_ret;
    sample_spec = run.Eo_ret;
    sky_wl = run.wl;
    sample_wl = run.wl;
    left = [758, 759];
    inner = [760, 761];
    right = [770, 775];
    EidxL  = sky_wl > left(1) & sky_wl < left(2); %outside / left shoulder
    Eidxi  = sky_wl > inner(1) & sky_wl < inner(2); %inside
    EidxR = sky_wl > right(1) & sky_wl < right(2); % right shoulder
    LidxL  = sample_wl > left(1) & sample_wl < left(2); %outside / left shoulder
    Lidxi  = sample_wl > inner(1) & sample_wl < inner(2); %inside
    LidxR = sample_wl > right(1) & sample_wl < right(2); % right shoulder

    fluo = zeros(nmeas_,1);
    refl = zeros(nmeas_,1);
    for i = 1:nmeas_
        sky_in = min(sky_spec(i,Eidxi));
        sky_left = mean(sky_spec(i,EidxL));
        sky_right = mean(sky_spec(i,EidxR));
        sample_in = min(sample_spec(i,Lidxi));
        sample_left = mean(sample_spec(i,LidxL));
        sample_right = mean(sample_spec(i,LidxR));
        wL = (mean(sky_wl(EidxR))-mean(sky_wl(Eidxi)))./(mean(sky_wl(EidxR))-mean(sky_wl(EidxL)));
        wR = (mean(sky_wl(Eidxi))-mean(sky_wl(EidxL)))./(mean(sky_wl(EidxR))-mean(sky_wl(EidxL)));
        fluo(i) = (sample_in - (sky_in./((wL.*sky_left) + (wR.*sky_right))) .* ((wL.*sample_left) + (wR.*sample_right))) ./ (1-(sky_in ./ ((wL.*sky_left) + (wR.*sky_right))));
        refl(i) = (sample_spec(i,494)-fluo(i)) / sky_spec(i,494);
    end
end

function [FR,xx] = doasFit(run, idx_doy) 
% Differential Optical Absorption Spectroscopy method (Grossmann et al. 2018) modified from Christian Frankenberg's scripts
    idxr = find(run.wl > min(run.doas_wlrange) & run.wl < max(run.doas_wlrange));
    nmeas_ = length(idx_doy);
    wavelength = run.wl;
    sky_spec = run.Ei_ret(idx_doy, :); 
    sample_spec = run.Eo_ret(idx_doy,:);
    hf = run.hf;
	xx = ((wavelength(idxr))-mean(wavelength(idxr)))/100; 

	idx_760 = (wavelength >= 760 & wavelength <= 760.05);
    hf_meanctr = (hf/mean(hf));
    h = hf_meanctr(idxr);
    FL760 = hf_meanctr(idx_760);

    xl = xx/max(abs(xx)); %normalize wavelengths into -1 to 1 range
    
	nn = length(idxr);
    le = nmeas_;
    FR.SIF = zeros(le,1)*NaN;
    FR.reduced_chisq = zeros(le,1)*NaN;
    FR.rmse = zeros(le,1)*NaN;
    FR.rrmse = zeros(le,1)*NaN;
    FR.SIF_1sigma = zeros(le,1)*NaN;
    FR.res = zeros(le,nn)*NaN;
    FR.spec_meas = zeros(le,nn)*NaN;
    FR.spec_mod = zeros(le,nn)*NaN;

    K = [];
    minBIC = 1000000000;
	legendre_num = 6; 
    
    for j = 0:legendre_num 
        K = [K legendreP(j,xl)];
    end

    FR.refl = zeros(le, length(wavelength))*NaN;

    for i=1:nmeas_

        diffSVD = [K h./sample_spec(i,idxr)']; %Philipp Kohler's Legendre polynomial solution

        y = log(sample_spec(i,idxr))-log(sky_spec(i,idxr));
        x = diffSVD\y';

        mod = diffSVD*x;
        resid = mod-y';

		Ssquared = (mod-y').^2; 
		FR.res(i,:) = resid; 
		FR.spec_meas(i,:) = y;
		FR.spec_mod(i,:) = mod;
		FR.wavewindow = wavelength(idxr);
		S = sum(Ssquared)/(length(y)-length(x))*inv(transpose(diffSVD)*diffSVD);
		FR.SIF_1sigma(i) = sqrt(S(end,end));
		FR.idxr = idxr;
		nparams = size(diffSVD,2);
		sv_end = x(end);
		[FR.SIF(i), FR.reduced_chisq(i)] = reduced_chisq_screen(y, mod', length(idxr), nparams, sv_end, FL760);
		SIF_range = sv_end .* hf_meanctr;
		FR.refl(i,:) = (sample_spec(i,:)-SIF_range') ./ sky_spec(i,:);
		[~, FR.rrmse(i), ~] = error_calc(y, mod');
    end
    
    [FR.SIF] = quality_filter(FR.SIF, FR.reduced_chisq, FR.rrmse, 1);
end

function [svd_] = rolling_svd_(spec_train, spec_veg, wavelengths, svd_settings, hf)
% Singular vector decomposition method, inspired by code from Ari Kornfeld and Christian Frankenberg
% if rolling > 0, abs(rolling) = the number of spectra in a centered moving window used for training the SVD
% if rolling = 0, it will use all spectra from the day for training the SVD
% if rolling < 0, it will use that number of spectra but evenly dispersed throughout the day for training the SVD
% Also includes automated optimization, see below. 
    
    % *** Step 1: define default parameters ***
    rolling = svd_settings.rolling;
    
    pixels = find(wavelengths <= max(svd_settings.wl_range) & wavelengths >= min(svd_settings.wl_range));
    
    %irrad coefs at 760nm
    FL760idx = wavelengths >= 760 & wavelengths <= 760.05; 
    
    %prescribed SIF shape (Guanter et al. 2013 RSE)
    FL760 = hf(FL760idx);
    hf_ = hf(pixels);
    
    spectra = spec_train(:,pixels);
    spec_sky = spectra;

    SIF2 = 0;
    refl = 0;
    FL760idx_ = wavelengths(pixels) <= 760.05 & wavelengths(pixels) >= 760;
    
    np1_combos = linspace(3,6,4);  
    np2_combos = linspace(3,6,4);
    sv_combos = combvec(1, np1_combos, np2_combos);
   
    np1_ = 0;
    np2_ = 0;
    numSV_ = 0;
    minBIC = 1000000000; %start with arbitrarily large number...
    svd_ = [];
    svd_.numSV = [];
    svd_.EV = [];
    svd_.rmse = [];
    svd_.rrmse = [];
    svd_.rmsp = [];
    
    if rolling > 1 %will use middle of rolling window as spec veg
        [n, ~] = size(spectra);
        spec_meas = zeros(n-rolling+1, length(pixels));
        spec_mod = zeros(n-rolling+1, length(pixels));
        sv_end = zeros(n-rolling+1, 1);
        
        for j=1:length(sv_combos)
            for k = 1:(n-rolling+1)
                spec_sky_window = spec_sky(k:(k+rolling-1),:);
                [K_frag, ~, FLidx, numSV, EV, all_V, norm_S, nparams] = train_svd_(spec_sky_window, wavelengths, pixels, hf_, sv_combos(2,j), sv_combos(3,j), norm_);
                [~, spec_mod(k,:), spec_meas(k,:), sv_frag] = solve_svd_(K_frag, spec_veg(k+floor(rolling/2), :), pixels, wavelengths);
                sv_end(k) = sv_frag(end);
            end
            [minBIC, svd_] = optimize_svd(spec_mod, spec_meas, pixels, FL760, minBIC, np1_, np2_, numSV_, sv_combos(2,j), sv_combos(3,j), numSV, rolling, EV, norm_S, svd_, nparams, sv_end, spec_train(1:(n-rolling+1),:), spec_veg(1:(n-rolling+1),:), hf, wavelengths, all_V);
        end

    elseif rolling < 0 %use x time points as training spectra
        n = abs(rolling);
        midx = (1:n)*floor(size(spectra,1)/n); %25 scattered throughout data
        for j=1:length(sv_combos)
            start_train = tic;
            [K, ~, FLidx, numSV, EV, all_V, norm_S, nparams] = train_svd_(spec_sky(midx,:), wavelengths, pixels, hf_, sv_combos(2,j), sv_combos(3,j), norm_);
            [~, spec_mod, spec_meas, sv] = solve_svd_(K, spec_veg, pixels, wavelengths);
            sv_end = sv(:,end);
            [minBIC, svd_] = optimize_svd(spec_mod, spec_meas, pixels, FL760, minBIC, np1_, np2_, numSV_, sv_combos(2,j), sv_combos(3,j), numSV, rolling, EV, norm_S, svd_, nparams, sv_end, spec_train, spec_veg, hf, wavelengths, all_V);
        end
    else

        for j = 1:length(sv_combos)
            [K, ~, FLidx, numSV, EV, all_V, norm_S, nparams] = train_svd_(spec_sky, wavelengths, pixels, hf_, sv_combos(2,j), sv_combos(3,j), norm_);
            [~, spec_mod, spec_meas, sv] = solve_svd_(K, spec_veg, pixels, wavelengths);
            sv_end = sv(:,end);
            [minBIC, svd_] = optimize_svd(spec_mod, spec_meas, pixels, FL760, minBIC, np1_, np2_, numSV_, sv_combos(2,j), sv_combos(3,j), numSV, rolling, EV, norm_S, svd_, nparams, sv_end, spec_train, spec_veg, hf, wavelengths, all_V);
        end 

    end
    
end
    
function [K, Ksd, FLidx, numSV, EV, all_V, norm_S, nparams] = train_svd_(training_spec, wavelengths, pixels, hf, np1, np2, norm_)
    % *** Step 2: train SVD ***
    [U,S,V] = svd(training_spec, 'econ');
    numSV = cutoff_threshold(training_spec, S);
    norm_S = diag(S)/sum(diag(S))*100;
    
    % Return at least the first 2 Eigenvectors passing the eigenvalue threshold :
    evToreturn = 1:min( max(2, numSV), size(V, 2)); % at least 2
    EV = -V(:,evToreturn);  % Singular vectors always come out -ve, so reverse them here. 
    all_V = -V;
    
    % Compute Matrices for the SV convolved with polynomial terms in lambda space which is needed to avoid spurious correlation of the power terms
    poly1 = zeros(size(EV,1),np1);
    poly2 = zeros(size(EV,1),np2);

    % 1st SV
    idx = 1;
    for i=np1:-1:0 %0
        poly1(:,idx) =  EV(:, 1).*(wavelengths(pixels) - mean(wavelengths(pixels))).^i;
        idx = idx + 1;
    end

    if size(EV,2) > 1
        % 2nd SV
        idx = 1;
        for i=np2:-1:0 %0
            poly2(:,idx) =  EV(:, 2).*(wavelengths(pixels) - mean(wavelengths(pixels))).^i;
            idx = idx + 1;
        end

        % remaining SV
        EV3plus = EV(: , 3:numSV);
        K = [poly1, poly2, EV3plus, hf];
    else
        K = [poly1, hf];
    end
    nparams = size(K,2);
    
    FLidx = (size(K, 2) - size(hf, 2)+1):size(K,2);
    Ksd = ones(1, size(K, 2));
    if norm_ == 1
        %Ari's version: normalize results so SD = 1 but do not center
        Ksd = std(K); %to disable, use: ones(1, size(R.K, 2));
        Kmean = zeros(1,size(K, 2)); %mean(K); 
        K = MxV(@minus, K, Kmean);
        K = MxV(@rdivide, K, Ksd);
    end
end

function [result, spec_mod, spec_meas, sv] = solve_svd_(K, spec_veg, pixels, wavelengths)
% ** Step 3: solve SVD *********************************************
    spec_meas = spec_veg(:,pixels);
    result = K\spec_meas';
    spec_mod = (K*result)'; %modeled spectrum
    sv = result';
end

function [minBIC, svd_] = optimize_svd(spec_mod, spec_meas, pixels, FL760, minBIC, np1_, np2_, numSV_, np1, np2, numSV, rolling, EV, norm_S, svd_, nparams, sv_end, spec_sky, spec_veg, hf, wavelengths, all_V)
% %% Step 4: optimize SVD using BIC ***************************
    norm_spec_mod = zeros(size(spec_mod));
    norm_spec_meas = zeros(size(spec_meas));
    save_time = 0;
    SIF = 0;
    rmse = 0;
    rrmse = 0;
    rmsp = 0;
    refl = 0;
    rolling_ = "see prev";
    
    for i=1:size(spec_mod,1)
        norm_spec_mod(i,:) = (spec_mod(i,:) - min(spec_mod(i,:))) / (max(spec_mod(i,:) - min(spec_mod(i,:))));
        norm_spec_meas(i,:) = (spec_meas(i,:) - min(spec_meas(i,:))) / (max(spec_meas(i,:) - min(spec_meas(i,:))));
    end
    
    BIC = zeros(1,size(spec_mod,1));
    AIC = zeros(1,size(spec_mod,1));
    for i=1:size(spec_mod,1)
        resids = norm_spec_meas(i,:)-norm_spec_mod(i,:);
        mu = nanmean(resids);
        sigma = nanstd(resids); 
        L = (1/(sqrt(2*pi)*sigma))^size(resids,2) * exp(-sum(((resids-mu).^2)/2*sigma.^2));
        LL = -(size(resids,2))*log(sqrt(2*pi))-(size(resids,2))*log(sigma) - sum(((resids-mu).^2)/(2*sigma.^2));
        BIC(i) = nparams*log(length(pixels)) - (2.*(LL));
        AIC(i) = 2*nparams - 2.*(LL);
    end
    SIF_range = zeros(length(sv_end), length(hf));
    if nanmean(BIC, 'all') < minBIC
        start_save = tic;
        [SIF, reduced_chisq] = reduced_chisq_screen(spec_meas, spec_mod, size(spec_meas,2), nparams, sv_end, FL760);
        [rmse, rrmse, rmsp] = error_calc(spec_meas, spec_mod);
        [SIF2] = quality_filter(SIF, reduced_chisq, rrmse, 2);
        minBIC = nanmean(BIC, 'all');
        np1_ = np1;
        np2_ = np2;
        numSV_ = numSV;
        rolling_ = rolling;
        save_time = toc(start_save);
        SIF2 = [zeros(floor(rolling/2),1) * NaN; SIF2; zeros(floor(rolling/2),1) * NaN]; %kludge to insert NaNs if using rollingsvd
        for i = 1:length(sv_end)
            SIF_range(i,:) = sv_end(i) .* hf;
        end
        svd_.refl = (spec_veg-SIF_range) ./ spec_sky; 
        svd_.refl = [zeros(floor(rolling/2),size(svd_.refl,2))*NaN; svd_.refl; zeros(floor(rolling/2),size(svd_.refl,2))*NaN];
        svd_.spec_meas = spec_meas;
        svd_.spec_mod = spec_mod;
        svd_.wavewindow = wavelengths(pixels);
        svd_.svd_fs0 = SIF2;
        svd_.numSV = [svd_.numSV, numSV_];
        svd_.EV = EV;
        svd_.all_V = all_V;
        svd_.norm_S = norm_S;
        svd_.reduced_chisq = reduced_chisq;
        svd_.rmse = [svd_.rmse, rmse];
        svd_.rrmse = [svd_.rrmse, rrmse];
        svd_.rmsp = [svd_.rmsp, rmsp];
        svd_.np1 = np1_;
        svd_.np2 = np2_;
    end
    msg = "BIC: "+string(minBIC)+", rolling = "+rolling_+"; numSV = ["+string(min(svd_.numSV))+"-"+string(min(svd_.numSV))+"]; np1 = "+string(np1_)+"; np2 = "+string(np2_);
    disp(msg);
end



%%
% Support functions
function numSV = cutoff_threshold(training_spec, s)
    %  HARD THRESHOLD FROM GAVISH & DONOHO 2014. use s from [u,s,v] for input.
    beta = min(size(training_spec)) / max(size(training_spec));
    if min(size(training_spec)) <= 25
        omega = (0.56*beta^3) - (0.95*beta^2) + (1.82*beta) + 1.43; % approximation, Gavish & Donoho 2014; this can save some time if it's really slow. only use this for the local window
    else
        omega = optimal_SVHT_coef(beta, 0);
    end
    cutoff = median(diag(s)) * omega;
    idx_keep = find(s(s > cutoff));
    numSV = max(idx_keep);
end

function [SIF, reduced_chisq] = reduced_chisq_screen(spec_meas, spec_mod, npx, nparams, sv_end, FL760)
% calculate reduced chi squared
    degf = npx - nparams; 
    reduced_chisq = zeros(size(spec_meas,1),1)*NaN;
    SIF = sv_end * FL760;
    for i = 1:size(spec_meas,1)
        var_resid = var(spec_meas(i,:) - spec_mod(i,:), 0, 2);   
        chisq = nansum((spec_meas(i,:)- spec_mod(i,:)).^2 ./ var_resid); 
        reduced_chisq(i) = chisq/degf;
    end
end

function [SIF] = quality_filter(SIF, reduced_chisq)
% Check for goodness of fit using reduced chi squared (Guanter et al. 2012, Sun et al. 2018)
    chisq_max = 2;
    chisq_min = 0.5;
    for i=1:size(SIF,1)
        if (reduced_chisq(i) > chisq_max) || (reduced_chisq(i) < chisq_min) 
            SIF(i) = NaN;
        end
    end
end

function [rmse, rrmse, rmsp] = error_calc(spec_meas, spec_mod)
    rmse = sqrt(nanmean((spec_mod - spec_meas).^2, 2));
    rrmse = (rmse ./ abs(nanmean(spec_meas,2)) .* 100);
    rmsp = rmse ./ spec_meas .* 100;
end