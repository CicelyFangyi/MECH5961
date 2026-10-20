function run_submission
%RUN_SUBMISSION Recompute the assignment results and regenerate report figures.
% Base MATLAB only. The unchanged course functions are local to this file.
root=fileparts(mfilename('fullpath'));
oldPath=path; restorePath=onCleanup(@()path(oldPath)); addpath(root);
parameters=readtable(fullfile(root,'data','design_parameters.csv'),'TextType','string','Delimiter',',');
saved=load(fullfile(root,'data','figure_data.mat'),'figureData');
figureData=saved.figureData;
temporaryFolder=tempname(root);
assert(startsWith(string(temporaryFolder),string([root filesep])));
mkdir(temporaryFolder);
removeTemporary=onCleanup(@()cleanup_temporary(temporaryFolder));
summaries=cell(3,1);
names=["baseline","best_two_type","final"];
for i=1:3
    name=names(i); p=parameters(parameters.design_id==name,:);
    geometry_mm=[p.cavity_depth_mm,p.neck_radius_mm,p.neck_extension_mm];
    out=course_tmm(geometry_mm,p.cell_count);
    [peak,index]=max(out.response.absorption); % first maximum chooses lowest frequency
    summaries{i}=table(name,out.J,out.frequency_Hz(index),peak,'VariableNames', ...
        {'design_id','J','peak_frequency_Hz','peak_absorption'});
    if name=="best_two_type", continue; end
    if name=="baseline"
        assert(isequal(geometry_mm,figureData.parametric.baseline.geometry_mm), ...
            'The saved parameter-study figures use a different baseline.');
        script=fullfile(temporaryFolder,'MECH5961_TMM_design.m');
        write_variant(fullfile(root,'MECH5961_TMM_design.m'),script,geometry_mm,p.cell_count);
    else
        script=fullfile(root,'MECH5961_TMM_design.m');
    end
    official=execute_script(script); % isolates the teacher's clear statement
    assert(isequal(official.geometry_mm,geometry_mm) && isequal(official.counts(:),p.cell_count(:)), ...
        'The parameter table and teacher script disagree.');
    differences=[max(abs(official.response.absorption-out.response.absorption)), ...
        max(abs(official.response.normalized_admittance-out.response.normalized_admittance)), ...
        max(abs(official.response.reflection-out.response.reflection)),abs(official.J-out.J)];
    assert(all(differences<=1e-9),'Teacher script and reused model disagree.');
    folder=fullfile(root,name);
    if ~isfolder(folder), mkdir(folder); end
    absorption=table(official.frequency_Hz,official.response.absorption, ...
        'VariableNames',{'frequency_Hz','absorption'});
    writetable(absorption,fullfile(folder,'absorption.csv'));
    % Export the Figure returned by the actual teacher script.
    set(official.fig,'Visible','off');
    axesHandles=findall(official.fig,'Type','axes');
    for k=1:numel(axesHandles), axesHandles(k).Toolbar.Visible='off'; end
    exportgraphics(official.fig,fullfile(folder,'absorption.png'),'Resolution',600);
    close(official.fig);
    figureData.(name)=plot_data(out);
end
writetable(vertcat(summaries{:}),fullfile(root,'data','performance.csv'));
final=figureData.final;
refined=course_tmm(final.geometry_mm,final.counts,(200:.25:2500)',80);
verification=table([1;.25],[40;80],[final.J;refined.J],[0;refined.J-final.J], ...
    'VariableNames',{'frequency_step_Hz','end_correction_terms','J','delta_J'});
writetable(verification,fullfile(root,'data','numerical_verification.csv'));
% Sixteen-identical-cell comparisons are separate counterfactual panels.
figureData.counterfactual_absorption=zeros(numel(final.frequency_Hz),size(final.geometry_mm,1));
for k=1:size(final.geometry_mm,1)
    identical=course_tmm(final.geometry_mm(k,:),16);
    figureData.counterfactual_absorption(:,k)=identical.response.absorption;
end
save(fullfile(root,'data','figure_data.mat'),'figureData','-v7');
make_report_figures;
fprintf('Final J = %.12f; refined J = %.12f; change = %+.12g\n',final.J,refined.J,refined.J-final.J);
end

function data=plot_data(out)
% Retain only arrays needed to draw the report figures.
data.geometry_mm=out.geometry_mm;
data.counts=out.counts;
data.frequency_Hz=out.frequency_Hz;
data.J=out.J;
data.response.absorption=out.response.absorption;
data.response.normalized_admittance=out.response.normalized_admittance;
data.Y_type=out.Y_type;
data.eta=out.eta;
end

function out=execute_script(scriptFile)
% The teacher script clears this small workspace, not the outer runner.
run(scriptFile);
out.frequency_Hz=frequency_Hz;
out.response=response;
out.J=bandMean;
out.fig=fig;
out.geometry_mm=[cavityDepth_mm(:),neckRadius_mm(:),neckExtension_mm(:)];
out.counts=cellCount(:);
end

function write_variant(sourceFile,destination,geometry_mm,counts)
% Change only the five permitted assignments for the baseline reference run.
text=fileread(sourceFile);
first=strfind(text,'%% ================= STUDENT DESIGN PARAMETERS');
last=strfind(text,'%% ================ COMMON COURSE SETTINGS');
assert(isscalar(first) && isscalar(last) && first<last);
student=text(first:last-1);
names='"Type 1"';
for k=2:size(geometry_mm,1), names=[names ', "Type ' num2str(k) '"']; end %#ok<AGROW>
fields={'typeName','cellCount','cavityDepth_mm','neckRadius_mm','neckExtension_mm'};
values={['[' names ']'],mat2str(counts(:)',17),mat2str(geometry_mm(:,1)',17), ...
    mat2str(geometry_mm(:,2)',17),mat2str(geometry_mm(:,3)',17)};
for k=1:5
    pattern=['(?m)^' fields{k} '[ \t]*=[^;]*;'];
    assert(numel(regexp(student,pattern))==1);
    student=regexprep(student,pattern,[fields{k} ' = ' values{k} ';'],'once');
end
fid=fopen(destination,'w','n','UTF-8'); assert(fid>0); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',[text(1:first-1) student text(last:end)]);
end

function cleanup_temporary(folder)
% Remove only the known temporary baseline file and its now-empty directory.
file=fullfile(folder,'MECH5961_TMM_design.m');
if isfile(file), delete(file); end
if isfolder(folder), rmdir(folder); end
end

% The course_tmm entry and all nine original local model functions follow.

function out = course_tmm(geometry_mm, counts, frequency, terms, cellsOnly)
%COURSE_TMM Research entry to the supplied MECH5961 model [C].
% geometry_mm: K-by-3 [cavity depth, neck radius, neck extension] in mm.
% counts: positive integer counts summing to 16. Frequency is a column in Hz.
% terms: 40 for official evaluation; 80 only for independent verification.
% cellsOnly skips the panel calculation when caching single-cell admittances.
% All nine local model functions below are copied verbatim from [C].
if nargin < 2 || isempty(counts), counts = 16; end
if nargin < 3, frequency = []; end
if nargin < 4 || isempty(terms), terms = 40; end
if nargin < 5, cellsOnly = false; end
validateattributes(geometry_mm, {'double'}, {'real','finite','2d','ncols',3});
K = size(geometry_mm,1);
validateattributes(counts, {'double'}, {'real','finite','vector','numel',K});
validateattributes(terms, {'double'}, {'scalar','integer','positive','finite'});
if ~isempty(frequency)
    validateattributes(frequency, {'double'}, {'real','finite','positive','column'});
    assert(numel(frequency)>1 && all(diff(frequency)>0), 'Frequency must increase.');
end
typeName = "Type " + string((1:K)');
cellCount = counts;
cavityDepth_mm = geometry_mm(:,1);
neckRadius_mm = geometry_mm(:,2);
neckExtension_mm = geometry_mm(:,3);
panelArea_m2          = 0.0025;       % 50 mm x 50 mm
cellFootprintArea_m2  = 0.00015625;   % 12.5 mm x 12.5 mm
cavityWidth_m         = 0.011;
cavityHeight_m        = 0.011;
panelThickness_m      = 0.0015;
backingThickness_m    = 0.0015;
extendedNeckWall_m    = 0.0006;

frequencyMin_Hz       = 200.0;
frequencyMax_Hz       = 2500.0;
frequencyStep_Hz      = 1.0;
endCorrectionTerms    = 40;
targetBand_Hz         = [800.0, 1200.0];

% Common air properties, matching the Python implementation.
air.density_kg_m3              = 1.204;
air.pressure_pa                = 101325.0;
air.gamma                      = 1.4;
air.dynamic_viscosity_pa_s     = 1.81e-5;
air.heat_capacity_j_kg_k       = 1005.0;
air.thermal_conductivity_w_m_k = 0.0257;

typeName = string(typeName(:));
cellCount = cellCount(:);
cavityDepth_m = cavityDepth_mm(:) * 1e-3;
neckRadius_m = neckRadius_mm(:) * 1e-3;
neckExtension_m = neckExtension_mm(:) * 1e-3;

nTypes = numel(typeName);
if any([numel(cellCount), numel(cavityDepth_m), numel(neckRadius_m), ...
        numel(neckExtension_m)] ~= nTypes)
    error('All student design vectors must have the same number of entries.');
end

cells = repmat(struct(), nTypes, 1);
for k = 1:nTypes
    cells(k).name = typeName(k);
    cells(k).count = cellCount(k);
    cells(k).footprint_area_m2 = cellFootprintArea_m2;
    cells(k).cavity_width_m = cavityWidth_m;
    cells(k).cavity_height_m = cavityHeight_m;
    cells(k).cavity_depth_m = cavityDepth_m(k);
    cells(k).panel_thickness_m = panelThickness_m;
    cells(k).backing_thickness_m = backingThickness_m;
    cells(k).neck_radius_m = neckRadius_m(k);
    cells(k).neck_extension_m = neckExtension_m(k);
    if neckExtension_m(k) > 0
        cells(k).neck_wall_thickness_m = extendedNeckWall_m;
    else
        cells(k).neck_wall_thickness_m = 0;
    end
end

validateCourseDesign(cells, panelArea_m2);

if isempty(frequency)
    frequency = makeFrequencyGrid(frequencyMin_Hz, frequencyMax_Hz, frequencyStep_Hz);
end
out.frequency_Hz = frequency;
out.geometry_mm = geometry_mm;
out.counts = counts(:);
out.cells = cells;
out.endCorrectionTerms = terms;
out.air = air;
out.panelArea_m2 = panelArea_m2;
if ~cellsOnly
    out.response = panelResponse(cells, panelArea_m2, frequency, air, terms);
    out.J = bandAverage(frequency, out.response.absorption, targetBand_Hz);
end
% p/U impedances are already converted to z using internal cavity area.
% A single cell contributes (A_c/A_panel)/z, NOT (1/16)/z.
out.z = complex(zeros(numel(frequency),K));
for k = 1:K
    out.z(:,k) = cellImpedance(cells(k),frequency,air,terms);
end
out.y = (cavityWidth_m*cavityHeight_m/panelArea_m2)./out.z;
if ~cellsOnly
    out.Y_type = out.y .* counts(:)';
    out.eta = 4*real(out.Y_type)./abs(1+out.response.normalized_admittance).^2;
    out.lengths_m = zeros(K,3); % front correction, cavity correction, effective neck
    for k = 1:K
        front = endCorrection(cells(k).neck_radius_m,sqrt(cellFootprintArea_m2/pi),terms);
        cavity = endCorrection(cells(k).neck_radius_m,sqrt(cavityWidth_m*cavityHeight_m/pi),terms);
        out.lengths_m(k,:) = [front,cavity,panelThickness_m+cells(k).neck_extension_m+front+cavity];
    end
end
end

function f = makeFrequencyGrid(lo, hi, step)
    if ~(isfinite(lo) && isfinite(hi) && isfinite(step) && lo > 0 && hi > lo && step > 0)
        error('Frequency bounds and step must be finite, positive and ordered.');
    end
    f = (lo:step:(hi - eps(hi))).';
    if isempty(f) || f(end) < hi
        f(end+1, 1) = hi;
    else
        f(end) = hi;
    end
end

function validateCourseDesign(cells, panelArea)
    if isempty(cells) || numel(cells) > 3
        error('Use between one and three cell types.');
    end
    names = string({cells.name});
    if numel(unique(names)) ~= numel(names)
        error('Cell type names must be unique.');
    end
    totalCount = 0;
    occupiedArea = 0;
    for k = 1:numel(cells)
        c = cells(k);
        if ~(isscalar(c.count) && isfinite(c.count) && c.count >= 1 && c.count == fix(c.count))
            error('%s: count must be a positive integer.', c.name);
        end
        if c.cavity_depth_m < 0.010 || c.cavity_depth_m > 0.020
            error('%s: cavity depth must be between 10 and 20 mm.', c.name);
        end
        if c.neck_radius_m < 0.00050 || c.neck_radius_m > 0.00125
            error('%s: neck radius must be between 0.50 and 1.25 mm.', c.name);
        end
        if c.neck_extension_m < 0 || c.neck_extension_m > 0.010 || ...
                c.neck_extension_m >= c.cavity_depth_m
            error('%s: neck extension must be 0-10 mm and shorter than the cavity depth.', c.name);
        end
        if c.neck_extension_m > 0 && abs(c.neck_wall_thickness_m - 0.0006) > 1e-12
            error('%s: an extended neck must use the common 0.60 mm wall thickness.', c.name);
        end
        cavityArea = c.cavity_width_m * c.cavity_height_m;
        if cavityArea > c.footprint_area_m2 * (1 + 1e-12)
            error('%s: internal cavity area exceeds its footprint.', c.name);
        end
        outerRadius = c.neck_radius_m + (c.neck_extension_m > 0) * c.neck_wall_thickness_m;
        if 2 * outerRadius >= min(c.cavity_width_m, c.cavity_height_m)
            error('%s: the centered neck does not fit inside the cavity.', c.name);
        end
        totalCount = totalCount + c.count;
        occupiedArea = occupiedArea + c.count * c.footprint_area_m2;
    end
    if totalCount ~= 16
        error('The cell counts must sum to 16.');
    end
    if occupiedArea > panelArea * (1 + 1e-12)
        error('The cell footprints exceed the panel area.');
    end
end

function response = panelResponse(cells, panelArea, f, air, endTerms)
    if isempty(cells) || ~(isscalar(panelArea) && isfinite(panelArea) && panelArea > 0)
        error('Provide at least one cell type and a positive panel area.');
    end
    admittance = complex(zeros(size(f)));
    for k = 1:numel(cells)
        c = cells(k);
        cavityArea = c.cavity_width_m * c.cavity_height_m;
        weight = c.count * cavityArea / panelArea;
        admittance = admittance + weight ./ cellImpedance(c, f, air, endTerms);
    end
    reflection = (1 - admittance) ./ (1 + admittance);
    alpha = 1 - abs(reflection).^2;
    if any(~isfinite(alpha)) || any(alpha < -1e-9) || any(alpha > 1 + 1e-9)
        error('Absorption lies outside the passive range; review the design.');
    end
    response.frequencies_hz = f;
    response.normalized_admittance = admittance;
    response.reflection = reflection;
    response.absorption = min(max(alpha, 0), 1);
end

function normalized = cellImpedance(cellData, f, air, endTerms)
    omega = 2 * pi * f;
    area = cellData.cavity_width_m * cellData.cavity_height_m;
    perimeter = 2 * (cellData.cavity_width_m + cellData.cavity_height_m);
    viscousRadius = 2 * area / perimeter;
    areaRadius = sqrt(area / pi);
    frontRadius = sqrt(cellData.footprint_area_m2 / pi);

    effectiveLength = cellData.panel_thickness_m + cellData.neck_extension_m;
    effectiveLength = effectiveLength + endCorrection(cellData.neck_radius_m, frontRadius, endTerms);
    effectiveLength = effectiveLength + endCorrection(cellData.neck_radius_m, areaRadius, endTerms);

    [zNeck, kNeck] = effectiveMedium(omega, cellData.neck_radius_m, ...
        pi * cellData.neck_radius_m^2, air);
    [zCavity, kCavity] = effectiveMedium(omega, viscousRadius, area, air);

    lowerDepth = cellData.cavity_depth_m - cellData.neck_extension_m;
    loadAdmittance = 1i * tan(kCavity * lowerDepth) ./ zCavity;
    if cellData.neck_extension_m > 0
        outerNeckArea = pi * (cellData.neck_radius_m + cellData.neck_wall_thickness_m)^2;
        surroundingArea = area - outerNeckArea;
        zSurrounding = zCavity * area / surroundingArea;
        loadAdmittance = loadAdmittance + ...
            1i * tan(kCavity * cellData.neck_extension_m) ./ zSurrounding;
    end

    cosineTerm = cos(kNeck * effectiveLength);
    sineTerm = sin(kNeck * effectiveLength);
    inputImpedance = (cosineTerm + 1i * zNeck .* sineTerm .* loadAdmittance) ./ ...
        (1i * sineTerm ./ zNeck + cosineTerm .* loadAdmittance);
    soundSpeed = sqrt(air.gamma * air.pressure_pa / air.density_kg_m3);
    normalized = inputImpedance * area / (air.density_kg_m3 * soundSpeed);
    if any(~isfinite(normalized)) || any(real(normalized) < -1e-9)
        error('Nonfinite or nonpassive cell impedance; review geometry and resolution.');
    end
end

function [zCharacteristic, waveNumber] = effectiveMedium(omega, radius, area, air)
    nu = air.dynamic_viscosity_pa_s / air.density_kg_m3;
    thermalDiffusivity = air.thermal_conductivity_w_m_k / ...
        (air.density_kg_m3 * air.heat_capacity_j_kg_k);
    effectiveDensity = air.density_kg_m3 ./ circularLossFunction(omega, nu, radius);
    compressibility = (air.gamma - (air.gamma - 1) .* ...
        circularLossFunction(omega, thermalDiffusivity, radius)) / ...
        (air.gamma * air.pressure_pa);
    zCharacteristic = sqrt(effectiveDensity ./ compressibility) / area;
    waveNumber = omega .* sqrt(effectiveDensity .* compressibility);
end

function F = circularLossFunction(omega, diffusivity, radius)
    q = radius * sqrt(-1i * omega / diffusivity);
    small = abs(q) < 1e-3;
    safeQ = q;
    safeQ(small) = 1 + 0i;
    % MATLAB's scale=1 applies the same exponential factor to J0 and J1;
    % it therefore cancels in their ratio, matching scipy.special.jve.
    exact = 1 - 2 * besselj(1, safeQ, 1) ./ (safeQ .* besselj(0, safeQ, 1));
    series = -q.^2 / 8 - q.^4 / 48 - 11 * q.^6 / 3072;
    F = exact;
    F(small) = series(small);
end

function correction = endCorrection(neckRadius, receivingRadius, terms)
    roots = j1Roots(terms);
    x = roots * neckRadius / receivingRadius;
    correction = 4 * neckRadius * sum(besselj(1, x).^2 ./ ...
        (x .* (roots .* besselj(0, roots)).^2));
end

function roots = j1Roots(count)
    persistent cachedRoots
    if isempty(cachedRoots) || numel(cachedRoots) < count
        rootsToFind = max(count, 80);
        cachedRoots = zeros(rootsToFind, 1);
        for m = 1:rootsToFind
            cachedRoots(m) = fzero(@(x) besselj(1, x), [m*pi, (m + 0.5)*pi]);
        end
    end
    roots = cachedRoots(1:count);
end

function value = bandAverage(f, alpha, band)
    low = band(1);
    high = band(2);
    if low >= high || low < f(1) || high > f(end)
        error('The frequency grid must cover the complete target band.');
    end
    interior = f > low & f < high;
    bandF = [low; f(interior); high];
    bandAlpha = interp1(f, alpha, bandF, 'linear');
    value = trapz(bandF, bandAlpha) / (high - low);
end
