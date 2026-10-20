%% MECH5961 acoustic TMM design tool
% Edit only the STUDENT DESIGN PARAMETERS section. This single file contains
% the complete model and displays the predicted absorption response.

clear;
clc;
fprintf('Running: %s\n', mfilename('fullpath'));

%% ================= STUDENT DESIGN PARAMETERS =================
% Use one, two or three cell types. All vectors must have the same length.
% Lengths in this section are entered in millimetres for convenience.

typeName = ["Type 1", "Type 2", "Type 3"];
cellCount = [7 5 4];       % positive integers; total must be 16
cavityDepth_mm = [20 17.228999999999999 18.468];     % allowed range: 10 to 20 mm
neckRadius_mm = [0.73970000000000002 0.8266 1.0278]; % allowed range: 0.50 to 1.25 mm
neckExtension_mm = [0 0 0];   % allowed range: 0 to 10 mm and < cavity depth

% neckExtension_mm is the distance that the neck projects from the inner surface of
% the front panel into the cavity; it occupies part of the cavity depth and does not
% add to the total panel depth.

%% ================ COMMON COURSE SETTINGS ====================
% These settings are common to every group and must not be optimized.

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

%% Convert the student inputs to SI units and assemble the cells

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

%% Calculate the complete-panel response

frequency_Hz = makeFrequencyGrid(frequencyMin_Hz, frequencyMax_Hz, frequencyStep_Hz);
fprintf('Calculating response at %d frequencies...\n', numel(frequency_Hz));
response = panelResponse(cells, panelArea_m2, frequency_Hz, air, endCorrectionTerms);
bandMean = bandAverage(frequency_Hz, response.absorption, targetBand_Hz);

%% Display results

fprintf('Mean absorption, 800-1200 Hz: %.9f\n', bandMean);
fprintf('Constraint status: course limits passed.\n');

fprintf('Displaying absorption plot...\n');
fig = figure('Color', 'white', 'Name', 'MECH5961 Absorption Prediction', ...
    'NumberTitle', 'off', 'WindowStyle', 'normal', ...
    'Units', 'normalized', 'Position', [0.12, 0.12, 0.76, 0.76], ...
    'WindowState', 'normal', 'Visible', 'on');
ax = axes('Parent', fig, 'Visible', 'on');
hold(ax, 'on');
yl = [0, 1];
patch(ax, [targetBand_Hz(1), targetBand_Hz(2), targetBand_Hz(2), targetBand_Hz(1)], ...
      [yl(1), yl(1), yl(2), yl(2)], [0.69, 0.77, 0.87], ...
      'FaceAlpha', 0.45, 'EdgeColor', 'none', 'DisplayName', 'Target band');
plot(ax, frequency_Hz, response.absorption, 'Color', [0, 0, 0.5], ...
    'LineWidth', 1.8, 'DisplayName', sprintf('Band mean = %.4f', bandMean));
xlabel(ax, 'Frequency (Hz)');
ylabel(ax, 'Normal-incidence absorption coefficient');
xlim(ax, [frequency_Hz(1), frequency_Hz(end)]);
ylim(ax, yl);
grid(ax, 'on');
legend(ax, 'Location', 'best');
hold(ax, 'off');

% Select the figure and render it before the script returns.
set(fig, 'Visible', 'on');
figure(fig);
drawnow;
fprintf('Plot ready.\n');

%% ======================= MODEL FUNCTIONS =====================
% Students do not need to edit anything below this line.

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
