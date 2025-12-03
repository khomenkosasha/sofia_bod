
# sofia_bod
Code for Sofia burden of disease (BoD) analysis within the UBD Policy project (https://ubdpolicy.eu/). This repository contains the R scripts used in the study *Health burden and inequities of urban environmental stressors in Sofia, Bulgaria* (2025). DOI: https://doi.org/10.1016/j.envres.2025.121782

## Overview

This project quantifies the mortality and morbidity burden attributable to major urban environmental stressors in Sofia: air pollution (PM2.5, NO2), insufficient green space, road-traffic noise, and the urban heat island (UHI) effect.

The analyses follow a comparative risk assessment framework, estimating health impacts by comparing current exposure levels with counterfactual “optimal” scenarios (e.g., WHO Air Quality Guidelines, recommended noise levels, sufficient access to green space, and no UHI effect).

The study also examines inequalities by socioeconomic status (SES) and uses spatial analyses (global and local Moran’s I) to identify clusters of high exposure, high burden, and vulnerability.

## Key Findings (from the published study)

* All residents live in areas exceeding WHO guidelines for PM2.5 and NO2.
* Green space, noise, and UHI exposures are suboptimal across the city.
* PM2.5 contributes the largest mortality burden, followed by NO2, noise, lack of green space, and UHI.
* Environmental exposures and associated burdens particularly affect lower-SES areas in the north, northeast, and northwest of Sofia.
* Results highlight environmental injustice and support the need for integrated planning addressing environmental and social inequalities.

## Repository Structure

* `bod_airp_montecarlo.R` – Monte Carlo health burden estimation for PM2.5/NO2.
* `bod_airp_montecarlo_SES.R` – Air pollution burden stratified by SES.
* `bod_noise_montecarlo.R` – Burden estimation for road-traffic noise.
* `bod_noise_montecarlo_SES.R` – Noise burden stratified by SES.
* `bod_gs_montecarlo.R` – Burden estimation for lack of green space.
* `bod_gs_montecarlo_SES.R` – Green space burden stratified by SES.
* `bod_uhi_montecarlo.R` – Burden estimation for urban heat island exposure.
* `bod_ses_analysis.R` – SES stratification, inequalities, and spatial (Moran’s I) analyses.

## Purpose

This repository supports:

* Baseline environmental health impact assessment for Sofia.
* Evaluation of inequalities and environmental justice.
* A transferable methodological framework for other cities.
* Evidence-based urban planning and health policy.

## License

This project is licensed under the **GPL-3.0** license.

## Citation

If using this code, please cite the study:

*Khomenko, S., Burov, A., Dzhambov, A. M., de Hoogh, K., Helbich, M., Mijling, B., Hlebarov, I., Popov, I., Dimitrova, D., Dimitrova, R., Markevych, I., Germanova, N., Brezov, D., Iungman, T., Montana, F., Chen, X., Gehring, U., Khreis, H., Mueller, N., Zapata-Diomedi, B., … Nieuwenhuijsen, M. (2025). Health burden and inequities of urban environmental stressors in Sofia, Bulgaria. Environmental research, 279(Pt 1), 121782. https://doi.org/10.1016/j.envres.2025.121782*

## Contact

For questions or collaboration, please contact the repository author.
