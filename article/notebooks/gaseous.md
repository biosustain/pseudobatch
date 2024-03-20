# Description of gaseous species handling
Gaseous species are either not or only partially affected by sample withdrawal, therefore they require special treatment when using the pseudo batch transformation. To derive the special treatment for gaseous species we start at a general mass balance over a bioreactor for any one species.
$$
\frac{dM_{species}(t)}{dt}=In_{gas,\: species}(t)+In_{liquid,\: species}(t)-Out_{gas,\: species}(t)+Metabolism_{species}(t)-Sampled_{species}(t)
$$
where $M_{species}$ is the total amount of species in the reactor in both gas and aqueous phase. The term $Metabolism_{species}$ is the net result of the microbial metabolism and covers both production (positive) and consumption (negative). The $In$ variables describe the amount of species entering the bioreactor through inlet stream of gas or liquid. The $Out_{gas,\: species}$ describes the amount of species leaving the bioreactor through outlet gas stream. We only consider fed-batch operations, thus liquid outlet stream is omitted. Finally, the term $Sampled_{species}$ describes the amount of species removed through sampling. All the the terms are functions over time, $t$. If we integrate the mass balance we get the following:

$$
M_{species}(t)=\int_0^t In_{gas,\: species}(t) dt +\int_0^t In_{liquid,\: species}(t) dt - \int_0^t Out_{gas,\: species}(t) dt + \int_0^t Metabolism_{species}(t) dt - \int_0^t Sampled_{species}(t) dt 
$$

This shows that the concentration of the species at a give time is result of the integrated mass balance, i.e. the accumulated "events", e.g. input of medium, gas, sampling etc. We are interested in the microbes contribution to the mass balance thus we will isolate the metabolism term.

$$
\int_0^t Metabolism_{species}(t) dt = M_{species}(t) - \int_0^t In_{gas,\: species}(t) dt +\int_0^t In_{liquid,\: species}(t) dt + \int_0^t Out_{gas,\: species}(t) dt + \int_0^t Sampled_{species}(t) dt 
$$

In a real world setting we do not have access to the continuous integrals. Instead we have measurements of the which estimates the terms. The liquid concentration, $M_{species}(t)$ is commonly measured through light spectrometry, HPLC or mass spectrometry methods. The inlet integrals are estimated through flow rates and known concentrations in inlet gas and liquid. The outlet integral is estimated using an off-gas measuring device and the mass removal due to sampling can be calculated the from the sample volume and the liquid concentration measurement. However all these are discrete measurements which estimates each integral at a specific time point. For the rest of this section we will use measured integral values which we will denote as **Accumulated quantities**. To indicate that these contain discrete values we will use t in square brackets, $[t]$, to denote the value at a specific time point. In this discrete notation the mass balance look as follows.

$$
AccumMetabolism_{species}[t] = M_{species}[t] - AccumIn_{gas,\: species}[t] + AccumIn_{liquid,\: species}[t] + AccumOut_{gas,\: species}[t] + AccumSampled_{species}[t] 
$$

The pseudo batch Python package provide a utility function which calculate the accumulated metabolism using the mass balance above. The function is called `metabolised_amount`.

## Hypothetical concentration
Before applying the pseudo batch transformation to the measurements of the gaseous species, the measurements needs some preprocessing. We say that we calculate the hypothetical concentration of the species. This is the concentration is the concentration that we would have observed if the species did not evaporate. Hypothetical concentration have to take into account that the if some the species did not evaporate more of the species would be removed through sample withdrawal. We calculate the hypothetical concentration using the following equation.

$$
Hypothetical \: concentration_{species}[t] = (AccumMetabolism[t] - AccumSampleLoss[t]) * (V[t] - Sample\:volume[t])^{-1}
$$

Where 
$$
AccumMetabolism[t] = \int_0^t Metabolism_{species}(t) dt \\
AccumSampleLoss[t] = AccumSampleLoss[t-1] + \sum_1^t(AccumMetabolism[t] - AccumSampleLoss[t-1]) * \frac{Sample\:volume[t]}{V[t]}
$$

The calculation of the hypothetical concentration is provide through in the Python function `hypothetical_concentration`. The hypothetical concentrations can be transform through the pseudo batch transformation. To see examples of how to this please see the [tutorial](https://biosustain.github.io/pseudobatch/Tutorials/8%20-%20Special%20case%20gaseous%20species.html) in our documentation.
