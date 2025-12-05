# Description
This repository is a modified version of https://github.com/mjhmilla/Millard2023VexatMuscle.

Modifications include the addition of several MATLAB scripts and data files necessary to reproduce our simulations.

In our simulations, the muscle (the Human Soleus model you created) holds a certain level of activation against a fixed load (isometric conditions), before being either stretched or compressed by 1.6 mm (within 40 ms), from its optimal length. In our first runs, we simulate a muscle with a rigid tendon. Prior to perturbation, muscle length is constant over time, while muscle activation is non-zero.

# Simulation details
Muscle activation begins at zero and saturates approaching 100%. Equations of activation and derivative:
$$
a = 1 - e^{-12t}
$$
$$
\dot{a} = -12e^{-12t}
$$

From 0 to 500 ms, path length = (tendon slack length + muscle optimal length * cos(pennation angle)). From 500 ms to 540 ms, the path length either increases or decreases by 1.6 mm. From 540 ms to 1000 ms, the path length remains constant at initial +/- 1.6 mm. Equations of path length and velocity during the perturbation (500 to 540 ms):
$$
l^{P} = P(3\psi^{2}(t-0.5)^{2} - 2\psi^{3}(t-0.5)^{3}) + l^{P}(t=0)
$$
$$
v^{P} = 6P(\psi^{2}(t-0.5) - \psi^{3}(t-0.5)^{2})
$$
$$
P = \pm 1.6 mm
$$
$$
\psi = 25
$$
$$
l^{P}(t=0) = l^{T}_{s} + l^{M}_{0}cos(\alpha_{0})
$$

# How to run the simulations
Execute the MATLAB script 'main_PS_OuterLoop.m' to run the simulations.

To choose whether to run the simulations with either an elastic or rigid tendon, set the variable 'select_useElasticTendon' equal to either 1 (Elastic) or 0 (Rigid).

Running the simulations will produce figures of various model outputs and the variable 'mtInfoTimeSeries' which contains the model outputs. Model inputs are stored in the variable 'PS_inputState'. 'mtInfoTimeSeries' is a struct with the same fields as in the 'mtInfo' variables used throughout the original repo, however each field contains the full time series of that value rather than its value at a single instant.

# Modifications to repository
- Modifications to Millard's original code
    - Changed 'rootDirName' in 'getRootProjectDirectory.m' to 'Millard2023VexatMuscle_PetersSergiCustomPert'
- Added MATLAB scripts (.m)
    - main_PS_OuterLoop.m
        - Location: main folder
        - Purpose: top layer of simulation. analagous to 'main_OuterLoop.m' from the original repo.
    - main_PS_CreateModels_OuterLoop
        - Location: subfolder 'PetersSergiCustomPert'
        - Purpose: creates muscle model for simulation. analagous to 'main_CreateModels_OuterLoop' from the original repo. the human soleus model here is identical to the one in the original repo.
    - main_PS_SingleAntagonistPert_HumanSoleus_OuterLoop
        - Location: subfolder 'PetersSergiCustomPert'
        - Purpose: analagous to 'main_KirschBoskovRymer1994_OuterLoop' from the original repo.
    - main_PS_SingleAntagonistPert_HumanSoleus
        - Location: subfolder 'PetersSergiCustomPert'
        - Purpose: analagous to 'main_KirschBoskovRymer1994' from the original repo.
    - PS_runSingleAgonistPert_HumanSoleus
        - Location: subfolder 'PetersSergiCustomPert'
        - Purpose: calls model wrapper function 'calcMillard2023VexatMuscleInfo.m' to run simulation. somewhat analagous to 'runKirschBoskovRymer1994SimulationsVexat.m' from the original repo, but with much more deviation from that file than the previous analagous pairs
    - PS_plotSingleHSolPerturbation
        - Location: subfolder 'PetersSergiCustomPert'
        - Purpose: generate plots of model outputs.
- Added MATLAB data files (.mat)
    - PS_inputState
        - Location: subfolder 'PetersSergiCustomPert'
        - Purpose: contains the time series of model inputs (path length, path velocity, activation, activation derivative)
        - Overview:
            - .activationState
                - .a : muscle activation
                - .da : derivative of muscle activation
            - .pathState
                - .lp_stretch : path length in the perturbation where the MTU is lengthened
                - .lp_shorten : path length in the perturbation where the MTU is shortened
                - .dlp_stretch : path velocity in the perturbation where the MTU is lengthened
                - .dlp_shorten : path velocity in the perturbation where the MTU is shortened
    - mtInfoTimeSeries_ElasticTendon
        - Location: subfolder 'PetersSergiCustomPert'
        - Purpose: contains the time series of model outputs generated on my machine for both perturbations (stretch and shorten) with elastic tendon
        - Overview: a struct with the same fields as in the 'mtInfo' variables used throughout the original repo, however each field contains the full time series of that value rather than its value at a single instant
    - mtInfoTimeSeries_RigidTendon
        - Location: subfolder 'PetersSergiCustomPert'
        - Purpose: contains the time series of model outputs generated on my machine for both perturbations (stretch and shorten) with rigid tendon
        - Overview: a struct with the same fields as in the 'mtInfo' variables used throughout the original repo, however each field contains the full time series of that value rather than its value at a single instant

