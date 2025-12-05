%% PS_plotSingleHSolPerturbation

%% time
simT = 0 : 0.001 : 0.6 ;
pertStyle = {'-', ':'} ;

%% plot states
figure ;

% CrossBridge attachment velocity v^S
subplot(4, 1, 1) ;
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.state.value(1 : 601, ip, 1), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[m/s]") ; title('CrossBridge Attachment Velocity', "Fontsize", 12) ;
xline(0.5) ; legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([-0.005 0.005]) ;

% CrossBridge attachment position l^S
subplot(4, 1, 2) ;
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.state.value(1 : 601, ip, 2), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[m]") ; title('CrossBridge Attachment Position', "Fontsize", 12) ;
xline(0.5) ; legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
ylim([0 0.01]) ;

% titin segment 1 length l^1
subplot(4, 1, 3) ;
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.state.value(1 : 601, ip, 3), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[m]") ; title('Titin Segment 1 Length', "Fontsize", 12) ;
xline(0.5) ; legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([0 0.01]) ;

% activation
subplot(4, 1, 4) ;
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.muscleDynamicsInfo.activation(1 : 601, ip), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[n.u.]") ; title('Activation', "Fontsize", 12) ;
xline(0.5) ; legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([]) ;

sgtitle("State Values") ;

%% plot state derivatives
figure ;

% CrossBridge attachment velocity v^S
subplot(4, 1, 1) ;
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.state.derivative(1 : 601, ip, 1), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[m/s/s]") ; title('CrossBridge Attachment Acceleration', "Fontsize", 12) ;
xline(0.5) ; legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([]) ;

% CrossBridge attachment position l^S
subplot(4, 1, 2) ;
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.state.derivative(1 : 601, ip, 2), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[m/s]") ; title('CrossBridge Attachment Velocity', "Fontsize", 12) ;
xline(0.5) ; legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([-0.005 0.005]) ;

% titin segment 1 length l^1
subplot(4, 1, 3) ;
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.state.derivative(1 : 601, ip, 3), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[m/s]") ; title('Titin Segment 1 Velocity', "Fontsize", 12) ;
xline(0.5) ; legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([]) ;

% activation
subplot(4, 1, 4) ;
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.muscleDynamicsInfo.activationDerivative(1 : 601, ip), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[1/s]") ; title('Activation Derivative', "Fontsize", 12) ;
xline(0.5) ; legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([]) ;

sgtitle("State Derivatives") ;

%% plot VEXAT component lengths
figure ;

% muscle length
subplot(3, 2, 1) ;
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.muscleLengthInfo.fiberLength(1 : 601, ip), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[m]") ; title('Fiber Length', "Fontsize", 12) ;
xline(0.5) ; % legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([]) ;

% path length l^P
subplot(3, 2, 2) ;
for ip = 1 : 2
    if ip == 1
        plot(simT, PS_inputState.pathState.lp_stretch(1 : 601), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
    else
        plot(simT, PS_inputState.pathState.lp_shorten(1 : 601), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
    end
end
xlabel("Time [s]") ; ylabel("[m]") ; title('Path Length', "Fontsize", 12) ;
xline(0.5) ; legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([]) ;

% CrossBridge length
subplot(3, 2, 3) ;
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.muscleLengthInfo.crossBridgeLength(1 : 601, ip), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[m]") ; title('CrossBridge Length', "Fontsize", 12) ;
xline(0.5) ; % legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([]) ;

% muscle length along tendon l^M * cos(alpha)
subplot(3, 2, 4) ;
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.muscleLengthInfo.fiberLengthAlongTendon(1 : 601, ip), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[m]") ; title('Fiber Length Along Tendon', "Fontsize", 12) ;
xline(0.5) ; legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([]) ;

% titin segment 2 length
subplot(3, 2, 5) ;
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.muscleLengthInfo.titin2Length(1 : 601, ip), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[m]") ; title('Titin Segment 2 Length', "Fontsize", 12) ;
xline(0.5) ; % legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([]) ;

% tendon length
subplot(3, 2, 6) ;
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.muscleLengthInfo.tendonLength(1 : 601, ip), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[m]") ; title('Tendon Length', "Fontsize", 12) ;
xline(0.5) ; legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([0.0401 - 1e-15 0.0401 + 1e-15]) ;

sgtitle("Muscle-Tendon Component Lengths") ;

%% plot VEXAT component velocities
figure ;

% muscle velo
subplot(3, 2, 1) ;
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.fiberVelocityInfo.fiberVelocity(1 : 601, ip), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[m/s]") ; title('Fiber Velocity', "Fontsize", 12) ;
xline(0.5) ; % legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([]) ;

% path velo v^P
subplot(3, 2, 2) ;
for ip = 1 : 2
    if ip == 1
        plot(simT, PS_inputState.pathState.dlp_stretch(1 : 601), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
    else
        plot(simT, PS_inputState.pathState.dlp_shorten(1 : 601), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
    end
end
xlabel("Time [s]") ; ylabel("[m/s]") ; title('Path Velocity', "Fontsize", 12) ;
xline(0.5) ; legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([]) ;

% CrossBridge velo
subplot(3, 2, 3) ;
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.muscleLengthInfo.crossBridgeVelocity(1 : 601, ip), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[m/s]") ; title('CrossBridge Velocity', "Fontsize", 12) ;
xline(0.5) ; % legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([]) ;

% muscle velo along tendon
subplot(3, 2, 4) ;
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.fiberVelocityInfo.fiberVelocityAlongTendon(1 : 601, ip), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[m/s]") ; title('Fiber Velocity Along Tendon', "Fontsize", 12) ;
xline(0.5) ; legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([]) ;

% titin segment 2 velo
subplot(3, 2, 5) ;
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.fiberVelocityInfo.titin2Velocity(1 : 601, ip), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[m/s]") ; title('Titin Segment 2 Velocity', "Fontsize", 12) ;
xline(0.5) ; % legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([]) ;

% tendon velo
subplot(3, 2, 6) ;
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.fiberVelocityInfo.tendonVelocity(1 : 601, ip), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[m/s]") ; title('Tendon Velocity', "Fontsize", 12) ;
xline(0.5) ; legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
ylim([-1e-15 1e-15]) ;

sgtitle("Muscle-Tendon Component Velocities") ;

%% plot VEXAT component forces
figure ;

% muscle force
subplot(3, 2, 1) ;
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.muscleDynamicsInfo.fiberForce(1 : 601, ip), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[N]") ; title('Fiber Force', "Fontsize", 12) ;
xline(0.5) ; % legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([0 1.1e-15]) ;

% fiber active force
subplot(3, 2, 3) ;
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.muscleDynamicsInfo.activeFiberForce(1 : 601, ip), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[N]") ; title('Fiber Active Force', "Fontsize", 12) ;
xline(0.5) ; % legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([0 1.1e-15]) ;

% fiber passive force
subplot(3, 2, 5) ;
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.muscleDynamicsInfo.passiveFiberForce(1 : 601, ip), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[N]") ; title('Fiber Passive Force', "Fontsize", 12) ;
xline(0.5) ; % legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([0 1.1e-15]) ;

% crossbridge force
subplot(3, 2, 2) ;
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.muscleDynamicsInfo.crossBridgeForce(1 : 601, ip), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[N]") ; title('CrossBridge Force', "Fontsize", 12) ;
xline(0.5) ; % legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([0 1.1e-15]) ;

% titin segment 1 force
subplot(3, 2, 4) ;
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.muscleDynamicsInfo.titin1Force(1 : 601, ip), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[N]") ; title('Titin Segment 1 Force', "Fontsize", 12) ;
xline(0.5) ; % legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([0 1.1e-15]) ;

% titin segment 2 force
subplot(3, 2, 6) ;
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.muscleDynamicsInfo.titin2Force(1 : 601, ip), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[N]") ; title('Titin Segment 2 Force', "Fontsize", 12) ;
xline(0.5) ; % legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([0 1.1e-15]) ;

sgtitle("Muscle-Tendon Component Forces") ;

figure;

% tendon force
for ip = 1 : 2
    plot(simT, mtInfoTimeSeries.muscleDynamicsInfo.tendonForce(1 : 601, ip), "b", "LineWidth", 2, "LineStyle", pertStyle{ip}) ; hold on ;
end
xlabel("Time [s]") ; ylabel("[N]") ; title('Tendon Force', "Fontsize", 12) ;
xline(0.5) ; legend({'Stretch Perturbation', 'Shorten Perturbation', 'Perturbation Onset'}, "Location", "eastoutside") ;
% ylim([0 1.1e-15]) ;

