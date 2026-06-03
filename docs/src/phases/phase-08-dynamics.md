# Phase 8 — Multibody + robotics

**Status:** 🟢 In progress — MuJoCo adapter live for MJCF / URDF playback with constant-control policy and trajectory capture.

## Goal

Rigid-body multibody dynamics with contact joins the simulation
palette — for mechanisms, robot design, and as the dynamics backbone
for reinforcement-learning + control-systems flows.

## Capability inventory

- Rigid-body dynamics with Featherstone-algorithm integrators.
- Accurate contact with friction (convex Coulomb, Stewart-Trinkle).
- Joint library: revolute, prismatic, spherical, planar, 6-DOF,
  custom constraints.
- Actuators: PD / torque / velocity control, position servos.
- Sensors: IMU, force/torque, rangefinder, camera.
- URDF + MJCF import for existing robots (Franka, UR, Spot, etc.).
- Integration with trajectory-planning libraries (MoveIt, OMPL).

## Integrated tools graduating to Implemented

| Tool    | Adapter crate                  | Role                                     |
|---------|--------------------------------|------------------------------------------|
| MuJoCo  | `valenx-adapter-mujoco`        | Primary multibody + contact engine       |

## Acceptance checklist

- [ ] Load a URDF (Franka arm), visualise in the viewport.
- [ ] Simulate a pick-and-place task with contact.
- [ ] Record and replay trajectories; export to CSV / ROS bag.
- [ ] Control UI: slider-per-joint + a Python scripting hook.
- [ ] Couple with FEA for joint-stress estimation via preCICE.

## Success metrics

| Metric                                              | Target      |
|-----------------------------------------------------|-------------|
| 1 kHz simulation rate for a 7-DOF arm               | real-time   |
| Joint trajectory export round-trip losslessly       | yes         |

## Leads into

[Phase 9 — Multi-physics coupling](./phase-09-coupling.md).
