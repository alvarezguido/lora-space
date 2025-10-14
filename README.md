# lora-space
LoRa Direct-to-Satellite IoT Simulator

LoRa-Space is a simulator for Direct-to-Satellite IoT (DtS-IoT) networks based on LoRa technology. It is implemented in Python using the SimPy discrete-event simulation library.

The tool enables realistic modeling of LEO satellite-based LoRa networks, allowing users to evaluate performance under a variety of orbital, link, and traffic configurations.

## Scenarios
- narrow scenario: Optimized for dense device deployments within a limited coverage area.
- wider scenario: Suitable for large-scale or regional networks, simulating broader satellite footprints and more diverse device distributions.

## Features
- LoRa Network Simulation: Model and test fully customizable LoRa-based networks for data collection onboard LEO satellites.
- Orbital Configuration: Modify satellite orbital parameters such as altitude, inclination, and coverage footprint.
- Device and Traffic Modeling: Define the number of LoRa devices, their geographical distribution, and the traffic injection model (periodic, random, etc.).
- Parameter Customization: Configure key LoRa end-device parameters, including: Spreading Factor (SF), Transmission Power (TX Power), Data Payload Size, Bandwidth and Frequency Channel
- Performance Evaluation: Collect metrics such as packet delivery ratio, collision rate, and channel utilization.

## Requirements
- Python 3.x
- Libraries: simpy, numpy, math, sys

## Installation
> git clone https://github.com/alvarezguido/lora-space
- Once donwload, simple run the relevant script
- Edit the script parameters to adjust network size, traffic model, or orbital configuration.


This work was published on https://ieeexplore.ieee.org/document/9825690

![lora-space](https://github.com/alvarezguido/lora-space/assets/47746423/937fd5be-e15c-4f36-a779-5ee41bfb121e)
